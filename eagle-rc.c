/*
EAGLE: explicit alternative genome likelihood evaluator
Utility program that classifies reads based on EAGLE calculated likelihoods (from verbose output)

Run EAGLE with --rc which sets: --omega=1e-40 --mvh --verbose --isc
  where main variant evaluation output goes to STDIN
  while per read likelihoods (verbose) go to STDERR

Other program options are situational

ex) eagle -t 2 -v var.vcf -a align.bam -r ref.fa --rc 1> out.txt 2> readinfo.txt

Copyright 2016 Tony Kuo
This program is distributed under the terms of the GNU General Public License
*/

#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <getopt.h>

#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/khash.h"
#include "util.h"
#include "calc.h"
#include "vector.h"

/* Constants */
#define ALPHA 1.3     // Factor to account for longer read lengths lowering the probability a sequence matching an outside paralogous source
#define LOGALPHA (log(ALPHA))

/* Command line arguments */
static int debug;
static int listonly;
static int readlist;
static int reclassify;
static int refonly;
static int paired;
static int pao;
static int isc;
static int nodup;
static int splice;
static int ngi;
static int phred64;
static int bisulfite;
static int const_qual;
static double omega, lgomega;

/* Time info */
static time_t now; 
static struct tm *time_info; 
#define print_status(M, ...) time(&now); time_info = localtime(&now); fprintf(stderr, M, ##__VA_ARGS__);

KHASH_MAP_INIT_STR(rh, vector_t) // hashmap: string key, vector value
static khash_t(rh) *read_hash; // pointer to hashmap, read

KHASH_MAP_INIT_STR(orh, vector_t) // hashmap: string key, vector value
static khash_t(orh) *other_read_hash; // pointer to hashmap, reads from other_bam files

KHASH_MAP_INIT_STR(vh, vector_t) // hashmap: string key, vector value
static khash_t(vh) *var_hash;  // pointer to hashmap, variant

KHASH_MAP_INIT_STR(rsh, vector_t)   // hashmap: string key, vector value
static khash_t(rsh) *refseq_hash; // pointer to hashmap

static void add2var_list(vector_t *var_list, char *set) {
    char var[strlen(set)];

    int n;
    char *s;
    size_t i;
    for (s = set + 1; sscanf(s, "%[^;];%n", var, &n) == 1; s += n) { // scan variant set
        char *v = strdup(var);
        for (i = 0; i < var_list->len; i++) {
            if (strcmp(v, (char *)var_list->data[i]) == 0) break;
        }
        if (i == var_list->len) { vector_add(var_list, v); }
        else { free(v); v = NULL; } // if already exists then skip
        if (*(s + n) == ']') break;
    }
}

static int var_read(const char* filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) { exit_err("failed to open file %s\n", filename); }

    int nvars = 0;
    char *line = NULL;
    ssize_t read_file = 0;
    size_t line_length = 0;
    while ((read_file = getline(&line, &line_length, file)) != -1) {
        if (line_length <= 0 || line[strspn(line, " \t\v\r\n")] == '\0') continue; // blank line
        if (line[0] == '#') continue;

        int pos, n;
        char chr[line_length], ref[line_length], alt[line_length], set[line_length];

        int t = sscanf(line, "%s %d %s %s %*[^\t] %*[^\t] %*[^\t] %*[^\t] %*[^\t] %s", chr, &pos, ref, alt, set);
        if (t < 4) { exit_err("bad fields in EAGLE output file\n%s\n", line); }

        if (strcmp(set, "[]") == 0) snprintf(set, line_length, "[%s,%d,%s,%s;]", chr, pos, ref, alt);

        char *s;
        for (s = set + 1; sscanf(s, "%[^;];%n", chr, &n) == 1; s += n) { // scan variant set and add to hash
            char *v = strdup(chr);

            vector_t *node;
            khiter_t k = kh_get(vh, var_hash, chr);
            if (k != kh_end(var_hash)) {
                node = &kh_val(var_hash, k);
            }
            else {
                int absent;
                k = kh_put(vh, var_hash, v, &absent);
                node = &kh_val(var_hash, k);
                if (absent) vector_init(node, 8, VOID_T);
            }
            vector_add(node, v);

            nvars++;
            if (*(s + n) == ']') break;
        }

    }
    free(line); line = NULL;
    fclose(file);
    return nvars;
}

static int readinfo_read(const char* filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) { exit_err("failed to open file %s\n", filename); }

    int nreads = 0;
    char *line = NULL;
    ssize_t read_file = 0;
    size_t line_length = 0;
    while ((read_file = getline(&line, &line_length, file)) != -1) {
        if (line_length <= 0 || line[strspn(line, " \t\v\r\n")] == '\0') continue; // blank line
        if (line[0] == '#') continue;

        int pos;
        double prgu, prgv, pout;
        char name[line_length], chr[line_length], flag[line_length], var[line_length];
        int t = sscanf(line, "%s %s %d %lf %lf %lf %*[^\t] %*[^\t] %s %s", name, chr, &pos, &prgu, &prgv, &pout, flag, var);
        if (t < 8) {
            t = sscanf(line, "%s %s %d %lf %lf %lf %*[^\t] %*[^\t] %s", name, chr, &pos, &prgu, &prgv, &pout, var);
            if (t < 7) { exit_err("bad fields in EAGLE read info file\n%s\n", line); }
            flag[0] = '\0';
        }

        if (prgu == prgv) continue; // if ref and alt probabilities are equal, it didn't "align" and was in a splice zone

        int is_read2 = 0;
        int n;
        char *s, token[strlen(flag) + 1];
        for (s = flag; sscanf(s, "%[^,]%n", token, &n) == 1; s += n + 1) {
            if (strcmp("READ2", token) == 0) is_read2 = 1;
            if (*(s + n) != ',') break;
        }

        char key[256];
        snprintf(key, 256, "%s\t%d", name, is_read2);

        khiter_t k = kh_get(rh, read_hash, key);
        if (k != kh_end(read_hash)) {
            vector_t *node = &kh_val(read_hash, k);
            read_t **r = (read_t **)node->data;
            size_t i;
            for (i = 0; i < node->len; i++) {
                if (strcmp(r[i]->name, name) == 0 && strcmp(r[i]->qseq, key) == 0) {
                    r[i]->prgu = log_add_exp(r[i]->prgu, prgu);
                    r[i]->prgv = log_add_exp(r[i]->prgv, prgv);
                    add2var_list(r[i]->var_list, var);
                    break;
                }
            }
            if (i == node->len) { exit_err("failed to find %s in hash key %d\n", key, k); }
        }
        else {
            read_t *r = read_create(name, 0, chr, pos);
            r->prgu = prgu;
            r->prgv = prgv;
            r->pout = pout;
            r->flag = strdup(flag);
            r->qseq = strdup(key);
            add2var_list(r->var_list, var);

            int absent;
            k = kh_put(rh, read_hash, r->qseq, &absent);
            vector_t *node = &kh_val(read_hash, k);
            if (absent) vector_init(node, 8, READ_T);
            vector_add(node, r);
            nreads++;
        }
    }
    free(line); line = NULL;
    fclose(file);
    return nreads;
}

static void readinfo_classify() {
    /* read->index set to classification:
    reference = 0
    alternative = 1
    reference-like-allele, multi-allelic, closer to the reference = 2;
    multi-allelic that are undifferentiateable = 3
    unknown, ambiguous with equal likelihoods for reference and alternative = 4
    */

    size_t readi, i, j;
	khiter_t k;
    for (k = kh_begin(read_hash); k != kh_end(read_hash); k++) {
		if (kh_exist(read_hash, k)) {
            vector_t *node = &kh_val(read_hash, k);
            read_t **r = (read_t **)node->data;
            for (readi = 0; readi < node->len; readi++) {
                char **v = (char **)r[readi]->var_list->data;
                size_t nvariants = r[readi]->var_list->len;

                int multiallele;
                if (nvariants > 1) { // check if only multi-allelic variants at the same position & no hope of differentiating ref vs alt
                    multiallele = 1;
                    int prev_pos = -1;
                    for (i = 0; i < nvariants; i++) {
                        int pos;
                        int t = sscanf(v[i], "%*[^,],%d,%*[^,],%*[^,],", &pos);
                        if (t < 1) { exit_err("bad fields in %s\n", v[i]); }
                        if (prev_pos != -1 && pos - prev_pos != 0) {
                            multiallele = 0;
                            break;
                        }
                        prev_pos = pos;
                    }
                    if (multiallele) {
                        r[readi]->index = 3;
                        fprintf(stdout, "%s\tMUL\t%s\t%d\t%f\t%f\t%f\t%s\t", r[readi]->name, r[readi]->chr, r[readi]->pos, r[readi]->prgu, r[readi]->prgv, r[readi]->pout, r[readi]->flag);
                        for (i = 0; i < nvariants; i++) { fprintf(stdout, "%s;", v[i]); } fprintf(stdout, "\n");
                        continue;
                    }
                }

                // if any variant in list is not in EAGLE output, it suggests variants are from a different phase if alt wins
                multiallele = 0;
                for (i = 0; i < nvariants; i++) {
                    khiter_t k = kh_get(vh, var_hash, v[i]);
                    if (k != kh_end(var_hash)) {
                        vector_t *var_hash_node = &kh_val(var_hash, k);
                        char **vv = (char **)var_hash_node->data;
                        for (j = 0; j < var_hash_node->len; j++) {
                            if (strcmp(v[i], vv[j]) == 0) break;
                        }
                        if (j == var_hash_node->len) multiallele = 1;
                    }
                    else {
                        multiallele = 1;
                    }
                }

                fprintf(stdout, "%s\t", r[readi]->name);
                if (r[readi]->prgu > r[readi]->prgv && r[readi]->prgu - r[readi]->prgv >= 0.69) { // ref wins
                    r[readi]->index = 0;
                    fprintf(stdout, "REF\t");
                }
                else if (r[readi]->prgv > r[readi]->prgu && r[readi]->prgv - r[readi]->prgu >= 0.69) { // alt wins
                    if (multiallele) { // EAGLE outputs the set with highest likelihood ratio, i.e. most different from reference, leaving the "reference-like-allele"
                        r[readi]->index = 2;
                        fprintf(stdout, "RLA\t");
                    }
                    else {
                        r[readi]->index = 1;
                        fprintf(stdout, "ALT\t");
                    }
                }
                else { // unknown
                    r[readi]->index = 4;
                    fprintf(stdout, "UNK\t");
                }
                fprintf(stdout, "%s\t%d\t%f\t%f\t%f\t%s\t", r[readi]->chr, r[readi]->pos, r[readi]->prgu, r[readi]->prgv, r[readi]->pout, r[readi]->flag);
                for (i = 0; i < nvariants; i++) { fprintf(stdout, "%s;", v[i]); } fprintf(stdout, "\n");
            }
        }
    }
    fflush(stdout);
    print_status("# Reads Classified:\t%s", asctime(time_info));
}

static void bam_write(const char *bam_file, const char *output_prefix, char *other_bam, int reverse) {
    size_t i;
    other_read_hash = kh_init(orh);
    if (other_bam != NULL) {
        int n;
        char *f, fn[strlen(other_bam) + 1];
        for (f = other_bam; sscanf(f, "%[^,]%n", fn, &n) == 1; f += n + 1) {
            samFile *sam_in = sam_open(fn, "r"); // open bam file
            if (sam_in == NULL) { exit_err("failed to open BAM file %s\n", fn); }
            bam_hdr_t *bam_header = sam_hdr_read(sam_in); // bam header
            if (bam_header == 0) { exit_err("bad header %s\n", fn); }

            bam1_t *aln = bam_init1(); // initialize an alignment
            while (sam_read1(sam_in, bam_header, aln) >= 0) {
                if (aln->core.tid < 0) continue; // not mapped
                int is_secondary = 0;
                char *flag = bam_flag2str(aln->core.flag);

                int n2;
                char *s, token[strlen(flag) + 1];
                for (s = flag; sscanf(s, "%[^,]%n", token, &n2) == 1; s += n2 + 1) {
                    if (strcmp("SECONDARY", token) == 0 || strcmp("SUPPLEMENTARY", token) == 0) is_secondary = 1;
                    if (*(s + n2) != ',') break;
                }
                free(flag); flag = NULL;
                if (pao && is_secondary) continue;

                int found = 0;
                char *name = strdup((char *)aln->data);

                vector_t *node;
                khiter_t k = kh_get(orh, other_read_hash, name);
                if (k != kh_end(other_read_hash)) {
                    node = &kh_val(other_read_hash, k);
                }
                else {
                    int absent;
                    k = kh_put(orh, other_read_hash, name, &absent);
                    node = &kh_val(other_read_hash, k);
                    if (absent) vector_init(node, 8, VOID_T);
                }
                for (i = 0; i < node->len; i++) {
                    if (strcmp(node->data[i], name) == 0) {
                        found = 1;
                        break;
                    }
                }
                if (!found) vector_add(node, name);
            }
            bam_destroy1(aln);
            bam_hdr_destroy(bam_header);
            sam_close(sam_in);
            print_status("# Read other bam:\t%s\t%s", fn, asctime(time_info));
            if (*(f + n) != ',') break;
        }
    }

    samFile *sam_in = sam_open(bam_file, "r"); // open bam file
    if (sam_in == NULL) { exit_err("failed to open BAM file %s\n", bam_file); }
    bam_hdr_t *bam_header = sam_hdr_read(sam_in); // bam header
    if (bam_header == 0) { exit_err("bad header %s\n", bam_file); }

    /* Split input bam files into respective categories' output bam files */
    samFile *out;
    char out_fn[256];

    snprintf(out_fn, 256, "%s.ref.bam", output_prefix);
    samFile *ref_out = sam_open(out_fn, "wb"); // write bam
    if (ref_out == NULL) { exit_err("failed to open BAM file %s\n", out_fn); }
    if (sam_hdr_write(ref_out, bam_header) != 0) { exit_err("bad header write %s\n", out_fn); } // write bam header

    snprintf(out_fn, 256, "%s.alt.bam", output_prefix);
    samFile *alt_out = sam_open(out_fn, "wb"); // write bam
    if (alt_out == NULL) { exit_err("failed to open BAM file %s\n", out_fn); }
    if (sam_hdr_write(alt_out, bam_header) != 0) { exit_err("bad header write %s\n", out_fn); } // write bam header

    snprintf(out_fn, 256, "%s.mul.bam", output_prefix);
    samFile *mul_out = sam_open(out_fn, "wb"); // write bam
    if (mul_out == NULL) { exit_err("failed to open BAM file %s\n", out_fn); }
    if (sam_hdr_write(mul_out, bam_header) != 0) { exit_err("bad header write %s\n", out_fn); } // write bam header

    snprintf(out_fn, 256, "%s.unk.bam", output_prefix);
    samFile *unk_out = sam_open(out_fn, "wb"); // write bam
    if (unk_out == NULL) { exit_err("failed to open BAM file %s\n", out_fn); }
    if (sam_hdr_write(unk_out, bam_header) != 0) { exit_err("bad header write %s\n", out_fn); } // write bam header

    bam1_t *aln = bam_init1(); // initialize an alignment
    while (sam_read1(sam_in, bam_header, aln) >= 0) {
        if (aln->core.tid < 0) continue; // not mapped

        /* Mapped & Primary alignments only */
        int n;
        int is_secondary = 0;
        int is_read2 = 0;
        char *flag = bam_flag2str(aln->core.flag);
        char *s, token[strlen(flag) + 1];
        for (s = flag; sscanf(s, "%[^,]%n", token, &n) == 1; s += n + 1) {
            if (strcmp("READ2", token) == 0) is_read2 = 1;
            else if (strcmp("SECONDARY", token) == 0 || strcmp("SUPPLEMENTARY", token) == 0) is_secondary = 1;
            if (*(s + n) != ',') break;
        }
        if (pao && is_secondary) {
            free(flag); flag = NULL;
            continue;
        }
        /* Write reads to appropriate file */
        out = NULL;
        char *name = (char *)aln->data;

        if (paired) is_read2 = 0;
        char key[256];
        snprintf(key, 256, "%s\t%d", name, is_read2);

        khiter_t k = kh_get(rh, read_hash, key);
        if (k != kh_end(read_hash)) {
            vector_t *node = &kh_val(read_hash, k);
            read_t **r = (read_t **)node->data;                                                                                                                          
            for (i = 0; i < node->len; i++) {
                if (strcmp(r[i]->name, name) == 0 && strcmp(r[i]->qseq, key) == 0) {
                    if (r[i]->index == 0 && !reverse) out = ref_out;
                    else if (r[i]->index == 1 && reverse) out = ref_out; // reverse, ALT writes to ref.bam
                    else if (!refonly && r[i]->index == 1 && !reverse) out = alt_out;
                    else if (!refonly && r[i]->index == 0 && reverse) out = alt_out; // reverse, REF writes to alt.bam
                    else if (!refonly && r[i]->index == 2) out = ref_out;
                    else if (!refonly && r[i]->index == 3) out = mul_out;
                    else if (r[i]->index == 4) out = unk_out;
                    if (debug >= 1) {
                        fprintf(stderr, "%f\t%f\t%f\t%d\t", r[i]->prgu, r[i]->prgv, r[i]->pout, r[i]->index);
                        fprintf(stderr, "%s\t%s\t%d\t%d\t", r[i]->name, r[i]->chr, r[i]->pos, r[i]->end);
                        fprintf(stderr, "%s\t%s\t%p\n", r[i]->flag, r[i]->qseq, out);
                    }
                    break;
                }
            }
            if (i == node->len) { exit_err("failed to find %s in hash key %d\n", key, k); }
        }
        free(flag); flag = NULL;

        if (other_bam != NULL && out == NULL) {
            int unique = 1;
            khiter_t k = kh_get(orh, other_read_hash, name);
            if (k != kh_end(other_read_hash)) {
                vector_t *node = &kh_val(other_read_hash, k);
                for (i = 0; i < node->len; i++) {
                    if (strcmp(node->data[i], name) == 0) {
                        unique = 0;
                        break;
                    }
                }
                if (i == node->len) { exit_err("failed to find %s in hash key %d\n", key, k); }
            }
            if (unique) out = ref_out;
        }

        if (out != NULL) {
            int r = sam_write1(out, bam_header, aln);
            if (r < 0) { exit_err("Bad program call"); }
        }
    }

    khiter_t k;
    for (k = kh_begin(other_read_hash); k != kh_end(other_read_hash); k++) {
        if (kh_exist(other_read_hash, k)) vector_destroy(&kh_val(other_read_hash, k));
    }
    kh_destroy(orh, other_read_hash);

    bam_destroy1(aln);
    bam_hdr_destroy(bam_header);
    sam_close(sam_in);

    out = NULL;
    sam_close(ref_out);
    sam_close(alt_out);
    sam_close(mul_out);
    sam_close(unk_out);
    print_status("# BAM Processed:\t%s\t%s", bam_file, asctime(time_info));
}

static int type2ind(char *type) {
    if (strcmp("REF", type) == 0) return 0;
    else if (strcmp("ALT", type) == 0) return 1;
    else if (strcmp("RLA", type) == 0) return 2;
    else if (strcmp("MUL", type) == 0) return 3;
    else if (strcmp("UNK", type) == 0) return 4;
    return -1;
}

static void readlist_read(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) { exit_err("failed to open file %s\n", filename); }

    int nreads = 0;
    char *line = NULL;
    ssize_t read_file = 0;
    size_t line_length = 0;
    while ((read_file = getline(&line, &line_length, file)) != -1) {
        if (line_length <= 0 || line[strspn(line, " \t\v\r\n")] == '\0') continue; // blank line
        if (line[0] == '#') continue;

        int pos;
        double prgu, prgv, pout;
        char type[line_length], name[line_length], flag[line_length], chr[line_length];
        int t = sscanf(line, "%s %s %s %d %lf %lf %lf %s %*[^\n]", name, type, chr, &pos, &prgu, &prgv, &pout, flag);
        if (t < 5) { exit_err("bad fields in read classified list file\n"); }

        int n;
        int is_read2 = 0;
        if (flag[0] == '-') paired = 1;
        if (!paired) {
            char *s, token[strlen(flag) + 1];
            for (s = flag; sscanf(s, "%[^,]%n", token, &n) == 1; s += n + 1) {
                if (strcmp("READ2", token) == 0) is_read2 = 1;
                if (*(s + n) != ',') break;
            }
        }

        char key[256];
        snprintf(key, 256, "%s\t%d", name, is_read2);

        khiter_t k = kh_get(rh, read_hash, key);
        if (k != kh_end(read_hash)) {
            size_t i;
            vector_t *node = &kh_val(read_hash, k);
            read_t **r = (read_t **)node->data;                                                                                                                          
            for (i = 0; i < node->len; i++) {
                if (strcmp(r[i]->name, name) == 0 && strcmp(r[i]->qseq, key) == 0) {
                    if (log_add_exp(prgu, prgv) > log_add_exp(r[i]->prgu, r[i]->prgv)) r[i]->index = type2ind(type);
                    r[i]->prgu = log_add_exp(r[i]->prgu, prgu);
                    r[i]->prgv = log_add_exp(r[i]->prgv, prgv);
                    r[i]->pout = log_add_exp(r[i]->pout, pout);
                    break;
                }
            }
            if (i == node->len) { exit_err("failed to find %s in hash key %d\n", key, k); }
        }
        else {
            read_t *r = read_create(name, 0, chr, pos);
            r->prgu = prgu;
            r->prgv = prgv;
            r->pout = pout;
            r->flag = strdup(flag);
            r->qseq = strdup(key);
            r->index = type2ind(type);

            int absent;
            khiter_t k = kh_put(rh, read_hash, r->qseq, &absent);
            vector_t *node = &kh_val(read_hash, k);
            if (absent) vector_init(node, 8, READ_T);
            vector_add(node, r);
            nreads++;
        }
    }
    free(line); line = NULL;
    fclose(file);
    print_status("# Classified list: %s\t%i reads\t%s", filename, nreads, asctime(time_info));
}

static void fasta_read(const char *fa_file) {
    faidx_t *fai = fai_load(fa_file);
    if (fai == NULL) {
        errno = fai_build(fa_file);
        if (errno == 0) { fai = fai_load(fa_file); }
        else { exit_err("failed to build and open FA index %s\n", fa_file); }
    }

    char *filename = malloc((strlen(fa_file) + 5) * sizeof (*filename));
    filename[0] = '\0';
    strcat(filename, fa_file);
    strcat(filename, ".fai");
    FILE *file = fopen(filename, "r");
    if (file == NULL) { exit_err("failed to open FA index for parsing %s\n", filename); }
    free(filename); filename = NULL;

    char *line = NULL;
    ssize_t read_file = 0;
    size_t line_length = 0;
    while ((read_file = getline(&line, &line_length, file)) != -1) {
        if (line_length <= 0 || line[strspn(line, " \t\v\r\n")] == '\0') continue; // blank line
        if (line[0] == '#') continue;

        char name[line_length];
        int t = sscanf(line, "%s%*[^ \t\v\r\n]", name);
        if (t < 1) { exit_err("bad fields in FA index file\n"); }
        if (!faidx_has_seq(fai, name)) { exit_err("failed to find %s in reference %s\n", name, fa_file); }

        fasta_t *f = fasta_create(name);
        f->seq = fai_fetch(fai, name, &f->seq_length);
        char *s;
        for (s = f->seq; *s != '\0'; s++) *s = toupper(*s);

        //u_int32_t hash = fnv_32a_str(name);
        //fprintf(stdout, "%s\t%u\n", name, hash);
        //continue;

        int absent;
        khiter_t k = kh_put(rsh, refseq_hash, f->name, &absent);
        vector_t *node = &kh_val(refseq_hash, k); // point to the bucket associated to k
        if (absent) vector_init(node, 8, FASTA_T);
        vector_add(node, f);
    }
    free(line); line = NULL;
    fclose(file);
    fai_destroy(fai);
    print_status("# Read reference genome: %s\t%s", fa_file, asctime(time_info));
}

static fasta_t *refseq_fetch(char *name) {
    size_t i;
	khiter_t k = kh_get(rsh, refseq_hash, name);
    if (k != kh_end(refseq_hash)) {
        vector_t *node = &kh_val(refseq_hash, k);
        fasta_t **f = (fasta_t **)node->data;
        for (i = 0; i < node->len; i++) {
            if (strcmp(f[i]->name, name) == 0) return f[i];
        }
        exit_err("failed to find %s in hash key %d\n", name, k);
    }
    return NULL;
}

static void bam_read(const char *bam_file, int ind) {
    samFile *sam_in = sam_open(bam_file, "r"); // open bam file
    if (sam_in == NULL) { exit_err("failed to open BAM file %s\n", bam_file); }
    bam_hdr_t *bam_header = sam_hdr_read(sam_in); // bam header
    if (bam_header == 0) { exit_err("bad header %s\n", bam_file); }

    int nreads = 0;
    bam1_t *aln = bam_init1(); // initialize an alignment
    while (sam_read1(sam_in, bam_header, aln) >= 0) {
        size_t i, j;
        if (aln->core.tid < 0) continue; // not mapped
        read_t *read = read_create((char *)aln->data, aln->core.tid, bam_header->target_name[aln->core.tid], aln->core.pos);
        fasta_t *f = refseq_fetch(read->chr);
        if (f == NULL) {
            read_destroy(read); free(read); read = NULL;
            continue;
        }
        char *refseq = f->seq;
        int refseq_length = f->seq_length;

        char *flag = bam_flag2str(aln->core.flag);
        if (flag != NULL) read->flag = strdup(flag);
        else read->flag = NULL;
        free(flag); flag = NULL;

        int n;
        char *s, token[strlen(read->flag) + 1];
        for (s = read->flag; sscanf(s, "%[^,]%n", token, &n) == 1; s += n + 1) {
            if (strcmp("DUP", token) == 0) read->is_dup = 1;
            else if (strcmp("REVERSE", token) == 0) read->is_reverse = 1;
            else if (strcmp("SECONDARY", token) == 0 || strcmp("SUPPLEMENTARY", token) == 0) read->is_secondary = 1;
            else if (strcmp("READ2", token) == 0) read->is_read2 = 1;
            if (*(s + n) != ',') break;
        }
        if ((nodup && read->is_dup) || (pao && read->is_secondary)) {
            read_destroy(read); free(read); read = NULL;
            continue;
        }

        int start_align = 0;
        int s_offset = 0; // offset for softclip at start
        int e_offset = 0; // offset for softclip at end

        u_int32_t *cigar = bam_get_cigar(aln);
        read->n_cigar = aln->core.n_cigar;
        read->cigar_oplen = malloc(read->n_cigar * sizeof (read->cigar_oplen));
        read->cigar_opchr = malloc((read->n_cigar + 1) * sizeof (read->cigar_opchr));
        read->splice_pos = malloc(read->n_cigar * sizeof (read->splice_pos));
        read->splice_offset = malloc(read->n_cigar * sizeof (read->splice_offset));

        j = 0;
        int splice_pos = 0; // track splice position in reads
        for (i = 0; i < read->n_cigar; i++) {
            read->cigar_oplen[i] = bam_cigar_oplen(cigar[i]);
            read->cigar_opchr[i] = bam_cigar_opchr(cigar[i]);
            read->splice_pos[i] = 0;
            read->splice_offset[i] = 0;

            if (read->cigar_opchr[i] == 'M' || read->cigar_opchr[i] == '=' || read->cigar_opchr[i] == 'X') start_align = 1;
            else if (start_align == 0 && read->cigar_opchr[i] == 'S') s_offset = read->cigar_oplen[i];
            else if (start_align == 1 && read->cigar_opchr[i] == 'S') e_offset = read->cigar_oplen[i];

            if (splice && read->cigar_opchr[i] == 'N') {
                read->splice_pos[j] = (isc) ? splice_pos - s_offset : splice_pos;
                read->splice_offset[j] = read->cigar_oplen[i];
                j++;
            }
            else if (splice && read->cigar_opchr[i] != 'D') {
                splice_pos += read->cigar_oplen[i];
            }

            if (read->cigar_opchr[i] != 'I') read->end += read->cigar_oplen[i];
        }
        read->cigar_opchr[read->n_cigar] = '\0';
        read->inferred_length = bam_cigar2qlen(read->n_cigar, cigar);
        read->n_splice = j;

        if (!isc) {
            read->pos -= s_offset; // compensate for soft clip in mapped position
            s_offset = 0;
            e_offset = 0;
        }
        else {
            read->end -= e_offset; // compensate for soft clip in mapped position
        }
        read->length = aln->core.l_qseq - (s_offset + e_offset);
        read->qseq = malloc((read->length + 1) * sizeof (read->qseq));
        read->qual = malloc(read->length  * sizeof (read->qual));
        uint8_t *qual = bam_get_qual(aln);
        for (i = 0; i < read->length; i++) {
            read->qseq[i] = toupper(seq_nt16_str[bam_seqi(bam_get_seq(aln), i + s_offset)]); // get nucleotide id and convert into IUPAC id.
            if (const_qual > 0) read->qual[i] = const_qual;
            else read->qual[i] = (phred64) ? qual[i] - 31 : qual[i]; // account for phred64
        }
        read->qseq[read->length] = '\0';

        read->multimapXA = NULL;
        if (bam_aux_get(aln, "XA")) read->multimapXA = strdup(bam_aux2Z(bam_aux_get(aln, "XA")));

        read->multimapNH = 1;
        if (bam_aux_get(aln, "NH")) read->multimapNH = bam_aux2i(bam_aux_get(aln, "NH"));

        double is_match[read->length], no_match[read->length];
        for (i = 0; i < read->length; i++) {
            is_match[i] = p_match[read->qual[i]];
            no_match[i] = p_mismatch[read->qual[i]];
        }
        /* Read probability matrix */
        double readprobmatrix[NT_CODES * read->length];
        set_prob_matrix(readprobmatrix, read, is_match, no_match, seqnt_map, bisulfite);
        /* Outside Paralog */
        double delta[read->length];
        for (i = 0; i < read->length; i++) delta[i] = no_match[i] - is_match[i];
        double a = sum_d(is_match, read->length);
        double elsewhere = log_add_exp(a, a + log_sum_exp(delta, read->length)) - (LOGALPHA * (read->length - read->inferred_length));
        double prgu = -DBL_MAX;
        double prgv = -DBL_MAX;
        double pout = elsewhere;

        if (ind == 0) prgu = calc_prob(readprobmatrix, read->length, refseq, refseq_length, read->pos, read->splice_pos, read->splice_offset, read->n_splice, seqnt_map);
        else if (ind == 1) prgv = calc_prob(readprobmatrix, read->length, refseq, refseq_length, read->pos, read->splice_pos, read->splice_offset, read->n_splice, seqnt_map);

        // Assuming all secondary alignments, corresponding to multi-map tags, are outputted to the bam file and will be processed eventually
        pout += lgomega;
        prgu = log_add_exp(pout, prgu);
        prgv = log_add_exp(pout, prgv);

        char key[256];
        snprintf(key, 256, "%s\t%d", read->name, read->is_read2);

        if (debug >= 2) {
            fprintf(stderr, "%f\t%f\t%f\t", prgu, prgv, pout);
            fprintf(stderr, "%s\t%s\t%d\t%d\t", read->name, read->chr, read->pos, read->end);
            for (i = 0; i < read->n_cigar; i++) fprintf(stderr, "%d%c ", read->cigar_oplen[i], read->cigar_opchr[i]);
            fprintf(stderr, "\t");
            if (read->multimapXA != NULL) fprintf(stderr, "%s\t", read->multimapXA);
            else fprintf(stderr, "%d\t", read->multimapNH);
            if (read->flag != NULL) fprintf(stderr, "%s\t", read->flag);
            fprintf(stderr, "%s\n", key);
        }

        khiter_t k = kh_get(rh, read_hash, key);
        if (k != kh_end(read_hash)) {
            vector_t *node = &kh_val(read_hash, k);
            read_t **r = (read_t **)node->data;
            for (i = 0; i < node->len; i++) {
                if (strcmp(r[i]->name, read->name) == 0 && strcmp(r[i]->qseq, key) == 0) { // seen before
                    r[i]->prgu = log_add_exp(r[i]->prgu, prgu);
                    r[i]->prgv = log_add_exp(r[i]->prgv, prgv);
                    r[i]->pout = log_add_exp(r[i]->pout, pout);
                    nreads++;
                    break;
                }
            }
            if (i == node->len) { exit_err("failed to find %s in hash key %d\n", key, k); }
        }
        else {
            read_t *r = read_create(read->name, read->tid, read->chr, read->pos);
            r->prgu = prgu;
            r->prgv = prgv;
            r->pout = pout;
            r->flag = strdup(read->flag);
            r->qseq = strdup(key);

            int absent;
            k = kh_put(rh, read_hash, r->qseq, &absent);
            vector_t *node = &kh_val(read_hash, k);
            if (absent) vector_init(node, 8, READ_T);
            vector_add(node, r);
            nreads++;
        }
        read_destroy(read); free(read); read = NULL;
        if (nreads % 1000000 == 0 ) { print_status("# Read bam:\t%s\t%d reads\t%s", bam_file, nreads, asctime(time_info)); }
    }
    bam_destroy1(aln);
    bam_hdr_destroy(bam_header);
    sam_close(sam_in);
    print_status("# Read bam:\t%s\t%d reads\t%s", bam_file, nreads, asctime(time_info));
}

static void combine_pe() {
	khiter_t k;
    size_t i, readi;
    other_read_hash = kh_init(orh);
    for (k = kh_begin(read_hash); k != kh_end(read_hash); k++) {
		if (kh_exist(read_hash, k)) {
            vector_t *node = &kh_val(read_hash, k);
            read_t **r = (read_t **)node->data;
            for (readi = 0; readi < node->len; readi++) {
                char key[256];
                snprintf(key, 256, "%s\t0", r[readi]->name);

                khiter_t k2 = kh_get(orh, other_read_hash, key);
                if (k2 != kh_end(other_read_hash)) {
                    vector_t *node2 = &kh_val(other_read_hash, k2);
                    read_t **r2 = (read_t **)node2->data;
                    for (i = 0; i < node2->len; i++) {
                        if (strcmp(r2[i]->name, r[readi]->name) == 0 && strcmp(r2[i]->qseq, key) == 0) {
                            r2[i]->prgu += r[readi]->prgu;
                            r2[i]->prgv += r[readi]->prgv;
                            r2[i]->pout += r[readi]->pout;
                            char flag[256];
                            snprintf(flag, 256, "%s;%s", r2[i]->flag, r[readi]->flag);
                            free(r2[i]->flag); r2[i]->flag = NULL;
                            r2[i]->flag = strdup(flag);
                            break;
                        }
                    }
                    if (i == node2->len) { exit_err("failed to find %s in hash key %d\n", key, k2); }
                }
                else {
                    read_t *r2 = read_create(r[readi]->name, r[readi]->tid, r[readi]->chr, r[readi]->pos);
                    r2->prgu = r[readi]->prgu;
                    r2->prgv = r[readi]->prgv;
                    r2->pout = r[readi]->pout;
                    r2->flag = strdup(r[readi]->flag);
                    r2->qseq = strdup(key);

                    int absent;
                    k2 = kh_put(orh, other_read_hash, r2->qseq, &absent);
                    vector_t *node2 = &kh_val(other_read_hash, k2);
                    if (absent) vector_init(node2, 8, READ_T);
                    vector_add(node2, r2);
                }
            }
            vector_destroy(&kh_val(read_hash, k));
        }
    }
    kh_destroy(rh, read_hash);

    read_hash = kh_init(rh);
    for (k = kh_begin(other_read_hash); k != kh_end(other_read_hash); k++) {
        if (kh_exist(other_read_hash, k)) {
            vector_t *node = &kh_val(other_read_hash, k);
            read_t **r = (read_t **)node->data;
            for (readi = 0; readi < node->len; readi++) {
                char key[256];
                snprintf(key, 256, "%s\t0", r[readi]->name);

                khiter_t k2 = kh_get(rh, read_hash, key);
                if (k2 == kh_end(read_hash)) {
                    read_t *r2 = read_create(r[readi]->name, r[readi]->tid, r[readi]->chr, r[readi]->pos);
                    r2->prgu = r[readi]->prgu;
                    r2->prgv = r[readi]->prgv;
                    r2->pout = r[readi]->pout;
                    r2->flag = strdup(r[readi]->flag);
                    r2->qseq = strdup(key);

                    int absent;
                    k2 = kh_put(rh, read_hash, r2->qseq, &absent);
                    vector_t *node2 = &kh_val(read_hash, k2);
                    if (absent) vector_init(node2, 8, READ_T);
                    vector_add(node2, r2);
                }
            }
            vector_destroy(&kh_val(other_read_hash, k));
        }
    }
    kh_destroy(orh, other_read_hash);
}

static void print_usage() {
    printf("\n");
    printf("Usage:\n");
    printf("(Default mode) with genotype info: eagle-rc [options] -a align.bam -o out -v eagle.out.txt eagle.readinfo.txt > classified_reads.list\n");
    printf("Options:\n");
    printf("  -v --var       FILE             EAGLE output text with variant likelihood ratios\n");
    printf("  -a --bam       FILE             Alignment data BAM file whose reads are to be classified\n");
    printf("  -o --out       String           Prefix for output BAM files\n");
    printf("  -u --unique    FILE1,FILE2,...  Optionally, also output reads that are unique against other BAM files (comma separated list)\n");
    printf("     --listonly                   Print classified read list only (stdout) without processing BAM file\n");
    printf("     --readlist                   Read from classified read list file instead of EAGLE outputs and proccess BAM file\n");
    printf("     --reclassify                 Reclassify after reading in classified read list file\n");
    printf("     --refonly                    Write REF classified reads only when processing BAM file\n");
    printf("     --paired                     Consider paired-end reads together.\n");
    printf("     --pao                        Primary alignments only.\n");
    printf("\nNo genotype info mode: eagle-rc [options] --ngi --ref1=ref1.fa --ref2=ref2.fa --bam1=align1.bam --bam2=align2.bam -o out > classified_reads.list\n");
    printf("Options (the above default mode options are also applicable):\n");
    printf("     --ngi                         No genotype information (i.e. vcf).  Directly classify read alignments mapped to two different reference (sub)genomes.\n");
    printf("     --ref1      FILE             --ngi mode: Reference genome 1 fasta file\n");
    printf("     --bam1      FILE             --ngi mode: Alignments to reference genome 1 bam file\n");
    printf("     --ref2      FILE             --ngi mode: Reference genome 2, recommend that sequence names are different from ref1\n");
    printf("     --bam2      FILE             --ngi mode: Alignments to reference genome 2 bam file\n");
    printf("     --isc                        --ngi mode: Ignore soft-clipped bases.\n");
    printf("     --nodup                      --ngi mode: Ignore marked duplicate reads (based on SAM flag).\n");
    printf("     --splice                     --ngi mode: RNA-seq spliced reads.\n");
    printf("     --bs        INT              --ngi mode: Bisulfite treated reads. 0: off, 1: top/forward strand, 2: bottom/reverse strand, 3: both. [0]\n");
    printf("     --phred64                    --ngi mode: Read quality scores are in phred64.\n");
    printf("     --omega     FLOAT            --ngi mode: Prior probability of originating from outside paralogous source, between [0,1]. [1e-40]\n");
    printf("     --cq        INT              --ngi mode: Constant quality as a phred score, ignoring the quality field in SAM. [0 is off]\n");
}

int main(int argc, char **argv) {
    /* Command line parameters defaults */
    debug = 0;
    char *readinfo_file = NULL;
    char *var_file = NULL;
    char *bam_file = NULL;
    char *output_prefix = NULL;
    char *other_bam = NULL;
    listonly = 0;
    readlist = 0;
    reclassify = 0;
    refonly = 0;
    paired = 0;
    pao = 0;
    isc = 0;
    nodup = 0;
    splice = 0;
    phred64 = 0;
    bisulfite = 0;
    const_qual = 0;
    reclassify = 0;

    ngi = 0;
    omega = 1.0e-40;
    char *bam_file1 = NULL;
    char *bam_file2 = NULL;
    char *ref_file1 = NULL;
    char *ref_file2 = NULL;

    static struct option long_options[] = {
        {"debug", optional_argument, NULL, 'd'},
        {"var", optional_argument, NULL, 'v'},
        {"bam", optional_argument, NULL, 'a'},
        {"out", optional_argument, NULL, 'o'},
        {"unique", optional_argument, NULL, 'u'},
        {"listonly", no_argument, &listonly, 1},
        {"readlist", no_argument, &readlist, 1},
        {"reclassify", no_argument, &reclassify, 1},
        {"refonly", no_argument, &refonly, 1},
        {"paired", no_argument, &paired, 1},
        {"pao", no_argument, &pao, 1},
        {"ngi", no_argument, &ngi, 1},
        {"isc", no_argument, &isc, 1},
        {"nodup", no_argument, &nodup, 1},
        {"splice", no_argument, &splice, 1},
        {"phred64", no_argument, &phred64, 1},
        {"ref1", optional_argument, NULL, 981},
        {"ref2", optional_argument, NULL, 982},
        {"bam1", optional_argument, NULL, 983},
        {"bam2", optional_argument, NULL, 984},
        {"omega", optional_argument, NULL, 991},
        {"bs", optional_argument, NULL, 992},
        {"cq", optional_argument, NULL, 993},
        {0, 0, 0, 0}
    };

    int opt = 0;
    while ((opt = getopt_long(argc, argv, "d:v:a:o:u:", long_options, &opt)) != -1) {
        switch (opt) {
            case 0: 
                //if (long_options[option_index].flag != 0) break;
                break;
            case 'd': debug = parse_int(optarg); break;
            case 'v': var_file = optarg; break;
            case 'a': bam_file = optarg; break;
            case 'o': output_prefix = optarg; break;
            case 'u': other_bam = optarg; break;
            case 981: ref_file1 = optarg; break;
            case 982: ref_file2 = optarg; break;
            case 983: bam_file1 = optarg; break;
            case 984: bam_file2 = optarg; break;
            case 991: omega = parse_float(optarg); break;
            case 992: bisulfite = parse_int(optarg); break;
            case 993: const_qual = parse_int(optarg); break;
            default: exit_usage("Bad options");
        }
    }
    if (optind > argc) { exit_usage("Bad program call"); }

    if (!listonly && !ngi && bam_file == NULL) { exit_usage("Missing BAM file! -a bam"); }
    else if (!listonly && output_prefix == NULL) { exit_usage("Missing output prefix!"); }

    print_status("# Options: listonly=%d readlist=%d reclassify=%d refonly=%d paired=%d pao=%d\n", listonly, readlist, reclassify, refonly, paired, pao);
    print_status("#          ngi=%d isc=%d nodup=%d splice=%d bs=%d phred64=%d omega=%g cq=%d\n", ngi, isc, nodup, splice, bisulfite, phred64, omega, const_qual);
    print_status("# Start: \t%s", asctime(time_info));

    /* Start processing data */
    clock_t tic = clock();

    var_hash = kh_init(vh);
    read_hash = kh_init(rh);

    if (ngi) {
        init_seqnt_map(seqnt_map);
        init_q2p_table(p_match, p_mismatch, 50);

        if (ref_file1 == NULL || ref_file2 == NULL) { exit_usage("Missing reference FASTA file! ref1 or ref2!"); }
        if (bam_file1 == NULL || bam_file2 == NULL) { exit_usage("Missing BAM file! bam1 or bam2!"); }
        if (omega < 0 || omega > 1) omega = 1e-40;
        lgomega = (log(omega) - log(1.0-omega));

        khiter_t k;
        refseq_hash = kh_init(rsh);
        fasta_read(ref_file1);
        bam_read(bam_file1, 0);
        for (k = kh_begin(refseq_hash); k != kh_end(refseq_hash); k++) {
            if (kh_exist(refseq_hash, k)) vector_destroy(&kh_val(refseq_hash, k));
        }
        kh_destroy(rsh, refseq_hash);

        refseq_hash = kh_init(rsh);
        fasta_read(ref_file2);
        bam_read(bam_file2, 1);
        for (k = kh_begin(refseq_hash); k != kh_end(refseq_hash); k++) {
            if (kh_exist(refseq_hash, k)) vector_destroy(&kh_val(refseq_hash, k));
        }
        kh_destroy(rsh, refseq_hash);

        if (paired) combine_pe();
        readinfo_classify();

        if (!listonly) {
            char output_prefix1[256], output_prefix2[256];
            snprintf(output_prefix1, 256, "%s1", output_prefix);
            snprintf(output_prefix2, 256, "%s2", output_prefix);

            bam_write(bam_file1, output_prefix1, other_bam, 0);
            bam_write(bam_file2, output_prefix2, other_bam, 1);
        }
    }
    else if (readlist) {
        readlist_read(argv[optind]);
        if (paired) combine_pe();
        if (paired || reclassify) readinfo_classify();
        if (!listonly) bam_write(bam_file, output_prefix, other_bam, 0);
    }
    else {
        //var_file = argv[optind++];
        if (var_file == NULL) { exit_usage("Missing EAGLE output file!"); } 

        readinfo_file = argv[optind++];
        if (readinfo_file == NULL) { exit_usage("Missing EAGLE read info file!"); } 

        int nvars = var_read(var_file);
        print_status("# Read EAGLE variants file: %s\t%i entries\t%s", var_file, nvars, asctime(time_info));

        int nreads = readinfo_read(readinfo_file);
        print_status("# Read EAGLE read info file: %s\t%i reads\t%s", readinfo_file, nreads, asctime(time_info));

        if (paired) combine_pe();
        readinfo_classify();
        if (!listonly) bam_write(bam_file, output_prefix, other_bam, 0);
    }

    khiter_t k;
    for (k = kh_begin(var_hash); k != kh_end(var_hash); k++) {
        if (kh_exist(var_hash, k)) vector_destroy(&kh_val(var_hash, k));
    }
    kh_destroy(vh, var_hash);

    for (k = kh_begin(read_hash); k != kh_end(read_hash); k++) {
        if (kh_exist(read_hash, k)) vector_destroy(&kh_val(read_hash, k));
    }
    kh_destroy(rh, read_hash);

    clock_t toc = clock();
    print_status("# CPU time (hr):\t%f\n", (double)(toc - tic) / CLOCKS_PER_SEC / 3600);

    return 0;
}
