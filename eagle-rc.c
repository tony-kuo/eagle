/*
EAGLE: explicit alternative genome likelihood evaluator
Utility program that classifies reads based on EAGLE calculated likelihoods (from verbose output)

Run EAGLE with --rc which sets: --omega=1e-40 --mvh --verbose --isc
  where main variant evaluation output goes to STOUT
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
#define VERSION "1.1.1"
#define ALPHA 1.3     // Factor to account for longer read lengths lowering the probability a sequence matching an outside paralogous source
#define LOGALPHA (log(ALPHA))

/* Command line arguments */
static int debug;
static int listonly;
static int readlist;
static int reclassify;
static int refonly;
static int paired, already_paired;
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

KHASH_MAP_INIT_STR(rh, read_t *) // hashmap: string key, read_t * value
static khash_t(rh) *read_hash; // pointer to hashmap, read

KHASH_MAP_INIT_STR(orh, read_t *) // hashmap: string key, read_t * value
static khash_t(orh) *other_read_hash; // pointer to hashmap, reads from other_bam files

KHASH_MAP_INIT_STR(vh, char) // hashmap: string key, char value
static khash_t(vh) *var_hash;  // pointer to hashmap, variant

KHASH_MAP_INIT_STR(rsh, fasta_t *)   // hashmap: string key, fasta_t * value
static khash_t(rsh) *refseq_hash; // pointer to hashmap

static void add2var_list(vector_t *var_list, char *set) {
    char var[strlen(set)];

    int n;
    char *s;
    for (s = set + 1; sscanf(s, "%[^;];%n", var, &n) == 1; s += n) { // scan variant set
        char *v = strdup(var);
        size_t i;
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
        if (t < 5) { exit_err("bad fields in EAGLE output file\n%s\n", line); }

        if (strcmp(set, "[]") == 0) snprintf(set, line_length, "[%s,%d,%s,%s;]", chr, pos, ref, alt);

        char *s;
        for (s = set + 1; sscanf(s, "%[^;];%n", chr, &n) == 1; s += n) { // scan variant set and add to hash
            int absent;
            khiter_t k = kh_put(vh, var_hash, chr, &absent);
            if (absent) kh_key(var_hash, k) = strdup(chr);
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

        char key[strlen(name) + 3];
        snprintf(key, strlen(name) + 3, "%s\t%d", name, is_read2);

        int absent;
        khiter_t k = kh_put(rh, read_hash, key, &absent);
        if (absent) {
            read_t *r = read_create(name, 0, chr, pos);
            r->prgu = (float)prgu;
            r->prgv = (float)prgv;
            r->pout = (float)pout;
            r->flag = strdup(flag);
            r->qseq = strdup(key);
            add2var_list(r->var_list, var);
            kh_key(read_hash, k) = r->qseq;
            kh_val(read_hash, k) = r;
            nreads++;
        }
        else {
            read_t *r = kh_val(read_hash, k);
            r->prgu = (float)log_add_exp((double)r->prgu, prgu);
            r->prgv = (float)log_add_exp((double)r->prgv, prgv);
            add2var_list(r->var_list, var);
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

    size_t i;
    khiter_t k;
    for (k = kh_begin(read_hash); k != kh_end(read_hash); k++) {
		if (kh_exist(read_hash, k)) {
            read_t *r = kh_val(read_hash, k);
            char **v = (char **)r->var_list->data;
            size_t nvariants = r->var_list->len;

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
                    r->index = 3;
                    fprintf(stdout, "%s\tMUL\t%s\t%d\t%f\t%f\t%f\t%s\t", r->name, r->chr, r->pos, r->prgu, r->prgv, r->pout, r->flag);
                    for (i = 0; i < nvariants; i++) { fprintf(stdout, "%s;", v[i]); } fprintf(stdout, "\n");
                    continue;
                }
            }

            // if any variant in list is not in EAGLE output, it suggests variants are from a different phase if alt wins
            multiallele = 0;
            for (i = 0; i < nvariants; i++) {
                khiter_t k = kh_get(vh, var_hash, v[i]);
                if (k != kh_end(var_hash)) multiallele = 1;
            }

            fprintf(stdout, "%s\t", r->name);
            if (r->prgu > r->prgv && r->prgu - r->prgv >= 0.69) { // ref wins
                r->index = 0;
                fprintf(stdout, "REF\t");
            }
            else if (r->prgv > r->prgu && r->prgv - r->prgu >= 0.69) { // alt wins
                if (multiallele) { // EAGLE outputs the set with highest likelihood ratio, i.e. most different from reference, leaving the "reference-like-allele"
                    r->index = 2;
                    fprintf(stdout, "RLA\t");
                }
                else {
                    r->index = 1;
                    fprintf(stdout, "ALT\t");
                }
            }
            else { // unknown
                r->index = 4;
                fprintf(stdout, "UNK\t");
            }
            fprintf(stdout, "%s\t%d\t%f\t%f\t%f\t%s\t", r->chr, r->pos, r->prgu, r->prgv, r->pout, r->flag);
            for (i = 0; i < nvariants; i++) { fprintf(stdout, "%s;", v[i]); } fprintf(stdout, "\n");
        }
    }
    fflush(stdout);
    print_status("# Reads Classified:\t%s", asctime(time_info));
}

static void bam_write(const char *bam_file, const char *output_prefix, char *other_bam, int reverse) {
    khiter_t k;
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

                int absent;
                k = kh_put(orh, other_read_hash, (char *)aln->data, &absent);
                if (absent) {
                    read_t *r = read_create((char *)aln->data, aln->core.tid, bam_header->target_name[aln->core.tid], aln->core.pos);
                    kh_key(other_read_hash, k) = r->name;
                    kh_val(other_read_hash, k) = r;
                }
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
    print_status("# Open input bam:\t%s\t%s", bam_file, asctime(time_info));

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

        int is_secondary = 0;
        int is_read2 = 0;
        char *flag = bam_flag2str(aln->core.flag);

        int n;
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

        if (paired || already_paired) is_read2 = 0;
        char key[strlen(name) + 3];
        snprintf(key, strlen(name) + 3, "%s\t%d", name, is_read2);

        if (other_bam != NULL) {
            k = kh_get(orh, other_read_hash, name);
            if (k == kh_end(other_read_hash)) out = ref_out; // unique vs other bams
        }
        else {
            k = kh_get(rh, read_hash, key);
            if (k != kh_end(read_hash)) {
                read_t *r = kh_val(read_hash, k);
                if (r->index == 0 && !reverse) out = ref_out;
                else if (r->index == 1 && reverse) out = ref_out; // reverse, ALT writes to ref.bam
                else if (!refonly && r->index == 1 && !reverse) out = alt_out;
                else if (!refonly && r->index == 0 && reverse) out = alt_out; // reverse, REF writes to alt.bam
                else if (!refonly && r->index == 2) out = ref_out;
                else if (!refonly && r->index == 3) out = mul_out;
                else if (r->index == 4) out = unk_out;
                if (debug >= 1) {
                    fprintf(stderr, "%f\t%f\t%f\t%d\t", r->prgu, r->prgv, r->pout, r->index);
                    fprintf(stderr, "%s\t%s\t%d\t", r->name, r->chr, r->pos);
                    fprintf(stderr, "%s\t%s\t%p\n", r->flag, r->qseq, out);
                }
            }
        }
        free(flag); flag = NULL;

        if (out != NULL) {
            int r = sam_write1(out, bam_header, aln);
            if (r < 0) { exit_err("Bad program call"); }
        }
    }

    for (k = kh_begin(other_read_hash); k != kh_end(other_read_hash); k++) {
        if (kh_exist(other_read_hash, k)) {
            //print_status("# key to exclude, from other BAM:\t%s", kh_key(other_read_hash, k));
            read_destroy(kh_val(other_read_hash, k)); free(kh_val(other_read_hash, k)); kh_val(other_read_hash, k) = NULL;
        }
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

static int readlist_read(FILE *file) {
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
        if (t < 8) {
            t = sscanf(line, "%s %s %d %lf %lf %lf %*[^\t] %*[^\n]", name, chr, &pos, &prgu, &prgv, &pout);
            if (t < 7) { exit_err("bad fields in EAGLE read info file\n%s\n", line); }
            flag[0] = '\0';
        }

        int n;
        int is_read2 = 0;
        if (flag[0] == '-') already_paired = 1;

        char *s, token[strlen(flag) + 1];
        for (s = flag; sscanf(s, "%[^,]%n", token, &n) == 1; s += n + 1) {
            if (strcmp("READ2", token) == 0) is_read2 = 1;
            if (*(s + n) != ',') break;
        }

        char key[strlen(name) + 3];
        snprintf(key, strlen(name) + 3, "%s\t%d", name, is_read2);

        int absent;
        khiter_t k = kh_put(rh, read_hash, key, &absent);
        if (absent) {
            read_t *r = read_create(name, 0, chr, pos);
            r->prgu = (float)prgu;
            r->prgv = (float)prgv;
            r->pout = (float)pout;
            r->flag = strdup(flag);
            r->qseq = strdup(key);
            r->index = type2ind(type);
            kh_key(read_hash, k) = r->qseq;
            kh_val(read_hash, k) = r;
            nreads++;
        }
        else {
            read_t *r = kh_val(read_hash, k);
            if (log_add_exp(prgu, prgv) > log_add_exp((double)r->prgu, (double)r->prgv)) {
                free(r->chr); r->chr = NULL;
                r->chr = strdup(chr);
                r->pos = pos;
                free(r->flag); r->flag = NULL;
                r->flag = strdup(flag);
                r->index = type2ind(type);
            }
            r->prgu = (float)log_add_exp((double)r->prgu, prgu);
            r->prgv = (float)log_add_exp((double)r->prgv, prgv);
            r->pout = (float)log_add_exp((double)r->pout, pout);
        }
    }
    free(line); line = NULL;
    fclose(file);
    return nreads;
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
        //f->seq = fai_fetch(fai, f->name, &f->seq_length);
        f->seq = faidx_fetch_seq(fai, f->name, 0, faidx_seq_len(fai, f->name) - 1, &f->seq_length);
        char *s;
        for (s = f->seq; *s != '\0'; s++) *s = toupper(*s);

        //u_int32_t hash = fnv_32a_str(name);
        //fprintf(stdout, "%s\t%u\n", name, hash);
        //continue;

        int absent;
        khiter_t k = kh_put(rsh, refseq_hash, f->name, &absent);
        if (absent) { kh_val(refseq_hash, k) = f; }
        else { exit_err("# refseq_hash collision: %s", asctime(time_info)); }
    }
    free(line); line = NULL;
    fclose(file);
    fai_destroy(fai);
    print_status("# Read reference genome: %s\t%s", fa_file, asctime(time_info));
}

static fasta_t *refseq_fetch(char *name) {
	khiter_t k = kh_get(rsh, refseq_hash, name);
    if (k != kh_end(refseq_hash)) {
        return kh_val(refseq_hash, k);
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
        int i;
        read_t *read = read_fetch(bam_header, aln, pao, isc, nodup, splice, phred64, const_qual);
        if (read == NULL) continue;

        fasta_t *f = refseq_fetch(read->chr);
        if (strcmp(read->chr, f->name) != 0) { exit_err("# request: %s\tvs\tfetch: %s\t%s", read->chr, f->name, asctime(time_info)); }
        if (f == NULL) {
            read_destroy(read); free(read); read = NULL;
            continue;
        }
        char *refseq = f->seq;
        int refseq_length = f->seq_length;

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

        char key[strlen(read->name) + 3];
        snprintf(key, strlen(read->name) + 3, "%s\t%d", read->name, read->is_read2);

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

        int absent;
        khiter_t k = kh_put(rh, read_hash, key, &absent);
        if (absent) {
            read_t *r = read_create(read->name, read->tid, read->chr, read->pos);
            r->prgu = (float)prgu;
            r->prgv = (float)prgv;
            r->pout = (float)pout;
            r->flag = strdup(read->flag);
            r->qseq = strdup(key);
            kh_key(read_hash, k) = r->qseq;
            kh_val(read_hash, k) = r;
        }
        else {
            read_t *r = kh_val(read_hash, k);
            if (log_add_exp(prgu, prgv) > log_add_exp((double)r->prgu, (double)r->prgv)) {
                r->tid = read->tid;
                free(r->chr); r->chr = NULL;
                r->chr = strdup(read->chr);
                r->pos = read->pos;
                free(r->flag); r->flag = NULL;
                r->flag = strdup(read->flag);
            }
            r->prgu = (float)log_add_exp((double)r->prgu, prgu);
            r->prgv = (float)log_add_exp((double)r->prgv, prgv);
            r->pout = (float)log_add_exp((double)r->pout, pout);
        }
        nreads++;
        read_destroy(read); free(read); read = NULL;
        if (nreads % 1000000 == 0 ) { print_status("# Read bam:\t%s\t%d reads\t%s", bam_file, nreads, asctime(time_info)); }
    }
    bam_destroy1(aln);
    bam_hdr_destroy(bam_header);
    sam_close(sam_in);
    print_status("# Read bam:\t%s\t%d reads\t%s", bam_file, nreads, asctime(time_info));
}

static void combine_pe() {
    other_read_hash = kh_init(orh);
	khiter_t k;
    for (k = kh_begin(read_hash); k != kh_end(read_hash); k++) {
		if (kh_exist(read_hash, k)) {
            read_t *r = kh_val(read_hash, k);
            char key[strlen(r->name) + 3];
            snprintf(key, strlen(r->name) + 3, "%s\t0", r->name);

            int absent;
            khiter_t k2 = kh_put(orh, other_read_hash, key, &absent);
            if (absent) {
                read_t *r2 = read_create(r->name, r->tid, r->chr, r->pos);
                r2->prgu = r->prgu;
                r2->prgv = r->prgv;
                r2->pout = r->pout;
                r2->flag = strdup(r->flag);
                r2->qseq = strdup(key);
                kh_key(other_read_hash, k2) = r2->qseq;
                kh_val(other_read_hash, k2) = r2;
            }
            else {
                read_t *r2 = kh_val(other_read_hash, k2);
                if (log_add_exp((double)r->prgu, (double)r->prgv) > log_add_exp((double)r2->prgu, (double)r2->prgv)) {
                    r2->tid = r->tid;
                    free(r2->chr); r2->chr = NULL;
                    r2->chr = strdup(r->chr);
                    r2->pos = r->pos;
                }
                r2->prgu += r->prgu;
                r2->prgv += r->prgv;
                r2->pout += r->pout;
                char flag[strlen(r2->flag) + strlen(r->flag) + 2];
                snprintf(flag, strlen(r2->flag) + strlen(r->flag) + 2, "%s;%s", r2->flag, r->flag);
                free(r2->flag); r2->flag = NULL;
                r2->flag = strdup(flag);
            }
            read_destroy(r); free(r); r = NULL;
        }
    }
    kh_destroy(rh, read_hash);

    read_hash = kh_init(rh);
    for (k = kh_begin(other_read_hash); k != kh_end(other_read_hash); k++) {
        if (kh_exist(other_read_hash, k)) {
            int absent;
            khiter_t k2 = kh_put(rh, read_hash, kh_key(other_read_hash, k), &absent);
            kh_val(read_hash, k2) = kh_val(other_read_hash, k); 
        }
    }
    kh_destroy(orh, other_read_hash);
}

static void print_usage() {
    printf("\n");
    printf("Usage:\n");
    printf("(Default mode) with genotype info: eagle-rc [options] -a align.bam -o out -v eagle.out.txt eagle.readinfo.txt > classified_reads.list\n");
    printf("Options:\n");
    printf("  -v --var       FILE             EAGLE output text with variant likelihood ratios.\n");
    printf("  -a --bam       FILE             Alignment data BAM file whose reads are to be classified.\n");
    printf("  -o --out       String           Prefix for output BAM files.\n");
    printf("  -u --unique    FILE1,FILE2,...  Optionally, also output reads that are unique against other BAM files (comma separated list).\n");
    printf("     --listonly                   Print classified read list only (stdout) without processing BAM file.\n");
    printf("     --readlist                   Read from classified read list file instead of EAGLE outputs and proccess BAM file.\n");
    printf("     --reclassify                 Reclassify after reading in classified read list file.\n");
    printf("     --refonly                    Write REF classified reads only when processing BAM file.\n");
    printf("     --paired                     Consider paired-end reads together.\n");
    printf("     --pao                        Primary alignments only.\n");
    printf("     --version                    Display version.\n");
    printf("\nNo genotype info mode: eagle-rc [options] --ngi --ref1=ref1.fa --ref2=ref2.fa --bam1=align1.bam --bam2=align2.bam -o out > classified_reads.list\n");
    printf("Options (the above default mode options are also applicable):\n");
    printf("     --ngi                         No genotype information (i.e. vcf).  Directly classify read alignments mapped to two different reference (sub)genomes.\n");
    printf("     --ref1      FILE             --ngi mode: Reference genome 1 fasta file.\n");
    printf("     --bam1      FILE             --ngi mode: Alignments to reference genome 1 bam file.\n");
    printf("     --ref2      FILE             --ngi mode: Reference genome 2, recommend that sequence names are different from ref1.\n");
    printf("     --bam2      FILE             --ngi mode: Alignments to reference genome 2 bam file.\n");
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
    already_paired = 0;
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
        {"version", optional_argument, NULL, 999},
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
            case 991: omega = parse_double(optarg); break;
            case 992: bisulfite = parse_int(optarg); break;
            case 993: const_qual = parse_int(optarg); break;
            case 999: printf("EAGLE-RC %s\n", VERSION); exit(0);
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
        for (k = kh_begin(ref_hash); k != kh_end(refseq_hash); k++) {
            if (kh_exist(refseq_hash, k)) {
                fasta_destroy(kh_val(refseq_hash, k)); free(kh_val(refseq_hash, k)); kh_val(refseq_hash, k) = NULL;
            }
        }
        kh_destroy(rsh, refseq_hash);

        refseq_hash = kh_init(rsh);
        fasta_read(ref_file2);
        bam_read(bam_file2, 1);
        for (k = kh_begin(ref_hash); k != kh_end(refseq_hash); k++) {
            if (kh_exist(refseq_hash, k)) {
                fasta_destroy(kh_val(refseq_hash, k)); free(kh_val(refseq_hash, k)); kh_val(refseq_hash, k) = NULL;
            }
        }
        kh_destroy(rsh, refseq_hash);

        if (paired) combine_pe();
        readinfo_classify();

        if (!listonly) {
            char output_prefix1[strlen(output_prefix) + 2], output_prefix2[strlen(output_prefix) + 2];
            snprintf(output_prefix1, strlen(output_prefix) + 2, "%s1", output_prefix);
            snprintf(output_prefix2, strlen(output_prefix) + 2, "%s2", output_prefix);

            bam_write(bam_file1, output_prefix1, other_bam, 0);
            bam_write(bam_file2, output_prefix2, other_bam, 1);
        }
    }
    else if (readlist) {
        FILE *list_fh = stdin;
        char *filename = argv[optind];
        if (filename != NULL) {
            list_fh = fopen(filename, "r");
            if (list_fh == NULL) { exit_err("failed to open file %s\n", filename); }
        }
        else {
            filename = "stdin";
        }
        int nreads = readlist_read(list_fh);
        print_status("# Classified list: %s\t%i reads\t%s", filename, nreads, asctime(time_info));

        if (paired) combine_pe();
        if (reclassify) readinfo_classify();
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
        if (kh_exist(var_hash, k)) {
            free((char*)kh_key(var_hash, k)); kh_key(var_hash, k) = NULL;
        }
    }
    kh_destroy(vh, var_hash);

    for (k = kh_begin(read_hash); k != kh_end(read_hash); k++) {
        if (kh_exist(read_hash, k)) {
            read_destroy(kh_val(read_hash, k)); free(kh_val(read_hash, k)); kh_val(read_hash, k) = NULL;
        }
    }
    kh_destroy(rh, read_hash);

    clock_t toc = clock();
    print_status("# CPU time (hr):\t%f\n", (double)(toc - tic) / CLOCKS_PER_SEC / 3600);

    return 0;
}
