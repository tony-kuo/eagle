/*
EAGLE: explicit alternative genome likelihood evaluator
Utility program that classifies reads based on EAGLE calculated likelihoods (from verbose output)

Run EAGLE with required options: --omega=1e-40 --mvh --verbose
  where main variant evaluation output goes to STDIN
  while per read likelihoods (verbose) go to STDERR

Other program options are situational

ex) eagle -t 2 -v var.vcf -a align.bam -r ref.fa --omega=1.0e-40 --mvh --verbose 1> out.txt  2> readinfo.txt

Copyright 2016 Tony Kuo
This program is distributed under the terms of the GNU General Public License
*/

#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "htslib/sam.h"
#include "htslib/khash.h"
#include "util.h"
#include "vector.h"

/* Command line arguments */
static int listonly;
static int readlist;
static int refonly;

/* Time info */
static time_t now; 
static struct tm *time_info; 
#define print_status(M, ...) time(&now); time_info = localtime(&now); fprintf(stderr, M, ##__VA_ARGS__);

KHASH_MAP_INIT_STR(rh, Vector) // hashmap: string key, vector value
static khash_t(rh) *read_hash; // pointer to hashmap, read

KHASH_MAP_INIT_STR(vh, Vector) // hashmap: string key, vector value
static khash_t(vh) *var_hash;  // pointer to hashmap, variant

static void add2var_list(Vector *var_list, char *set) {
    char var[strlen(set)];

    int n;
    char *s;
    size_t i;
    for (s = set + 1; sscanf(s, "%[^;];%n", var, &n) == 1; s += n) { // scan variant set
        char *v = strdup(var);
        for (i = 0; i < var_list->size; ++i) {
            if (strcmp(v, (char *)var_list->data[i]) == 0) break;
        }
        if (i == var_list->size) { vector_add(var_list, v); }
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
        if (t < 4) { exit_err("bad fields in EAGLE output file\n"); }

        if (strcmp(set, "[]") == 0) {
            snprintf(set, line_length, "[%s,%d,%s,%s;]", chr, pos, ref, alt);
        }

        char *s;
        for (s = set + 1; sscanf(s, "%[^;];%n", chr, &n) == 1; s += n) { // scan variant set and add to hash
            char *v = strdup(chr);

            Vector *node;
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

            ++nvars;
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

    size_t i;
    int nreads = 0;
    char *line = NULL;
    ssize_t read_file = 0;
    size_t line_length = 0;
    while ((read_file = getline(&line, &line_length, file)) != -1) {
        if (line_length <= 0 || line[strspn(line, " \t\v\r\n")] == '\0') continue; // blank line
        if (line[0] == '#') continue;

        int pos;
        double prgu, prgv, pout;
        char name[line_length], chr[line_length], set[line_length];
        int t = sscanf(line, "%s %s %d %lf %lf %lf %*[^\t] %*[^\t] %*[^\t] %s", name, chr, &pos, &prgu, &prgv, &pout, set);
        if (t < 7) { exit_err("bad fields in EAGLE read info file\n"); }

        khiter_t k = kh_get(rh, read_hash, name);
        if (k != kh_end(read_hash)) {
            Vector *node = &kh_val(read_hash, k);
            Read **r = (Read **)node->data;                                                                                                                          
            for (i = 0; i < node->size; ++i) {
                if (strcmp(r[i]->name, name) == 0) {
                    r[i]->prgu += prgu;
                    r[i]->prgv += prgv;
                    add2var_list(r[i]->var_list, set);
                    break;
                }
            }
            if (i == node->size) { exit_err("failed to find %s in hash key %d\n", name, k); }
        }
        else {
            Read *r = read_create(name, 0, chr, pos);
            r->prgu = prgu;
            r->prgv = prgv;
            r->pout = pout;
            add2var_list(r->var_list, set);

            int absent;
            k = kh_put(rh, read_hash, r->name, &absent);
            Vector *node = &kh_val(read_hash, k);
            if (absent) vector_init(node, 8, READ_T);
            vector_add(node, r);
            ++nreads;
        }
    }
    free(line); line = NULL;
    fclose(file);
    return nreads;
}

static void process_bam(const char *bam_file, const char *output_prefix) {
    samFile *sam_in = sam_open(bam_file, "r"); // open bam file
    if (sam_in == NULL) { exit_err("failed to open BAM file %s\n", bam_file); }
    bam_hdr_t *bam_header = sam_hdr_read(sam_in); // bam header
    if (bam_header == 0) { exit_err("bad header %s\n", bam_file); }

    /* Split input bam files into respective categories' output bam files */
    samFile *out;
    char out_fn[999];

    snprintf(out_fn, 999, "%s.ref.bam", output_prefix);
    samFile *ref_out = sam_open(out_fn, "wb"); // write bam
    if (ref_out == NULL) { exit_err("failed to open BAM file %s\n", out_fn); }
    if (sam_hdr_write(ref_out, bam_header) != 0) { exit_err("bad header write %s\n", out_fn); } // write bam header

    snprintf(out_fn, 999, "%s.alt.bam", output_prefix);
    samFile *alt_out = sam_open(out_fn, "wb"); // write bam
    if (alt_out == NULL) { exit_err("failed to open BAM file %s\n", out_fn); }
    if (sam_hdr_write(alt_out, bam_header) != 0) { exit_err("bad header write %s\n", out_fn); } // write bam header

    snprintf(out_fn, 999, "%s.mul.bam", output_prefix);
    samFile *mul_out = sam_open(out_fn, "wb"); // write bam
    if (mul_out == NULL) { exit_err("failed to open BAM file %s\n", out_fn); }
    if (sam_hdr_write(mul_out, bam_header) != 0) { exit_err("bad header write %s\n", out_fn); } // write bam header

    snprintf(out_fn, 999, "%s.com.bam", output_prefix);
    samFile *com_out = sam_open(out_fn, "wb"); // write bam
    if (com_out == NULL) { exit_err("failed to open BAM file %s\n", out_fn); }
    if (sam_hdr_write(com_out, bam_header) != 0) { exit_err("bad header write %s\n", out_fn); } // write bam header

    bam1_t *aln = bam_init1(); // initialize an alignment
    while (sam_read1(sam_in, bam_header, aln) >= 0) {
        /* Mapped & Primary alignments only */
        int n;
        int is_unmap = 0;
        int is_secondary = 0;
        char *flag = bam_flag2str(aln->core.flag);
        char *s, token[strlen(flag) + 1];
        for (s = flag; sscanf(s, "%[^,]%n", token, &n) == 1; s += n + 1) {
            if (strcmp("UNMAP", token) == 0) is_unmap = 1;
            else if (strcmp("SECONDARY", token) == 0 || strcmp("SUPPLEMENTARY", token) == 0) is_secondary = 1;
            if (*(s + n) != ',') break;
        }
        free(flag); flag = NULL;
        if (is_unmap || is_secondary) continue;

        /* Write reads to appropriate file */
        out = NULL;
        size_t i;
        char *name = (char *)aln->data;

        khiter_t k = kh_get(rh, read_hash, name);
        if (k != kh_end(read_hash)) {
            Vector *node = &kh_val(read_hash, k);
            Read **r = (Read **)node->data;                                                                                                                          
            for (i = 0; i < node->size; ++i) {
                if (strcmp(r[i]->name, name) == 0) {
                    if (r[i]->index == 0) out = ref_out;
                    else if (!refonly && r[i]->index == 1) out = alt_out;
                    else if (!refonly && r[i]->index == 2) out = ref_out;
                    else if (!refonly && r[i]->index == 3) out = mul_out;
                    else if (r[i]->index == 4) out = NULL;
                    break;
                }
            }
            if (i == node->size) { exit_err("failed to find %s in hash key %d\n", name, k); }
        }
        else if (!refonly) {
            out = com_out;
        }

        if (out != NULL) {
            int r = sam_write1(out, bam_header, aln);
            if (r < 0) { exit_err("Bad program call"); }
        }
    }
    bam_destroy1(aln);
    bam_hdr_destroy(bam_header);
    sam_close(sam_in);

    out = NULL;
    sam_close(ref_out);
    sam_close(alt_out);
    sam_close(mul_out);
    sam_close(com_out);
}

static void classify_reads(const char *bam_file, const char *output_prefix) {
    /* read->index set to classification:
    reference = 0
    alternative = 1
    reference-like-allele, multi-allelic, closer to the reference = 2;
    multi-allelic that are undifferentiateable = 3
    unknown, ambiguous with equal likelihoods for reference and alternative = 4
    */

    size_t readi, i, j;
	khiter_t k;
    for (k = kh_begin(read_hash); k != kh_end(read_hash); ++k) {
		if (kh_exist(read_hash, k)) {
            Vector *node = &kh_val(read_hash, k);
            Read **r = (Read **)node->data;
            for (readi = 0; readi < node->size; ++readi) {
                char **v = (char **)r[readi]->var_list->data;
                size_t nvariants = r[readi]->var_list->size;
 
                int multiallele;
                if (nvariants > 1) { // check if only multi-allelic variants at the same position & no hope of differentiating ref vs alt
                    multiallele = 1;
                    int prev_pos = -1;
                    for (i = 0; i < nvariants; ++i) { 
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
                        fprintf(stdout, "%s\tMUL\t%s\t%d\t%f\t%f\t%f\t", r[readi]->name, r[readi]->chr, r[readi]->pos, r[readi]->prgu, r[readi]->prgv, r[readi]->pout);
                        for (i = 0; i < nvariants; ++i) { fprintf(stdout, "%s;", v[i]); } fprintf(stdout, "\n");
                        continue;
                    }
                }

                // if any variant in list is not in EAGLE output, it suggests variants are from a different phase if alt wins
                multiallele = 0;
                for (i = 0; i < nvariants; ++i) { 
                    khiter_t k = kh_get(vh, var_hash, v[i]);
                    if (k != kh_end(var_hash)) {
                        Vector *var_hash_node = &kh_val(var_hash, k);
                        char **vv = (char **)var_hash_node->data;                                                                                                                          
                        for (j = 0; j < var_hash_node->size; ++j) {
                            if (strcmp(v[i], vv[j]) == 0) break;
                        }
                        if (j == var_hash_node->size) multiallele = 1;
                    }
                    else { 
                        multiallele = 1; 
                    }
                }

                if (r[readi]->prgu > r[readi]->prgv && r[readi]->prgu - r[readi]->prgv >= 0.69) { // ref wins
                    r[readi]->index = 0;
                    fprintf(stdout, "%s\tREF\t%s\t%d\t%f\t%f\t%f\t", r[readi]->name, r[readi]->chr, r[readi]->pos, r[readi]->prgu, r[readi]->prgv, r[readi]->pout);
                }
                else if (r[readi]->prgv > r[readi]->prgu && r[readi]->prgv - r[readi]->prgu >= 0.69) { // alt wins
                    if (multiallele) { // EAGLE outputs the set with highest likelihood ratio, i.e. most different from reference, leaving the "reference-like-allele"
                        r[readi]->index = 2;
                        fprintf(stdout, "%s\tRLA\t%s\t%d\t%f\t%f\t%f\t", r[readi]->name, r[readi]->chr, r[readi]->pos, r[readi]->prgu, r[readi]->prgv, r[readi]->pout);
                    }
                    else {
                        r[readi]->index = 1;
                        fprintf(stdout, "%s\tALT\t%s\t%d\t%f\t%f\t%f\t", r[readi]->name, r[readi]->chr, r[readi]->pos, r[readi]->prgu, r[readi]->prgv, r[readi]->pout);
                    }
                }
                else { // unknown
                    r[readi]->index = 4;
                    fprintf(stdout, "%s\tUNK\t%s\t%d\t%f\t%f\t%f\t", r[readi]->name, r[readi]->chr, r[readi]->pos, r[readi]->prgu, r[readi]->prgv, r[readi]->pout);
                }
                for (i = 0; i < nvariants; ++i) { fprintf(stdout, "%s;", v[i]); } fprintf(stdout, "\n");
            }
        }
    }
    fflush(stdout);
    print_status("# Reads Classified:\t%s", asctime(time_info));

    if (!listonly) {
        process_bam(bam_file, output_prefix);
        print_status("# BAM Processed:\t%s", asctime(time_info));
    }
}

static int type2ind(char *type) {
    if (strcmp("REF", type) == 0) return 0;
    else if (strcmp("ALT", type) == 0) return 1;
    else if (strcmp("RLA", type) == 0) return 2;
    else if (strcmp("MUL", type) == 0) return 3;
    else if (strcmp("UNK", type) == 0) return 4;
    return -1;
}

static void process_list(const char *filename, const char *bam_file, const char *output_prefix) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) { exit_err("failed to open file %s\n", filename); }

    Vector *ref = vector_create(64, VOID_T); // reference
    Vector *alt = vector_create(64, VOID_T); // alternative
    Vector *mul = vector_create(64, VOID_T); // multi-allelic that are undifferentiateable
    Vector *unk = vector_create(64, VOID_T); // unknown, ambiguous with equal likelihoods for reference and alternative

    int nreads = 0;
    char *line = NULL;
    ssize_t read_file = 0;
    size_t line_length = 0;
    while ((read_file = getline(&line, &line_length, file)) != -1) {
        if (line_length <= 0 || line[strspn(line, " \t\v\r\n")] == '\0') continue; // blank line
        if (line[0] == '#') continue;

        double prgu, prgv;
        char type[line_length], name[line_length];
        int t = sscanf(line, "%s %s %*[^\t] %*[^\t] %lf %lf %*[^\n]", name, type, &prgu, &prgv);
        if (t < 4) { exit_err("bad fields in read classified list file\n"); }

        khiter_t k = kh_get(rh, read_hash, name);
        if (k != kh_end(read_hash)) {
            size_t i;
            Vector *node = &kh_val(read_hash, k);
            Read **r = (Read **)node->data;                                                                                                                          
            for (i = 0; i < node->size; ++i) {
                if (strcmp(r[i]->name, name) == 0) {
                    if ((prgu + prgv) > (r[i]->prgu + r[i]->prgv)) {
                        r[i]->prgu = prgu;
                        r[i]->prgv = prgv;
                        r[i]->index = type2ind(type);
                        break;
                    }
                }
            }
            if (i == node->size) { exit_err("failed to find %s in hash key %d\n", name, k); }
        }
        else {
            Read *r = read_create(name, 0, "", 0);
            r->prgu = prgu;
            r->prgv = prgv;
            r->index = type2ind(type);

            int absent;
            khiter_t k = kh_put(rh, read_hash, r->name, &absent);
            Vector *node = &kh_val(read_hash, k);
            if (absent) vector_init(node, 8, READ_T);
            vector_add(node, r);
            ++nreads;
        }
    }
    free(line); line = NULL;
    fclose(file);
    print_status("# Classified list: %s\t%i reads\t%s", filename, nreads, asctime(time_info));

    process_bam(bam_file, output_prefix);
    print_status("# BAM Processed:\t%s", asctime(time_info));

    vector_destroy(ref); free(ref); ref = NULL;
    vector_destroy(alt); free(alt); alt = NULL;
    vector_destroy(mul); free(mul); mul = NULL;
    vector_destroy(unk); free(unk); unk = NULL;
}

static void print_usage() {
    printf("\n");
    printf("Usage: eagle-rc [options] eagle.out.txt eagle.readinfo.txt > classified_reads.list\n\n");
    printf("*  EAGLE with runtime options --omega=1e-40 --mvh --verbose\n");
    printf("*  ex) eagle -t 2 -v var.vcf -a align.bam -r ref.fa --omega=1.0e-40 --mvh --pao --isc --verbose 1> out.txt  2> readinfo.txt\n\n");
    printf("Options:\n");
    printf("  -o --out=     String   prefix for output BAM files\n");
    printf("  -a --bam=     FILE     alignment data BAM file corresponding to EAGLE output to be grouped into classes\n");
    printf("     --listonly          print classified read list only (stdout) without processing BAM file\n");
    printf("     --readlist          read from classified read list file instead of EAGLE outputs and proccess BAM file\n");
    printf("     --refonly           write REF classified reads only when processing BAM file\n");
}

int main(int argc, char **argv) {
    /* Command line parameters defaults */
    char *readinfo_file = NULL;
    char *var_file = NULL;
    char *bam_file = NULL;
    char *output_prefix = NULL;
    listonly = 0;
    readlist = 0;
    refonly = 0;

    static struct option long_options[] = {
        {"var", required_argument, NULL, 'v'},
        {"bam", required_argument, NULL, 'a'},
        {"out", optional_argument, NULL, 'o'},
        {"listonly", no_argument, &listonly, 1},
        {"readlist", no_argument, &readlist, 1},
        {"refonly", no_argument, &refonly, 1},
        {0, 0, 0, 0}
    };

    int opt = 0;
    while ((opt = getopt_long(argc, argv, "v:a:o:", long_options, &opt)) != -1) {
        switch (opt) {
            case 0: 
                //if (long_options[option_index].flag != 0) break;
                break;
            case 'v': var_file = optarg; break;
            case 'a': bam_file = optarg; break;
            case 'o': output_prefix = optarg; break;
            default: exit_usage("Bad options");
        }
    }
    if (optind > argc) { exit_usage("Bad program call"); }

    if (!listonly && bam_file == NULL) { exit_usage("Missing BAM file!"); } 

    /* Start processing data */
    clock_t tic = clock();

    var_hash = kh_init(vh);
    read_hash = kh_init(rh);

    if (!readlist) {
        var_file = argv[optind++];
        if (var_file == NULL) { exit_usage("Missing EAGLE output file!"); } 

        readinfo_file = argv[optind];
        if (readinfo_file == NULL) { exit_usage("Missing EAGLE read info file!"); } 

        int nvars = var_read(var_file);
        print_status("# Read EAGLE: %s\t%i entries\t%s", var_file, nvars, asctime(time_info));

        int nreads = readinfo_read(readinfo_file);
        print_status("# Read EAGLE: %s\t%i reads\t%s", readinfo_file, nreads, asctime(time_info));

        classify_reads(bam_file, output_prefix);
    }
    else {
        process_list(argv[optind], bam_file, output_prefix);
    }

    khiter_t k;
    for (k = kh_begin(var_hash); k != kh_end(var_hash); ++k) {
        if (kh_exist(var_hash, k)) vector_destroy(&kh_val(var_hash, k));
    }
    kh_destroy(vh, var_hash);

    for (k = kh_begin(read_hash); k != kh_end(read_hash); ++k) {
        if (kh_exist(read_hash, k)) vector_destroy(&kh_val(read_hash, k));
    }
    kh_destroy(rh, read_hash);

    clock_t toc = clock();
    print_status("# CPU time (hr):\t%f\n", (double)(toc - tic) / CLOCKS_PER_SEC / 3600);

    return 0;
}
