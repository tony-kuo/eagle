/*
Utility program that classifies reads based on EAGLE calculated likelihoods

Copyright 2016 Tony Kuo
This program is distributed under the terms of the GNU General Public License
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include "readclassify.h"
#include "htslib/sam.h"
#include "htslib/khash.h"

/* Precalculated log values */
#define M_1_LOG10E (1.0/M_LOG10E)
#define M_1_LN10 (1.0/M_LN10)

/* Time info */
static time_t now; 
static struct tm *time_info; 
#define print_status(M, ...) time(&now); time_info = localtime(&now); fprintf(stderr, M, ##__VA_ARGS__);

KHASH_MAP_INIT_STR(rh, Vector) // hashmap: string key, vector value
static khash_t(rh) *read_hash; // pointer to hashmap

static inline int has_numbers(const char *str) {
    while (*str != '\0') {
        if (isdigit(*str++) == 1) return 1;
    }
    return 0;
}

static inline void str_resize(char **str, size_t size) {
    char *p = realloc(*str, size * sizeof *str);
    if (p == NULL) { exit_err("failed to realloc in str_resize\n"); }
    else { *str = p; }
}

static inline double log_add_exp(double a, double b) {
    double max_exp = a > b ? a : b;
    return log(exp(a - max_exp) + exp(b - max_exp)) + max_exp;
}

static int nat_sort_cmp(const void *a, const void *b, enum type var_type) {
    char *str1, *str2;
    switch (var_type) {
        case VARIANT_T: {
            Variant *c1 = *(Variant **)a;
            Variant *c2 = *(Variant **)b;
            if (strcasecmp(c1->chr, c2->chr) == 0) return (c1->pos > c2->pos) - (c1->pos < c2->pos);
            str1 = strdup(c1->chr);
            str2 = strdup(c2->chr);
            c1 = NULL;
            c2 = NULL;
            break;
        }
        case STR_T: {
            str1 = strdup((char *)a);
            str2 = strdup((char *)b);
            break;
        }
        default: // vector i.e. VOID_T
            str1 = strdup(*(char **)a);
            str2 = strdup(*(char **)b);
            break;
    }
    char *s1 = str1;
    char *s2 = str2;
    int cmp = 0;
    while (cmp == 0 && *s1 != '\0' && *s2 != '\0') {
        if (isspace(*s1) && isspace(*s2)) { // ignore whitespace
            s1 += 1;
            s2 += 1;
        }
        else if ((isalpha(*s1) && isalpha(*s2)) || (ispunct(*s1) && ispunct(*s2))) { // compare alphabet and punctuation
            *s1 = tolower(*s1);
            *s2 = tolower(*s2);
            cmp = (*s1 > *s2) - (*s1 < *s2);
            s1 += 1;
            s2 += 1;
        }
        else { // compare digits
            int i1, i2, n1, n2, t1, t2;
            t1 = sscanf(s1, "%d%n", &i1, &n1);
            if (t1 == 0) t1 = sscanf(s1, "%*[^0123456789]%d%n", &i1, &n1);
            t2 = sscanf(s2, "%d%n", &i2, &n2);
            if (t2 == 0) t2 = sscanf(s2, "%*[^0123456789]%d%n", &i2, &n2);

            if (t1 < 1 || t2 < 1) { // one string has no digits
                cmp = strcmp(s1, s2);
            }
            else {
                cmp = (i1 > i2) - (i1 < i2);
                if (cmp == 0) { // first set of digits are equal, check further
                    s1 += n1;
                    s2 += n2;
                }
            }
        }
    }
    free(str1); str1 = NULL;
    free(str2); str2 = NULL;
    return cmp;
}

static int nat_sort_vector(const void *a, const void *b) {
    return nat_sort_cmp(a, b, VOID_T);
}

static int nat_sort_var(const void *a, const void *b) {
    return nat_sort_cmp(a, b, VARIANT_T);
}

static inline int str_find(const Vector *a, const void *s) {
    int i = 0;
    int j = a->size - 1;
    int n = (i + j) / 2;
    while (i <= j) {
        int v = nat_sort_cmp(s, a->data[n], STR_T);
        if (v == 0) return n;
        if (v > 0) i = n + 1;
        else j = n - 1;
        n = (i + j) / 2;
    }
    return -1;
}

void variant_destroy(Variant *v) {
    if (v == NULL) return;
    v->pos = 0;
    free(v->chr); v->chr = NULL;      
    free(v->ref); v->ref = NULL;
    free(v->alt); v->alt = NULL;
}

void read_destroy(Read *r) {
    if (r == NULL) return;
    r->pos = 0;
    r->prgu = r->prgv = r->pout = 0;
    free(r->name); r->name = NULL;
    free(r->chr); r->chr = NULL;
}

void vector_init(Vector *a, size_t initial_size, enum type var_type) {
    a->size = 0;
    a->type = var_type;
    a->capacity = initial_size;
    a->data = malloc(initial_size * sizeof (void *));
}

void vector_add(Vector *a, void *entry) {
    if (a->size >= a->capacity) {
        a->capacity *= 2;
        void **p = realloc(a->data, a->capacity * sizeof (void *));
        if (p == NULL) { exit_err("failed to realloc in vector_add\n"); }
        else { a->data = p; }
    }
    a->data[a->size++] = entry;
}

void vector_del(Vector *a, int i) {
    a->data[i] = NULL;
    if (i == --a->size) return;
    memmove(&(a->data[i]), &(a->data[i + 1]), (a->size - i) * sizeof (void *));
    a->data[a->size] = NULL;
}

void *vector_pop(Vector *a) {
    if (a->size <= 0) return NULL;
    void *entry = a->data[--a->size];
    a->data[a->size] = NULL;
    return entry;
}

Vector *vector_create(size_t initial_size, enum type var_type) {
    Vector *a = malloc(sizeof (Vector));
    vector_init(a, initial_size, var_type);
    return a;
}

Vector *vector_dup(Vector *a) {
    Vector *v = vector_create(a->capacity, a->type);
    v->size = a->size;
    memcpy(&(v->data[0]), &(a->data[0]), a->size * sizeof (void *));
    return v;
}

void vector_free(Vector *a) {
    a->size = 0;
    free(a->data); a->data = NULL;
}

void vector_destroy(Vector *a) {
    size_t i;
    enum type var_type = a->type;
    for (i = 0; i < a->size; ++i) {
        switch (var_type) {
            case VARIANT_T:
                variant_destroy((Variant *)a->data[i]);
                break;
            case READ_T:
                read_destroy((Read *)a->data[i]);
                break;
            default:
                break;
        }
        free(a->data[i]); a->data[i] = NULL;
    }
    a->size = a->capacity = 0;
    free(a->data); a->data = NULL;
}

static inline int variant_find(const Vector *a, const Variant *v, const int sorted) {
    int i = 0;
    if (sorted) {
        int j = a->size - 1;
        int n = (i + j) / 2;
        while (i <= j) {
            Variant *curr = (Variant *)a->data[n];
            if (strcmp(v->chr, curr->chr) == 0 && strcmp(v->ref, curr->ref) == 0 && strcmp(v->alt, curr->alt) == 0 && v->pos == curr->pos) return n;
            if (v->pos > curr->pos) i = n + 1;
            else j = n - 1;
            n = (i + j) / 2;
        }
        return -1;
    }
    else {
        for (i = 0; i < a->size; ++i) {
            Variant *curr = (Variant *)a->data[i];
            if (strcmp(v->chr, curr->chr) == 0 && strcmp(v->ref, curr->ref) == 0 && strcmp(v->alt, curr->alt) == 0 && v->pos == curr->pos) return i;
        }
        return -1;
    }
}

static void add_to_varlist(Vector *var_list, char *set) {
    int pos;
    int line_length = strlen(set);
    char chr[line_length], ref[line_length], alt[line_length];

    int n;
    char *s;
    for (s = set + 1; sscanf(s, "%[^,],%d,%[^,],%[^;];%n", chr, &pos, ref, alt, &n) == 4; s += n) { // scan variant set
        Variant *v = malloc(sizeof (Variant));
        v->chr = strdup(chr);
        v->pos = pos;
        v->ref = strdup(ref);
        v->alt = strdup(alt);
        if (variant_find(var_list, v, 0) == -1) vector_add(var_list, v);
        if (*(s + n) == ']') break;
    }
}

static Vector *var_read(const char* filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) { exit_err("failed to open file %s\n", filename); }

    Vector *var_list = vector_create(64, VARIANT_T);

    char *line = NULL;
    ssize_t read_file = 0;
    size_t line_length = 0;
    while ((read_file = getline(&line, &line_length, file)) != -1) {
        if (line_length <= 0 || line[strspn(line, " \t\v\r\n")] == '\0') continue; // blank line
        if (line[0] == '#') continue;

        int pos;
        char chr[line_length], ref[line_length], alt[line_length], set[line_length];
        int t = sscanf(line, "%s %d %s %s %*[^\t] %*[^\t] %*[^\t] %*[^\t] %*[^\t] %s", chr, &pos, ref, alt, set);
        if (t < 4) { exit_err("bad fields in EAGLE output file\n"); }

        if (strcmp(set, "[]") == 0) {
            Variant *v = malloc(sizeof (Variant));
            v->chr = strdup(chr);
            v->pos = pos;
            v->ref = strdup(ref);
            v->alt = strdup(alt);
            vector_add(var_list, v);
        }
        else {
            add_to_varlist(var_list, set);
        }
    }
    free(line); line = NULL;
    fclose(file);
    qsort(var_list->data, var_list->size, sizeof (void *), nat_sort_var);
    return var_list;
}

static void readinfo_read(const char* filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) { exit_err("failed to open file %s\n", filename); }

    size_t i;
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
        //fprintf(stderr, "%d\t%s %s %d %f %f %f %s\n", t, name, chr, pos, prgu, prgv, pout, set);

        khiter_t k = kh_get(rh, read_hash, name);
        if (k != kh_end(read_hash)) {
            Vector *node = &kh_val(read_hash, k);
            Read **r = (Read **)node->data;                                                                                                                          
            int found = 0;
            for (i = 0; i < node->size; ++i) {
                if (strcmp(r[i]->name, name) == 0) {
                    r[i]->prgu = log_add_exp(r[i]->prgu, prgu);
                    r[i]->prgv = log_add_exp(r[i]->prgv, prgv);
                    add_to_varlist(r[i]->var_list, set);
                    found = 1;
                    break;
                }
            }
            if (!found) { exit_err("failed to find sequence %s in hash key %d\n", name, k); }
        }
        else {
            Read *r = malloc(sizeof (Read));
            r->name = strdup(name);
            r->chr = strdup(chr);
            r->pos = pos;
            r->prgu = prgu;
            r->prgv = prgv;
            r->pout = pout;
            r->var_list = vector_create(8, VARIANT_T);
            add_to_varlist(r->var_list, set);

            int absent;
            k = kh_put(rh, read_hash, r->name, &absent);
            Vector *node = &kh_val(read_hash, k);
            if(absent) vector_init(node, 8, READ_T);
            vector_add(node, r);
        }
    }
    free(line); line = NULL;
    fclose(file);
}

void classify_reads(Vector *var_list, const char* bam_file, const char *output_prefix) {
    Vector *ref = vector_create(64, VOID_T); // reference
    Vector *alt = vector_create(64, VOID_T); // alternative
    Vector *mul = vector_create(64, VOID_T); // multi-allelic that are undifferentiateable
    Vector *unk = vector_create(64, VOID_T); // unknown, ambiguous with equal likelihoods for reference and alternative

    size_t readi, i;
	khiter_t k;
    for (k = kh_begin(read_hash); k != kh_end(read_hash); ++k) {
		if (kh_exist(read_hash, k)) {
            Vector *node = &kh_val(read_hash, k);
            Read **r = (Read **)node->data;
            for (readi = 0; readi < node->size; ++readi) {
                Variant **v = (Variant **)r[readi]->var_list->data;
                size_t nvariants = r[readi]->var_list->size;
 
                int multiallele;
                if (nvariants > 1) { // check if only multi-allelic variants at the same position & no hope of differentiating ref vs alt
                    multiallele = 1;
                    for (i = 1; i < nvariants; ++i) { 
                        if (v[i]->pos != v[i - 1]->pos) multiallele = 0;
                    }
                    if (multiallele) {
                        vector_add(mul, r[readi]->name);
                        fprintf(stdout, "MUL=\t%s\t%s\t%d\t%f\t%f\t%f\t", r[readi]->name, r[readi]->chr, r[readi]->pos, r[readi]->prgu, r[readi]->prgv, r[readi]->pout);
                        for (i = 0; i < nvariants; ++i) { fprintf(stdout, "%s,%d,%s,%s;", v[i]->chr, v[i]->pos, v[i]->ref, v[i]->alt); } fprintf(stdout, "\n");
                        continue;
                    }
                }
                multiallele = 0;
                for (i = 0; i < nvariants; ++i) { // check if variants are not in the EAGLE output, suggesting variants in a different phase, 
                    if (variant_find(var_list, v[i], 1) == -1) multiallele = 1;
                }
                if (multiallele) { // EAGLE outputs the set with highest likelihood ratio, i.e. most different from reference, leaving the "reference-like-allele"
                    vector_add(ref, r[readi]->name);
                    fprintf(stdout, "RLA=\t%s\t%s\t%d\t%f\t%f\t%f\t", r[readi]->name, r[readi]->chr, r[readi]->pos, r[readi]->prgu, r[readi]->prgv, r[readi]->pout);
                    for (i = 0; i < nvariants; ++i) { fprintf(stdout, "%s,%d,%s,%s;", v[i]->chr, v[i]->pos, v[i]->ref, v[i]->alt); } fprintf(stdout, "\n");
                    continue;
                }
                if (r[readi]->prgu > r[readi]->prgv && r[readi]->prgu - r[readi]->prgv >= 0.69) { // ref wins
                    vector_add(ref, r[readi]->name);
                    fprintf(stdout, "REF=\t%s\t%s\t%d\t%f\t%f\t%f\t", r[readi]->name, r[readi]->chr, r[readi]->pos, r[readi]->prgu, r[readi]->prgv, r[readi]->pout);
                    for (i = 0; i < nvariants; ++i) { fprintf(stdout, "%s,%d,%s,%s;", v[i]->chr, v[i]->pos, v[i]->ref, v[i]->alt); } fprintf(stdout, "\n");
                }
                else if (r[readi]->prgv > r[readi]->prgu && r[readi]->prgv - r[readi]->prgu >= 0.69) { // alt wins
                    vector_add(alt, r[readi]->name);
                    fprintf(stdout, "ALT=\t%s\t%s\t%d\t%f\t%f\t%f\t", r[readi]->name, r[readi]->chr, r[readi]->pos, r[readi]->prgu, r[readi]->prgv, r[readi]->pout);
                    for (i = 0; i < nvariants; ++i) { fprintf(stdout, "%s,%d,%s,%s;", v[i]->chr, v[i]->pos, v[i]->ref, v[i]->alt); } fprintf(stdout, "\n");
                }
                else { // unknown
                    vector_add(unk, r[readi]->name);
                    fprintf(stdout, "UNK=\t%s\t%s\t%d\t%f\t%f\t%f\t", r[readi]->name, r[readi]->chr, r[readi]->pos, r[readi]->prgu, r[readi]->prgv, r[readi]->pout);
                    for (i = 0; i < nvariants; ++i) { fprintf(stdout, "%s,%d,%s,%s;", v[i]->chr, v[i]->pos, v[i]->ref, v[i]->alt); } fprintf(stdout, "\n");
                }
            }
        }
    }
    qsort(ref->data, ref->size, sizeof (void *), nat_sort_vector);
    qsort(alt->data, alt->size, sizeof (void *), nat_sort_vector);
    qsort(mul->data, mul->size, sizeof (void *), nat_sort_vector);
    qsort(unk->data, unk->size, sizeof (void *), nat_sort_vector);

    //for (i = 0; i < ref->size; ++i) { fprintf(stdout, "%s\n", (char *)ref->data[i]); }

    samFile *sam_in = sam_open(bam_file, "r"); // open bam file
    if (sam_in == NULL) { exit_err("failed to open BAM file %s\n", bam_file); }
    bam_hdr_t *bam_header = sam_hdr_read(sam_in); // bam header
    if (bam_header == 0) { exit_err("bad header %s\n", bam_file); }

    /* Split input bam files into respectively categories' output bam files */
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
        char *s, token[strlen(bam_flag2str(aln->core.flag)) + 1];
        for (s = bam_flag2str(aln->core.flag); sscanf(s, "%[^,]%n", token, &n) == 1; s += n + 1) {
            if (strcmp("UNMAP", token) == 0) is_unmap = 1;
            else if (strcmp("SECONDARY", token) == 0 || strcmp("SUPPLEMENTARY", token) == 0) is_secondary = 1;
            if (*(s + n) != ',') break;
        }
        if (is_unmap || is_secondary) continue;

        /* Write reads to appropriate file */
        char *name = (char *)aln->data;
        if (str_find(ref, name) >= 0) {
            out = ref_out;
        }
        else if (str_find(alt, name) >= 0) {
            out = alt_out;
        }
        else if (str_find(mul, name) >= 0) {
            out = alt_out;
        }
        else {
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

static void print_usage() {
    printf("\n");
    printf("Usage: readclassify [options] eagle.out.txt eagle.readinfo.txt > classified_reads.list\n\n");
    printf("*  EAGLE with runtime options --omega=1e-40 --mvh -d -1\n");
    printf("*  ex) eagle -t 2 -v var.vcf -a align.bam -r ref.fa --omega=1.0e-40 --mvh --pao --isc -d -1 1> out.txt  2> readinfo.txt\n\n");
    printf("Options:\n");
    printf("  -o --out=    String   output prefix for sam files\n");
    printf("  -a --bam=    FILE     alignment data bam files corresponding to EAGLE output\n");
}

int main(int argc, char **argv) {
    /* Command line parameters defaults */
    char *readinfo_file = NULL;
    char *var_file = NULL;
    char *bam_file = NULL;
    char *output_prefix = NULL;

    static struct option long_options[] = {
        {"var", required_argument, NULL, 'v'},
        {"bam", required_argument, NULL, 'a'},
        {"out", optional_argument, NULL, 'o'},
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


    var_file = argv[optind++];
    if (var_file == NULL) { exit_usage("Missing EAGLE output file!"); } 

    readinfo_file = argv[optind];
    if (readinfo_file == NULL) { exit_usage("Missing EAGLE read info file!"); } 

    /* Start processing data */
    clock_t tic = clock();
    Vector *var_list = var_read(var_file);
    print_status("# Read EAGLE: %s\t%i entries\t%s", var_file, (int)var_list->size, asctime(time_info));

    read_hash = kh_init(rh);
    readinfo_read(readinfo_file);

    classify_reads(var_list, bam_file, output_prefix);

    clock_t toc = clock();
    print_status("# Done:\t%s\t%s", readinfo_file, asctime(time_info));
    print_status("# CPU time (hr):\t%f\n", (double)(toc - tic) / CLOCKS_PER_SEC / 3600);

	khiter_t k;
    for (k = kh_begin(read_hash); k != kh_end(read_hash); ++k) {
		if (kh_exist(read_hash, k)) vector_destroy(&kh_val(read_hash, k));
    }
	kh_destroy(rh, read_hash);
    vector_destroy(var_list); free(var_list); var_list = NULL;
    return 0;
}
