/*
EAGLE: explicit alternative genome likelihood evaluator
Given the sequencing data and candidate variant, explicitly test 
the alternative hypothesis against the reference hypothesis

Copyright 2016 Tony Kuo
This program is distributed under the terms of the GNU General Public License
*/

#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <pthread.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/khash.h"
#include "util.h"
#include "vector.h"

/* Constants */
#define ALPHA 1.3     // Factor to account for longer read lengths lowering the probability a sequence matching an outside paralogous source
#define NT_CODES 17   // Size of nucleotide code table

/* Precalculated log values */
#define M_1_LOG10E (1.0/M_LOG10E)
#define M_1_LN10 (1.0/M_LN10)
#define LG50 (log(0.5))
#define LG10 (log(0.1))
#define LG90 (log(0.9))
#define LGALPHA (log(ALPHA))

/* Command line arguments */
static char *vcf_file;
static char *bam_file;
static char *fa_file;
static char *out_file;
static int nthread;
static int sharedr;
static int distlim;
static int maxdist;
static int maxh;
static int mvh;
static int pao;
static int isc;
static int nodup;
static int splice;
static int dp;
static int verbose;
static int debug;
static double hetbias;
static double omega, lgomega;
static int match, mismatch, gap_op, gap_ex;
static double ref_prior, alt_prior, het_prior;

/* Time info */
static time_t now; 
static struct tm *time_info; 
#define print_status(M, ...) time(&now); time_info = localtime(&now); fprintf(stderr, M, ##__VA_ARGS__);

/* Mapping table */
static int seqnt_map[26];

/* Fastq quality to probability table */
static double p_match[50], p_mismatch[50];

KHASH_MAP_INIT_STR(rsh, Vector)   // hashmap: string key, vector value
static khash_t(rsh) *refseq_hash; // pointer to hashmap
static pthread_mutex_t refseq_lock; 

static Vector *vcf_read(FILE *file) {
    Vector *var_list = vector_create(64, VARIANT_T);

    char *line = NULL;
    ssize_t read_file = 0;
    size_t line_length = 0;
    while ((read_file = getline(&line, &line_length, file)) != -1) {
        if (line_length <= 0 || line[strspn(line, " \t\v\r\n")] == '\0') continue; // blank line
        if (line[0] == '#') continue;

        int pos;
        char chr[line_length], ref[line_length], alt[line_length];
        int t = sscanf(line, "%s %d %*[^\t] %s %s", chr, &pos, ref, alt);
        if (t < 4 || has_numbers(ref) || has_numbers(alt)) { exit_err("bad fields in VCF file\n"); }

        int n1, n2;
        char *s1, *s2, ref_token[strlen(ref) + 1], alt_token[strlen(alt) + 1];
        for (s1 = ref; sscanf(s1, "%[^, ]%n", ref_token, &n1) == 1 || sscanf(s1, "%[-]%n", ref_token, &n1) == 1; s1 += n1 + 1) { // heterogenenous non-reference (comma-delimited) as separate entries
            for (s2 = alt; sscanf(s2, "%[^, ]%n", alt_token, &n2) == 1 || sscanf(s2, "%[-]%n", alt_token, &n2) == 1; s2 += n2 + 1) {
                Variant *v = variant_create(chr, pos, ref_token, alt_token);
                vector_add(var_list, v);
                if (*(s2 + n2) != ',') break;
            }
            if (*(s1 + n1) != ',') break;
        }
    }
    free(line); line = NULL;
    fclose(file);
    qsort(var_list->data, var_list->size, sizeof (void *), nat_sort_variant);
    return var_list;
}

void fasta_read(const char *fa_file) {
    faidx_t *fai = fai_load(fa_file);
    if (fai == NULL) { 
        errno = fai_build(fa_file);
        if (errno == 0) {
            fai = fai_load(fa_file);
        }
        else {
            exit_err("failed to build and open FA index %s\n", fa_file);
        }
    }

    char *filename = malloc((strlen(fa_file) + 5) * sizeof *filename);
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

        Fasta *f = malloc(sizeof (Fasta));
        f->name = strdup(name);
        f->seq = fai_fetch(fai, name, &f->seq_length);
        char *s;
        for (s = f->seq; *s != '\0'; ++s) *s = toupper(*s);

        int absent;
        khiter_t k = kh_put(rsh, refseq_hash, f->name, &absent);
        Vector *node = &kh_val(refseq_hash, k); // point to the bucket associated to k
        if (absent) vector_init(node, 8, FASTA_T);
        vector_add(node, f);
    }
    free(line); line = NULL;
    fclose(file);
    fai_destroy(fai);
    print_status("# Read reference genome: %s\t%s", fa_file, asctime(time_info));
}

static int bam_fetch_last(const char *bam_file, const char *chr, const int pos1, const int pos2) {
    /* Reads in variant j = i + 1 region coordinates */
    samFile *sam_in = sam_open(bam_file, "r"); // open bam file
    if (sam_in == NULL) { exit_err("failed to open BAM file %s\n", bam_file); }
    bam_hdr_t *bam_header = sam_hdr_read(sam_in); // bam header
    if (bam_header == 0) { exit_err("bad header %s\n", bam_file); }
    hts_idx_t *bam_idx = sam_index_load(sam_in, bam_file); // bam index
    if (bam_idx == NULL) { exit_err("failed to open BAM index %s\n", bam_file); }

    int last = -1;
    int tid = bam_name2id(bam_header, chr);
    hts_itr_t *iter = sam_itr_queryi(bam_idx, tid, pos1-1, pos2); // read iterator
    if (iter != NULL) {
        bam1_t *aln = bam_init1(); // initialize an alignment
        while (sam_itr_next(sam_in, iter, aln) >= 0) {
            last = aln->core.pos + aln->core.l_qseq;
        }
        bam_destroy1(aln);
    }
    hts_itr_destroy(iter);
    hts_idx_destroy(bam_idx);
    bam_hdr_destroy(bam_header);
    sam_close(sam_in);
    return(last);
}

static Vector *bam_fetch(const char *bam_file, const char *chr, const int pos1, const int pos2) {
    /* Reads in variant region coordinates */
    Vector *read_list = vector_create(64, READ_T);

    samFile *sam_in = sam_open(bam_file, "r"); // open bam file
    if (sam_in == NULL) { exit_err("failed to open BAM file %s\n", bam_file); }
    bam_hdr_t *bam_header = sam_hdr_read(sam_in); // bam header
    if (bam_header == 0) { exit_err("bad header %s\n", bam_file); }
    hts_idx_t *bam_idx = sam_index_load(sam_in, bam_file); // bam index
    if (bam_idx == NULL) { exit_err("failed to open BAM index %s\n", bam_file); }

    int tid = bam_name2id(bam_header, chr);
    hts_itr_t *iter = sam_itr_queryi(bam_idx, tid, pos1-1, pos2); // read iterator
    if (iter != NULL) {
        bam1_t *aln = bam_init1(); // initialize an alignment
        while (sam_itr_next(sam_in, iter, aln) >= 0) {
            size_t i, j;
            Read *read = read_create((char *)aln->data, aln->core.tid, bam_header->target_name[aln->core.tid], aln->core.pos);

            int saw_M = 0;
            int s_offset = 0; // offset for softclip at start
            int e_offset = 0; // offset for softclip at end

            uint32_t *cigar = bam_get_cigar(aln);
            read->n_cigar = aln->core.n_cigar;
            read->cigar_oplen = malloc(read->n_cigar * sizeof read->cigar_oplen);
            read->cigar_opchr = malloc((read->n_cigar + 1) * sizeof read->cigar_opchr);
            read->splice_pos = malloc(read->n_cigar * sizeof read->splice_pos);
            read->splice_offset = malloc(read->n_cigar * sizeof read->splice_offset);

            j = 0;
            int splice_pos = 0; // track splice position in reads
            for (i = 0; i < read->n_cigar; ++i) {
                read->cigar_oplen[i] = bam_cigar_oplen(cigar[i]);
                read->cigar_opchr[i] = bam_cigar_opchr(cigar[i]);
                read->splice_pos[i] = 0;
                read->splice_offset[i] = 0;

                if (isc && read->cigar_opchr[i] == 'M') { 
                    saw_M = 1; 
                }
                else if (isc && saw_M == 0 && read->cigar_opchr[i] == 'S') { 
                    s_offset = read->cigar_oplen[i]; 
                }
                else if (isc && saw_M == 1 && read->cigar_opchr[i] == 'S') { 
                    e_offset = read->cigar_oplen[i]; 
                }

                if (splice && read->cigar_opchr[i] == 'N') {
                    read->splice_pos[j] = splice_pos - s_offset;
                    read->splice_offset[j] = read->cigar_oplen[i];
                    ++j;
                }
                else if (splice && (read->cigar_opchr[i] == 'M' || read->cigar_opchr[i] == 'I' || read->cigar_opchr[i] == 'S')) {
                    splice_pos += read->cigar_oplen[i];
                }
            }
            read->cigar_opchr[read->n_cigar] = '\0';
            read->inferred_length = bam_cigar2qlen(read->n_cigar, cigar);
            read->n_splice = j;

            read->pos += s_offset;
            read->length = aln->core.l_qseq - (s_offset + e_offset);
            read->qseq = malloc((read->length + 1) * sizeof read->qseq);
            read->qual = malloc(read->length  * sizeof read->qual);
            uint8_t *qual = bam_get_qual(aln);
            for (i = 0; i < read->length; ++i) {
                read->qseq[i] = toupper(seq_nt16_str[bam_seqi(bam_get_seq(aln), i + s_offset)]); // get nucleotide id and convert into IUPAC id.
                read->qual[i] = qual[i + s_offset];
            }
            read->qseq[read->length] = '\0';

            char *s;
            read->flag = NULL;
            s = bam_flag2str(aln->core.flag);
            if (s != NULL) read->flag = strdup(s);
            free(s); s = NULL;

            read->multimapXA = NULL;
            if (bam_aux_get(aln, "XA")) read->multimapXA = strdup(bam_aux2Z(bam_aux_get(aln, "XA")));

            read->multimapNH = 1;
            if (bam_aux_get(aln, "NH")) read->multimapNH = bam_aux2i(bam_aux_get(aln, "NH"));

            vector_add(read_list, read);
        }
        bam_destroy1(aln);
    }
    hts_itr_destroy(iter);
    hts_idx_destroy(bam_idx);
    bam_hdr_destroy(bam_header);
    sam_close(sam_in);
    return read_list;
}

static Fasta *refseq_fetch(char *name, const char *fa_file) {
    pthread_mutex_lock(&refseq_lock);
    size_t i;
	khiter_t k = kh_get(rsh, refseq_hash, name);
    if (k != kh_end(refseq_hash)) {
        Vector *node = &kh_val(refseq_hash, k);
        Fasta **f = (Fasta **)node->data;
        for (i = 0; i < node->size; ++i) {
            if (strcmp(f[i]->name, name) == 0) {
                pthread_mutex_unlock(&refseq_lock);
                return f[i];
            }
        }
        exit_err("failed to find %s in hash key %d\n", name, k);
    }

    faidx_t *fai = fai_load(fa_file);
    if (fai == NULL) { 
        errno = fai_build(fa_file);
        if (errno == 0) {
            fai = fai_load(fa_file);
        }
        else {
            exit_err("failed to build and open FA index %s\n", fa_file);
        }
    }
    if (!faidx_has_seq(fai, name)) { exit_err("failed to find %s in reference %s\n", name, fa_file); }

    Fasta *f = fasta_create(name);
    f->seq = fai_fetch(fai, name, &f->seq_length);
    char *s;
    for (s = f->seq; *s != '\0'; ++s) *s = toupper(*s);

    int absent;
    k = kh_put(rsh, refseq_hash, f->name, &absent);
    Vector *node = &kh_val(refseq_hash, k);
    if (absent) vector_init(node, 8, FASTA_T);
    vector_add(node, f);
    fai_destroy(fai);
    pthread_mutex_unlock(&refseq_lock);
    return f;
}

static char *construct_altseq(const char *refseq, int refseq_length, const Vector_Int *combo, Variant **var_data, int *altseq_length) {
    size_t i;
    int offset = 0;
    char *altseq = strdup(refseq);
    *altseq_length = refseq_length;
    for (i = 0; i < combo->size; ++i) {
        Variant *curr = var_data[combo->data[i]];
        size_t pos = curr->pos - 1 + offset;
        if (pos < 0 || pos > *altseq_length) { exit_err("Variant at %s:%d is out of bounds in reference\n", curr->chr, curr->pos); }

        char *var_ref, *var_alt;
        if (curr->ref[0] == '-') { // account for "-" variant representations
            ++pos;
            var_ref = "";
            var_alt = curr->alt;
        }
        else if (curr->alt[0] == '-') { // account for "-" variant representations
            var_ref = curr->ref;
            var_alt = "";
        }
        else {
            char *s1 = curr->ref;
            char *s2 = curr->alt;
            while (*s1 == *s2) { // account for and disregard common prefix in variant
                ++s1;
                ++s2;
                ++pos;
            }
            var_ref = s1;
            var_alt = s2;
        }
        size_t var_ref_length = strlen(var_ref);
        size_t var_alt_length = strlen(var_alt);
        int delta = var_alt_length - var_ref_length;
        offset += delta;
        if (delta == 0) { // snps, equal length haplotypes
            memcpy(altseq + pos, var_alt, var_alt_length * sizeof *var_alt);
        }
        else { // indels
            char *newalt = malloc((*altseq_length + delta + 1) * sizeof *newalt);
            memcpy(newalt, altseq, pos * sizeof *newalt);
            memcpy(newalt + pos, var_alt, var_alt_length * sizeof *newalt);
            memcpy(newalt + pos + var_alt_length, altseq + pos + var_ref_length, (*altseq_length - pos - var_ref_length) * sizeof *newalt);
            *altseq_length += delta;
            newalt[*altseq_length] = '\0';
            free(altseq); altseq = NULL;
            altseq = newalt;
        }
    }
    return altseq;
}

static inline int variant_find(const Vector_Int *a, int v) {
    int i = 0;
    int j = a->size - 1;
    int n = (i + j) / 2;
    while (i <= j) {
        if (v == a->data[n]) return n;
        if (v > a->data[n]) i = n + 1;
        else j = n - 1;
        n = (i + j) / 2;
    }
    return -1;
}

static inline void variant_print(char **output, const Vector *var_set, size_t i, int nreads, int not_alt_count, int has_alt_count, double total, double has_alt, double not_alt) {
    Variant **var_data = (Variant **)var_set->data;

    size_t n;
    char *token;
    double prob = (has_alt - total) * M_1_LN10;
    double odds = (has_alt - not_alt) * M_1_LN10;

    n = snprintf(NULL, 0, "%s\t%d\t%s\t%s\t%d\t%d\t%d\t%e\t%f\t", var_data[i]->chr, var_data[i]->pos, var_data[i]->ref, var_data[i]->alt, nreads, not_alt_count, has_alt_count, prob, odds) + 1;
    token = malloc(n * sizeof *token);
    snprintf(token, n, "%s\t%d\t%s\t%s\t%d\t%d\t%d\t%e\t%f\t", var_data[i]->chr, var_data[i]->pos, var_data[i]->ref, var_data[i]->alt, nreads, not_alt_count, has_alt_count, prob, odds);
    str_resize(output, strlen(*output) + n);
    strcat(*output, token);
    free(token); token = NULL;

    str_resize(output, strlen(*output) + 2);
    strcat(*output, "[");
    if (var_set->size > 1) {
        for (i = 0; i < var_set->size; ++i) {
            n = snprintf(NULL, 0, "%s,%d,%s,%s;", var_data[i]->chr, var_data[i]->pos, var_data[i]->ref, var_data[i]->alt) + 1;
            token = malloc(n * sizeof *token);
            snprintf(token, n, "%s,%d,%s,%s;", var_data[i]->chr, var_data[i]->pos, var_data[i]->ref, var_data[i]->alt);
            str_resize(output, strlen(*output) + n);
            strcat(*output, token);
            free(token); token = NULL;
        }
    }
    str_resize(output, strlen(*output) + 3);
    strcat(*output, "]\n");
}

static inline void set_prob_matrix(double *matrix, const char *seq, int read_length, const double *is_match, const double *no_match) {
    int i, b; // array[width * row + col] = value
    for (b = 0; b < read_length; ++b) {
        for (i = 0; i < NT_CODES; ++i) matrix[NT_CODES * b + i] = no_match[b];
        matrix[NT_CODES * b + seqnt_map[seq[b] - 'A']] = is_match[b];
        switch (seq[b]) {
        case 'A':
            matrix[NT_CODES * b + seqnt_map['M' - 'A']] = is_match[b];
            matrix[NT_CODES * b + seqnt_map['R' - 'A']] = is_match[b];
            matrix[NT_CODES * b + seqnt_map['V' - 'A']] = is_match[b];
            matrix[NT_CODES * b + seqnt_map['H' - 'A']] = is_match[b];
            matrix[NT_CODES * b + seqnt_map['D' - 'A']] = is_match[b];
            matrix[NT_CODES * b + seqnt_map['W' - 'A']] = is_match[b];
            matrix[NT_CODES * b + 9] = is_match[b]; // also W
            break;
        case 'T':
            matrix[NT_CODES * b + seqnt_map['K' - 'A']] = is_match[b];
            matrix[NT_CODES * b + seqnt_map['Y' - 'A']] = is_match[b];
            matrix[NT_CODES * b + seqnt_map['B' - 'A']] = is_match[b];
            matrix[NT_CODES * b + seqnt_map['H' - 'A']] = is_match[b];
            matrix[NT_CODES * b + seqnt_map['D' - 'A']] = is_match[b];
            matrix[NT_CODES * b + seqnt_map['W' - 'A']] = is_match[b];
            matrix[NT_CODES * b + 9] = is_match[b]; // also W
            break;
        case 'C':
            matrix[NT_CODES * b + seqnt_map['M' - 'A']] = is_match[b];
            matrix[NT_CODES * b + seqnt_map['Y' - 'A']] = is_match[b];
            matrix[NT_CODES * b + seqnt_map['B' - 'A']] = is_match[b];
            matrix[NT_CODES * b + seqnt_map['V' - 'A']] = is_match[b];
            matrix[NT_CODES * b + seqnt_map['H' - 'A']] = is_match[b];
            matrix[NT_CODES * b + seqnt_map['S' - 'A']] = is_match[b];
            matrix[NT_CODES * b + 10] = is_match[b]; // also S
            break;
        case 'G':
            matrix[NT_CODES * b + seqnt_map['K' - 'A']] = is_match[b];
            matrix[NT_CODES * b + seqnt_map['R' - 'A']] = is_match[b];
            matrix[NT_CODES * b + seqnt_map['B' - 'A']] = is_match[b];
            matrix[NT_CODES * b + seqnt_map['V' - 'A']] = is_match[b];
            matrix[NT_CODES * b + seqnt_map['D' - 'A']] = is_match[b];
            matrix[NT_CODES * b + seqnt_map['S' - 'A']] = is_match[b];
            matrix[NT_CODES * b + 10] = is_match[b]; // also S
            break;
        }
    }
}

static double smith_waterman_gotoh(const double *matrix, int read_length, const char *seq, int seq_length, int start) { /* short in long version */
    size_t i, j;

    double prev[read_length + 1], curr[read_length + 1];
    double a_gap_prev[read_length + 1], a_gap_curr[read_length + 1];
    double b_gap_prev[read_length + 1], b_gap_curr[read_length + 1];

    memset(a_gap_prev, -1e6, sizeof(a_gap_prev));
    memset(b_gap_prev, -1e6, sizeof(b_gap_prev));

    prev[0] = 0;
    prev[1] = 0 - gap_op;
    for (j = 2; j < read_length + 1; ++j) prev[j] = prev[j - 1] - gap_ex;

    int n = start + seq_length;
    double max_score = 0;
    for (i = start; i < n; ++i) {
        curr[0] = 0;
        a_gap_curr[0] = 0;
        b_gap_curr[0] = 0;
        double row_max = 0;
        for (j = 1; j < read_length + 1; ++j) {
            int t = seq[i] - 'A';
            if (t < 0 || t >= 26) { exit_err("Character %c at pos %d (%d) not in valid alphabet\n", seq[i], (int)i, seq_length); }

            double upleft = prev[j - 1] + matrix[NT_CODES * (j - 1) + seqnt_map[t]];

            double open = curr[j - 1] - gap_op;
            double extend = a_gap_curr[j - 1] - gap_ex;
            a_gap_curr[j] = (open >= extend) ? open : extend;

            open = prev[j] - gap_op;
            extend = b_gap_prev[j] - gap_ex;
            b_gap_curr[j] = (open >= extend) ? open : extend;

            curr[j] = upleft;
            if (a_gap_curr[j] >= curr[j]) curr[j] = a_gap_curr[j];
            if (b_gap_curr[j] >= curr[j]) curr[j] = b_gap_curr[j];
            if (curr[j] > row_max) row_max = curr[j];
        }
        if (row_max > max_score) max_score = row_max;

        memcpy(prev, curr, sizeof(prev));
        memcpy(a_gap_prev, a_gap_curr, sizeof(a_gap_prev));
        memcpy(b_gap_prev, b_gap_curr, sizeof(b_gap_prev));
    }
    return max_score;
}

static inline double calc_readmodel(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int *splice_pos, int *splice_offset, int n_splice, double baseline) {
    int i;
    int b; // array[width * row + col] = value
    int n = pos + read_length;
    double probability = 0;
    for (b = pos;  b < n; ++b) {
        if (b < 0) continue;

        int c = b;
        for (i = 0; i < n_splice; ++i) {
            if (b - pos > splice_pos[i]) c += splice_offset[i];
        }
        if (c >= seq_length) break;

        i = seq[c] - 'A';
        if (i < 0 || i >= 26) { exit_err("Character %c at pos %d (%d) not in valid alphabet\n", seq[c], c, seq_length); }

        probability += matrix[NT_CODES * (b - pos) + seqnt_map[i]]; 
        if (probability < baseline - 10) break; // stop if less than 1% contribution to baseline (best/highest) probability mass
    }
    return probability;
}

static inline double calc_prob(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int *splice_pos, int *splice_offset, int n_splice) {
    /* Get the sequence g in G and its neighborhood (half a read length flanking regions) */
    int n1 = pos - (read_length / 2);
    int n2 = pos + (read_length / 2);
    if (n1 < 0) n1 = 0;

    int i;
    double probability = 0;
    if (dp) {
        n2 += read_length;
        for (i = 0; i < n_splice; ++i) n2 += splice_offset[i];
        if (n2 >= seq_length) n2 = seq_length - 1;
        probability = smith_waterman_gotoh(matrix, read_length, seq, n2 - n1 + 1, n1);
    }
    else {
        probability = calc_readmodel(matrix, read_length, seq, seq_length, pos, splice_pos, splice_offset, n_splice, -1e6);
        double baseline = probability;
        for (i = n1; i < n2; ++i) {
            if (i != pos) {
                probability = log_add_exp(probability, calc_readmodel(matrix, read_length, seq, seq_length, i, splice_pos, splice_offset, n_splice, baseline));
                if (probability > baseline) baseline = probability;
            }
        }
    }
    return probability;
}

static void calc_likelihood(Variant **var_data, char *refseq, int refseq_length, Read **read_data, size_t nreads,
                            Vector_Int *c, double *ref, double *alt, double *het, int *ref_count, int *alt_count, size_t seti) {
    size_t i, readi;
    *ref = 0;
    *alt = 0;
    *het = 0;
    *ref_count = 0;
    *alt_count = 0;

    /* Alternative sequence */
    int altseq_length = 0;
    char *altseq = construct_altseq(refseq, refseq_length, c, var_data, &altseq_length);

    /* Aligned reads */
    for (readi = 0; readi < nreads; ++readi) {
        /* Only consider reads that map to all in current combination
        Variant *first = (Variant *)c->data[0];
        Variant *last = (Variant *)c->data[c->size - 1];
        if (read_data[readi]->pos > first->pos || read_data[readi]->pos + read_data[readi]->length < last->pos) continue;
        */
        int is_unmap = 0;
        int is_dup = 0;
        int is_reverse = 0;
        int is_secondary = 0;
        int n;
        char *s, token[strlen(read_data[readi]->flag) + 1];
        for (s = read_data[readi]->flag; sscanf(s, "%[^,]%n", token, &n) == 1; s += n + 1) {
            if (strcmp("UNMAP", token) == 0) is_unmap = 1;
            else if (strcmp("DUP", token) == 0) is_dup = 1;
            else if (strcmp("REVERSE", token) == 0) is_reverse = 1;
            else if (strcmp("SECONDARY", token) == 0 || strcmp("SUPPLEMENTARY", token) == 0) is_secondary = 1;
            if (*(s + n) != ',') break;
        }
        if (is_unmap) continue;
        if (nodup && is_dup) continue;
        if (pao && is_secondary) continue;

         /* Read probability matrix */
        double is_match[read_data[readi]->length], no_match[read_data[readi]->length];
        for (i = 0; i < read_data[readi]->length; ++i) {
            is_match[i] = p_match[read_data[readi]->qual[i]] + 1;
            no_match[i] = p_mismatch[read_data[readi]->qual[i]] + 1;
        }
        double readprobmatrix[read_data[readi]->length * NT_CODES];
        set_prob_matrix(readprobmatrix, read_data[readi]->qseq, read_data[readi]->length, is_match, no_match);
       
        /* 
        Exact Formuation:
        Probability that read is from an outside the reference paralogous "elsewhere", f in F.  Approximate the bulk of probability distribution P(r|f):
           a) perfect match = prod[ (1-e) ]
           b) hamming/edit distance 1 = prod[ (1-e) ] * sum[ (e/3) / (1-e) ]
        Length distribution, for reads with different lengths (hard clipped), where longer reads should have a relatively lower P(r|f):
           c) lengthfactor = alpha ^ (read length - expected read length)
        P(r|f) = (perfect + hamming_1) / lengthfactor 
        */
        double delta[read_data[readi]->length];
        for (i = 0; i < read_data[readi]->length; ++i) delta[i] = no_match[i] - is_match[i];
        double a = sum(is_match, read_data[readi]->length);
        double elsewhere = log_add_exp(a, a + log_sum_exp(delta, read_data[readi]->length)) - (LGALPHA * (read_data[readi]->length - read_data[readi]->inferred_length));
        /* Constant outside paralog term, testing, seems to perform near identically to the exact formulation above for fixed read length sequences */
        //double elsewhere = -2.4; // < 10%

        double pout = elsewhere;
        double prgu = calc_prob(readprobmatrix, read_data[readi]->length, refseq, refseq_length, read_data[readi]->pos, read_data[readi]->splice_pos, read_data[readi]->splice_offset, read_data[readi]->n_splice);
        double prgv = calc_prob(readprobmatrix, read_data[readi]->length, altseq, altseq_length, read_data[readi]->pos, read_data[readi]->splice_pos, read_data[readi]->splice_offset, read_data[readi]->n_splice);

        /* Multi-map alignments from XA tags: chr8,+42860367,97M3S,3;chr9,-44165038,100M,4; */
        if (read_data[readi]->multimapXA != NULL) {
            int xa_pos, n;
            char xa_chr[strlen(read_data[readi]->multimapXA) + 1];
            for (s = read_data[readi]->multimapXA; sscanf(s, "%[^,],%d,%*[^;]%n", xa_chr, &xa_pos, &n) == 2; s += n + 1) {
                Fasta *f = refseq_fetch(xa_chr, fa_file);
                if (f == NULL) continue;
                char *xa_refseq = f->seq;
                int xa_refseq_length = f->seq_length;

                double *p_readprobmatrix = readprobmatrix;
                double *newreadprobmatrix = NULL;
                if ((xa_pos < 0 && !is_reverse) || (xa_pos > 0 && is_reverse)) { // opposite of primary alignment strand
                    newreadprobmatrix = reverse(readprobmatrix, read_data[readi]->length * NT_CODES);
                    p_readprobmatrix = newreadprobmatrix;
                }

                xa_pos = abs(xa_pos);
                pout = log_add_exp(pout, elsewhere); // the more multi-mapped, the more likely it is the read is from elsewhere (paralogous), hence it scales (multiplied) with the number of multi-mapped locations
                double readprobability = calc_prob(p_readprobmatrix, read_data[readi]->length, xa_refseq, xa_refseq_length, xa_pos, read_data[readi]->splice_pos, read_data[readi]->splice_offset, read_data[readi]->n_splice);
                prgu = log_add_exp(prgu, readprobability);
                if (strcmp(xa_chr, read_data[readi]->chr) == 0) { // secondary alignments are in same chromosome as primary
                    readprobability = calc_prob(p_readprobmatrix, read_data[readi]->length, altseq, altseq_length, xa_pos, read_data[readi]->splice_pos, read_data[readi]->splice_offset, read_data[readi]->n_splice);
                }
                prgv = log_add_exp(prgv, readprobability);
                free(newreadprobmatrix); newreadprobmatrix = NULL;
                if (*(s + n) != ';') break;
            }
        }
        else if (read_data[readi]->multimapNH > 1) { // scale by the number of multimap positions
            double n = log(read_data[readi]->multimapNH - 1);
            double readprobability = prgu + n;
            pout = log_add_exp(pout, elsewhere + n);
            prgu = log_add_exp(prgu, readprobability);
            prgv = log_add_exp(prgv, readprobability);
        }

        /* Mixture model: probability that the read is from elsewhere, outside paralogous source */
        pout += lgomega;
        prgu = log_add_exp(pout, prgu);
        prgv = log_add_exp(pout, prgv);

        /* Mixture model: heterozygosity or heterogeneity as explicit allele frequency mu such that P(r|GuGv) = (mu)(P(r|Gv)) + (1-mu)(P(r|Gu)) */
        double phet   = log_add_exp(LG50 + prgv, LG50 + prgu);
        double phet10 = log_add_exp(LG10 + prgv, LG90 + prgu);
        double phet90 = log_add_exp(LG90 + prgv, LG10 + prgu);
        if (phet10 > phet) phet = phet10;
        if (phet90 > phet) phet = phet90;

        /* Read count incremented only when the difference in probability is not ambiguous, > ~log(2) difference */
        if (prgv > prgu && prgv - prgu > 0.69) *alt_count += 1;
        else if (prgu > prgv && prgu - prgv > 0.69) *ref_count += 1;

        if (verbose && prgv > read_data[readi]->prgv) {
            read_data[readi]->prgu = prgu;
            read_data[readi]->prgv = prgv;
            read_data[readi]->pout = pout;
            read_data[readi]->index = (int)seti;
        }

        /* Priors */
        *ref += prgu + ref_prior;
        *alt += prgv + alt_prior;
        *het += phet + het_prior;

        if (debug >= 2) {
            fprintf(stderr, ":%d:\t%f\t%f\t%f\t%f\t%d\t%d\t", (int)seti, prgu, phet, prgv, pout, *ref_count, *alt_count);
            fprintf(stderr, "%s\t%s\t%d\t", read_data[readi]->name, read_data[readi]->chr, read_data[readi]->pos);
            for (i = 0; i < read_data[readi]->n_cigar; ++i) fprintf(stderr, "%d%c ", read_data[readi]->cigar_oplen[i], read_data[readi]->cigar_opchr[i]);
            fprintf(stderr, "\t");
            for (i = 0; i < c->size; ++i) { Variant *curr = var_data[c->data[i]]; fprintf(stderr, "%s,%d,%s,%s;", curr->chr, curr->pos, curr->ref, curr->alt); }
            fprintf(stderr, "\t");
            if (read_data[readi]->multimapXA != NULL) fprintf(stderr, "%s\t", read_data[readi]->multimapXA);
            else fprintf(stderr, "%d\t", read_data[readi]->multimapNH);
            if (read_data[readi]->flag != NULL) fprintf(stderr, "%s\t", read_data[readi]->flag);
            fprintf(stderr, "%s\t", read_data[readi]->qseq);
            for (i = 0; i < read_data[readi]->length; ++i) fprintf(stderr, "%d ", read_data[readi]->qual[i]);
            fprintf(stderr, "\n");
        }
    }
    free(altseq); altseq = NULL;
    if (debug >= 1) {
        fprintf(stderr, "=%d=\t%f\t%f\t%f\t%d\t%d\t%d\t", (int)seti, *ref, *het, *alt, *ref_count, *alt_count, (int)nreads);
        for (i = 0; i < c->size; ++i) { Variant *curr = var_data[c->data[i]]; fprintf(stderr, "%s,%d,%s,%s;", curr->chr, curr->pos, curr->ref, curr->alt); } fprintf(stderr, "\n");
    }
}

static char *evaluate(const Vector *var_set) {
    size_t i, readi, seti;

    Variant **var_data = (Variant **)var_set->data;

    /* Reference sequence */
    Fasta *f = refseq_fetch(var_data[0]->chr, fa_file);
    if (f == NULL) return NULL;
    char *refseq = f->seq;
    int refseq_length = f->seq_length;

    /* Reads in variant region coordinates */
    Vector *read_list = bam_fetch(bam_file, var_data[0]->chr, var_data[0]->pos, var_data[var_set->size - 1]->pos);
    if (read_list->size == 0) {
        free(read_list); read_list = NULL;
        return NULL;
    }
    Read **read_data = (Read **)read_list->data;
    size_t nreads = read_list->size;

    /* Variant combinations as a Vector of Vectors */
    Vector *combo = powerset(var_set->size, maxh);
    Vector_Int **c = (Vector_Int **)combo->data;
    //for (seti = 0; seti < combo->size; ++seti) { fprintf(stderr, "%d\t", (int)seti); for (i = 0; i < c[seti]->size; ++i) { fprintf(stderr, "%d;", c[seti]->data[i]); } fprintf(stderr, "\t"); for (i = 0; i < c[seti]->size; ++i) { Variant *curr = var_data[c[seti]->data[i]]; fprintf(stderr, "%s,%d,%s,%s;", curr->chr, curr->pos, curr->ref, curr->alt); } fprintf(stderr, "\n"); }

    double ref = 0;
    Vector_Double *alt_list = vector_double_create(32);
    Vector_Double *het_list = vector_double_create(32);
    double *alt = alt_list->data;
    double *het = het_list->data;

    Vector_Int *ref_count_list = vector_int_create(32);
    Vector_Int *alt_count_list = vector_int_create(32);
    int *ref_count = ref_count_list->data;
    int *alt_count = alt_count_list->data;

    for (seti = 0; seti < combo->size; ++seti) { // all, singles, doubles
        vector_double_add(alt_list, 0);
        vector_double_add(het_list, 0);
        vector_int_add(ref_count_list, 0);
        vector_int_add(alt_count_list, 0);

        calc_likelihood(var_data, refseq, refseq_length, read_data, nreads, c[seti], &ref, &alt[seti], &het[seti], &ref_count[seti], &alt_count[seti], seti);
    }
    if (var_set->size > 2 && (int)combo->size - var_set->size - 1 < maxh) { // triples and beyond, if doubles doesn't already exceed the limit
        size_t j;
        seti = var_set->size;
        while (++seti < combo->size) {
            //fprintf(stderr, "seti %d\t%f\t%f\t%f\t%f\n", (int)seti, ref, alt[seti], het[seti], log_add_exp(alt[seti], het[seti]) - ref); int xx; for (xx = 0; xx < combo->size; ++xx) { fprintf(stderr, "%d\t", (int)xx); for (i = 0; i < c[xx]->size; ++i) { fprintf(stderr, "%d;", c[xx]->data[i]); } fprintf(stderr, "\t"); for (i = 0; i < c[xx]->size; ++i) { Variant *curr = var_data[c[xx]->data[i]]; fprintf(stderr, "%s,%d,%s,%s;", curr->chr, curr->pos, curr->ref, curr->alt); } fprintf(stderr, "\n"); } fprintf(stderr, "\n");
            j = combo->size; // index of potential first add to vector
            if (log_add_exp(alt[seti], het[seti]) - ref < -10) continue;
            derive_combo(combo, c[seti], var_set->size);
            for (i = j; i < combo->size; ++i) {
                calc_likelihood(var_data, refseq, refseq_length, read_data, nreads, c[i], &ref, &alt[i], &het[i], &ref_count[i], &alt_count[i], i);
                if (log_add_exp(alt[i], het[i]) - ref < -10) {
                    vector_del(combo, i);
                    vector_double_del(alt_list, i);
                    vector_double_del(het_list, i);
                    vector_int_del(ref_count_list, i);
                    vector_int_del(alt_count_list, i);
                    --i;
                }
            }
            if ((int)combo->size - var_set->size - 1 > maxh) break;
        }
    }

    if (verbose) {
        for (readi = 0; readi < nreads; ++readi) {
            if (read_data[readi]->prgu == 0 && read_data[readi]->prgv == 0 && read_data[readi]->pout == 0) continue; // unprocessed read
            flockfile(stderr);
            fprintf(stderr, "%s\t%s\t%d\t", read_data[readi]->name, read_data[readi]->chr, read_data[readi]->pos);
            fprintf(stderr, "%f\t%f\t%f\t", read_data[readi]->prgu, read_data[readi]->prgv, read_data[readi]->pout);
            for (i = 0; i < read_data[readi]->n_cigar; ++i) fprintf(stderr, "%d%c", read_data[readi]->cigar_oplen[i], read_data[readi]->cigar_opchr[i]);
            fprintf(stderr, "\t");
            if (read_data[readi]->multimapXA != NULL) fprintf(stderr, "%s\t", read_data[readi]->multimapXA);
            else fprintf(stderr, "%d\t", read_data[readi]->multimapNH);
            if (read_data[readi]->flag != NULL) fprintf(stderr, "%s\t", read_data[readi]->flag);
            else fprintf(stderr, "NONE\t");
            fprintf(stderr, "[");
            for (i = 0; i < c[read_data[readi]->index]->size; ++i) { Variant *curr = var_data[c[read_data[readi]->index]->data[i]]; fprintf(stderr, "%s,%d,%s,%s;", curr->chr, curr->pos, curr->ref, curr->alt); }
            fprintf(stderr, "]\n");
            funlockfile(stderr);
        }
    }

    double total = log_add_exp(ref, log_add_exp(alt[0], het[0]));
    for (seti = 1; seti < combo->size; ++seti) { total = log_add_exp(total, log_add_exp(ref, log_add_exp(alt[seti], het[seti]))); }

    char *output = malloc(sizeof *output);
    output[0] = '\0';
    if (mvh) { /* Max likelihood variant hypothesis */
        int max_seti = 0;
        double has_alt = log_add_exp(alt[0], het[0]);
        for (seti = 1; seti < combo->size; ++seti) { 
            double p = log_add_exp(alt[seti], het[seti]);
            if (p > has_alt) {
                has_alt = p;
                max_seti = seti;
            }
        }
        Vector *v = vector_create(var_set->size, VARIANT_T);
        for (i = 0; i < c[max_seti]->size; ++i) vector_add(v, var_data[c[max_seti]->data[i]]);
        variant_print(&output, v, 0, (int)nreads, ref_count[max_seti], alt_count[max_seti], total, has_alt, ref);
        vector_free(v);
    }
    else { /* Marginal probabilities & likelihood ratios*/
        for (i = 0; i < var_set->size; ++i) {
            double has_alt = 0;
            double not_alt = 0;
            int acount = -1;
            int rcount = -1;
            for (seti = 0; seti < combo->size; ++seti) {
                not_alt = (not_alt == 0) ? ref : log_add_exp(not_alt, ref);
                if (variant_find(c[seti], i) != -1) { // if variant is in this combination
                    has_alt = (has_alt == 0) ? log_add_exp(alt[seti], het[seti]) : log_add_exp(has_alt, log_add_exp(alt[seti], het[seti]));
                    if (alt_count[seti] > acount) {
                        acount = alt_count[seti];
                        rcount = ref_count[seti];
                    }
                }
                else {
                    not_alt = log_add_exp(not_alt, log_add_exp(alt[seti], het[seti]));
                }
            }
            variant_print(&output, var_set, i, (int)nreads, rcount, acount, total, has_alt, not_alt);
        }
    }
    vector_double_destroy(alt_list); alt_list = NULL;
    vector_double_destroy(het_list); het_list = NULL;
    vector_int_destroy(ref_count_list); ref_count_list = NULL;
    vector_int_destroy(alt_count_list); alt_count_list = NULL;
    vector_destroy(read_list); free(read_list); read_list = NULL;
    for (i = 0; i < combo->size; ++i) vector_int_destroy(combo->data[i]);
    vector_free(combo);
    return output;
}

typedef struct {
    Vector *queue, *results;
    pthread_mutex_t q_lock;
    pthread_mutex_t r_lock;
} Work;

static void *pool(void *work) {
    Work *w = (Work *)work;
    Vector *queue = (Vector *)w->queue;
    Vector *results = (Vector *)w->results;

    while (1) { //pthread_t ptid = pthread_self(); uint64_t threadid = 0; memcpy(&threadid, &ptid, min(sizeof(threadid), sizeof(ptid)));
        pthread_mutex_lock(&w->q_lock);
        Vector *var_set = (Vector *)vector_pop(queue);
        pthread_mutex_unlock(&w->q_lock);
        if (var_set == NULL) break;
        
        char *outstr = evaluate(var_set);
        if (outstr == NULL) continue;

        pthread_mutex_lock(&w->r_lock);
        vector_add(results, outstr);
        vector_free(var_set);
        pthread_mutex_unlock(&w->r_lock);
    }
    return NULL;
}

static void process(const Vector *var_list, FILE *out_fh) {
    size_t i, j;

    Variant **var_data = (Variant **)var_list->data;

    i = 0;
    Vector *var_set = vector_create(var_list->size, VOID_T);
    if (sharedr == 1) { /* Variants that share a read: shared with a given first variant */
        while (i < var_list->size) {
            Vector *curr = vector_create(8, VARIANT_T);
            vector_add(curr, var_data[i]);

            /* Reads in variant i region coordinates */
            int i_last = bam_fetch_last(bam_file, var_data[i]->chr, var_data[i]->pos, var_data[i]->pos);

            j = i + 1;
            while (j < var_list->size && strcmp(var_data[i]->chr, var_data[j]->chr) == 0) { // while last read in i will reach j
                if (var_data[j]->pos > i_last) break;
                vector_add(curr, var_data[j++]);
            }
            i = j;
            vector_add(var_set, curr);
        }
    }
    else if (sharedr == 2) { /* Variants that share a read: shared with any neighboring variant */
        while (i < var_list->size) {
            Vector *curr = vector_create(8, VARIANT_T);
            vector_add(curr, var_data[i]);

            j = i + 1;
            while (j < var_list->size && strcmp(var_data[i]->chr, var_data[j]->chr) == 0) { // while last read in i will reach j
                /* Reads in variant i region coordinates */
                int i_last = bam_fetch_last(bam_file, var_data[i]->chr, var_data[i]->pos, var_data[i]->pos);
                if (var_data[j]->pos > i_last) break;
                vector_add(curr, var_data[j]);
                ++i;
                ++j;
            }
            i = j;
            vector_add(var_set, curr);
        }
    }
    else { /* Variants that are close together as sets */
        while (i < var_list->size) {
            Vector *curr = vector_create(8, VARIANT_T);
            vector_add(curr, var_data[i]);
            size_t j = i + 1;
            while (distlim > 0 && j < var_list->size && strcmp(var_data[j]->chr, var_data[j - 1]->chr) == 0 && abs(var_data[j]->pos - var_data[j - 1]->pos) <= distlim) {
                if (maxdist > 0 && abs(var_data[j]->pos - var_data[i]->pos) > maxdist) break;
                vector_add(curr, var_data[j++]);
            }
            i = j;
            vector_add(var_set, curr);
        }
    }
    /* Heterozygous non-reference variants as separate entries */
    int flag_add = 1;
    while (flag_add) {
        flag_add = 0;
        size_t nsets = var_set->size;
        for (i = 0; i < nsets; ++i) {
            Vector *curr_set = (Vector *)var_set->data[i];
            if (curr_set->size == 1) continue;

            int flag_nonset = 1;
            for (j = 0; j < curr_set->size - 1; ++j) { // check if all entries have the same position
                Variant *curr = (Variant *)curr_set->data[j];
                Variant *next = (Variant *)curr_set->data[j + 1];
                if (curr->pos == next->pos && strcmp(curr->chr, next->chr) == 0 && strcmp(curr->ref, next->ref) == 0 && strcmp(curr->alt, next->alt) == 0) vector_del(curr_set, j + 1); // delete duplicate entries
                else if (curr->pos != next->pos) flag_nonset = 0;
            }
            if (flag_nonset) { // only 1 entry, with multiple heterozygous non-reference variants
                while (curr_set->size > 1) {
                    Variant *curr = (Variant *)vector_pop(curr_set);
                    Vector *dup = vector_create(8, VARIANT_T);
                    vector_add(dup, curr);
                    vector_add(var_set, dup);
                }
            }
            else { // multiple entries comprising a set
                for (j = 0; j < curr_set->size - 1; ++j) {
                    Variant *curr = (Variant *)curr_set->data[j];
                    Variant *next = (Variant *)curr_set->data[j + 1];
                    if (curr->pos == next->pos) {
                        flag_add = 1;
                        Vector *dup = vector_dup(curr_set);
                        vector_del(curr_set, j);
                        vector_del(dup, j + 1);
                        vector_add(var_set, dup);
                    }
                }
            }
        }
    } 
    if (sharedr == 1) { print_status("# Variants with shared reads to first in set: %i entries\t%s", (int)var_set->size, asctime(time_info)); }
    else if (sharedr == 2) { print_status("# Variants with shared reads to any in set: %i entries\t%s", (int)var_set->size, asctime(time_info)); }
    else { print_status("# Variants within %d (max window: %d) bp: %i entries\t%s", distlim, maxdist, (int)var_set->size, asctime(time_info)); }

    print_status("# Options: maxh=%d mvh=%d pao=%d isc=%d omega=%g\n", maxh, mvh, pao, isc, omega);
    print_status("# Start: %d threads \t%s\t%s", nthread, bam_file, asctime(time_info));

    Vector *queue = vector_create(var_set->size, VOID_T);
    Vector *results = vector_create(var_set->size, VOID_T);
    for (i = 0; i < var_set->size; ++i) {
        vector_add(queue, var_set->data[i]);
    }
    Work *w = malloc(sizeof (Work));
    w->queue = queue;
    w->results = results;

    pthread_mutex_init(&w->q_lock, NULL);
    pthread_mutex_init(&w->r_lock, NULL);

    pthread_t tid[nthread];
    for (i = 0; i < nthread; ++i) pthread_create(&tid[i], NULL, pool, w);
    for (i = 0; i < nthread; ++i) pthread_join(tid[i], NULL);

    pthread_mutex_destroy(&w->q_lock);
    pthread_mutex_destroy(&w->r_lock);

    free(w); w = NULL;
    vector_free(var_set);

    qsort(results->data, results->size, sizeof (void *), nat_sort_vector);
    fprintf(out_fh, "# SEQ\tPOS\tREF\tALT\tReads\tRefReads\tAltReads\tProb\tOdds\tSet\n");
    for (i = 0; i < results->size; ++i) fprintf(out_fh, "%s", (char *)results->data[i]);
    vector_destroy(queue); free(queue); queue = NULL;
    vector_destroy(results); free(results); results = NULL;
    print_status("# Done:\t%s\t%s", bam_file, asctime(time_info));
}

static void print_usage() {
    printf("\nUsage: eagle [options] -v variants.vcf -a alignment.bam -r reference.fasta\n\n");
    printf("Required:\n");
    printf("  -v --vcf      FILE   Variants VCF file. [stdin]\n");
    printf("  -a --bam      FILE   Alignment data bam files, ref-coord sorted with bai index file.\n");
    printf("  -r --ref      FILE   Reference sequence, fasta file with fai index file.\n");
    printf("Options:\n");
    printf("  -o --out      FILE   Output file. [stdout]\n");
    printf("  -t --nthread  INT    Number of threads. [1]\n");
    printf("  -s --sharedr  INT    Group nearby variants that share a read, 0:distance based/off, 1:shared with first, 2:shared with any. [0]\n");
    printf("  -n --distlim  INT    Group nearby variants within n bases, 0:off. [10]\n");
    printf("  -w --maxdist  INT    Maximum number of bases between any two variants in a set of hypotheses, 0:off. [0]\n");
    printf("  -m --maxh     INT    Maximum number of combinations in the set of hypotheses, instead of all 2^n. [1024]\n");
    printf("     --mvh             Output the maximum likelihood hypothesis in the set instead of marginal probabilities.\n");
    printf("     --pao             Primary alignments only.\n");
    printf("     --isc             Ignore soft-clipped bases.\n");
    printf("     --nodup           Ignore marked duplicate reads (based on SAM flag).\n");
    printf("     --splice          Allow spliced reads.\n");
    printf("     --dp              Use dynamic programming to calculate likelihood instead of the basic model.\n");
    //printf("     --match    INT    DP match score. [1]. Recommend 2 for long reads with indel errors.\n");
    //printf("     --mismatch INT    DP mismatch penalty. [4]. Recommend 5 for long reads with indel errors.\n");
    printf("     --gap_op   INT    DP gap open penalty. [6]. Recommend 2 for long reads with indel errors.\n");
    printf("     --gap_ex   INT    DP gap extend penalty. [1].\n");
    printf("     --verbose         Verbose mode, output likelihoods for each read seen for each hypothesis to stderr.\n");
    printf("     --hetbias  FLOAT  Prior probability bias towards non-homozygous mutations, between [0,1]. [0.5]\n");
    printf("     --omega    FLOAT  Prior probability of originating from outside paralogous source, between [0,1]. [1e-5]\n");
}

int main(int argc, char **argv) {
    /* Command line parameters defaults */
    vcf_file = NULL;
    bam_file = NULL;
    fa_file = NULL;
    out_file = NULL;
    nthread = 1;
    sharedr = 0;
    distlim = 10;
    maxdist = 0;
    maxh = 1024;
    mvh = 0;
    pao = 0;
    isc = 0;
    nodup = 0;
    dp = 0;
    verbose = 0;
    debug = 0;
    hetbias = 0.5;
    omega = 1.0e-5;

    match = 1;
    mismatch = 4;
    gap_op = 6;
    gap_ex = 1;

    static struct option long_options[] = {
        {"vcf", required_argument, NULL, 'v'},
        {"bam", required_argument, NULL, 'a'},
        {"ref", required_argument, NULL, 'r'},
        {"out", optional_argument, NULL, 'o'},
        {"nthread", optional_argument, NULL, 't'},
        {"sharedr", optional_argument, NULL, 's'},
        {"distlim", optional_argument, NULL, 'n'},
        {"maxdist", optional_argument, NULL, 'w'},
        {"maxh", optional_argument, NULL, 'm'},
        {"mvh", no_argument, &mvh, 1},
        {"pao", no_argument, &pao, 1},
        {"isc", no_argument, &isc, 1},
        {"nodup", no_argument, &nodup, 1},
        {"splice", no_argument, &splice, 1},
        {"dp", no_argument, &dp, 1},
        {"verbose", no_argument, &verbose, 1},
        {"debug", optional_argument, NULL, 'd'},
        {"match", optional_argument, NULL, 980},
        {"mismatch", optional_argument, NULL, 981},
        {"gap_op", optional_argument, NULL, 982},
        {"gap_ex", optional_argument, NULL, 983},
        {"hetbias", optional_argument, NULL, 990},
        {"omega", optional_argument, NULL, 991},
        {0, 0, 0, 0}
    };

    int opt = 0;
    while ((opt = getopt_long(argc, argv, "v:a:r:o:t:s:n:w:m:d:", long_options, &opt)) != -1) {
        switch (opt) {
            case 0: 
                //if (long_options[option_index].flag != 0) break;
                break;
            case 'v': vcf_file = optarg; break;
            case 'a': bam_file = optarg; break;
            case 'r': fa_file= optarg; break;
            case 'o': out_file = optarg; break;
            case 't': nthread = parse_int(optarg); break;
            case 's': sharedr = parse_int(optarg); break;
            case 'n': distlim = parse_int(optarg); break;
            case 'w': maxdist = parse_int(optarg); break;
            case 'm': maxh = parse_int(optarg); break;
            case 'd': debug = parse_int(optarg); break;
            case 980: match = parse_int(optarg); break;
            case 981: mismatch = parse_int(optarg); break;
            case 982: gap_op = parse_int(optarg); break;
            case 983: gap_ex = parse_int(optarg); break;
            case 990: hetbias = parse_float(optarg); break;
            case 991: omega = parse_float(optarg); break;
            default: exit_usage("Bad options");
        }
    }
    if (optind > argc) { exit_usage("Bad program call"); }

    FILE *vcf_fh = stdin;
    if (vcf_file != NULL) { // default vcf file handle is stdin unless a vcf file option is used
        vcf_fh = fopen(vcf_file, "r");
        if (vcf_fh == NULL) { exit_err("failed to open VCF file %s\n", vcf_file); }
    }
    else {
        vcf_file = "stdin";
    }
    if (bam_file == NULL) { exit_usage("Missing alignments given as BAM file!"); } 
    if (fa_file == NULL) { exit_usage("Missing reference genome given as Fasta file!"); }
    if (nthread < 1) nthread = 1;
    if (sharedr < 0 || sharedr > 2) sharedr = 0;
    if (distlim < 0) distlim = 10;
    if (maxdist < 0) maxdist = 0;
    if (maxh < 0) maxh = 0;
    if (match <= 0) match = 1;
    if (match <= 0) mismatch = 4;
    if (gap_op <= 0) gap_op = 6;
    if (gap_ex <= 0) gap_ex = 1;
    if (hetbias < 0 || hetbias > 1) hetbias = 0.5;
    if (omega < 0 || omega > 1) omega = 1e-5;
    lgomega = (log(omega) - log(1.0-omega));

    ref_prior = log(0.5);
    alt_prior = log(0.5 * (1 - hetbias));
    het_prior = log(0.5 * hetbias);

    FILE *out_fh = stdout;
    if (out_file != NULL) out_fh = fopen(out_file, "w"); // default output file handle is stdout unless output file option is used

    init_seqnt_map(seqnt_map);
    //if (dp) init_dp_q2p_table(p_match, p_mismatch, 50, match, mismatch);
    init_q2p_table(p_match, p_mismatch, 50);

    /* Start processing data */
    clock_t tic = clock();
    Vector *var_list = vcf_read(vcf_fh);
    print_status("# Read VCF: %s\t%i entries\t%s", vcf_file, (int)var_list->size, asctime(time_info));

    refseq_hash = kh_init(rsh);
    //fasta_read(fa_file);

    pthread_mutex_init(&refseq_lock, NULL);
    process(var_list, out_fh);
    if (out_file != NULL) fclose(out_fh);
    else fflush(stdout);
    pthread_mutex_destroy(&refseq_lock);

    khiter_t k;
    for (k = kh_begin(refseq_hash); k != kh_end(refseq_hash); ++k) {
        if (kh_exist(refseq_hash, k)) vector_destroy(&kh_val(refseq_hash, k));
    }
    kh_destroy(rsh, refseq_hash);
    vector_destroy(var_list); free(var_list); var_list = NULL;

    clock_t toc = clock();
    print_status("# CPU time (hr):\t%f\n", (double)(toc - tic) / CLOCKS_PER_SEC / 3600);

    return 0;
}
