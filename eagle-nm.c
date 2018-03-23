/*
EAGLE: explicit alternative genome likelihood evaluator
Utility program that calculates the likelihood that there is no mutation in the given genome regions

Given the sequencing data and genomic regions in the reference, 
explicitly test the reference hypothesis against total hamming 1

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
#define LG3 (log(3.0))
#define LG50 (log(0.5))
#define LG10 (log(0.1))
#define LG90 (log(0.9))
#define LGALPHA (log(ALPHA))

/* Command line arguments */
static char *bed_file;
static char *bam_file;
static char *fa_file;
static char *out_file;
static int nthread;
static int pao;
static int isc;
static int nodup;
static int splice;
static int verbose;
static int debug;
static double ref_prior, alt_prior, het_prior;
static double mut_prior, nomut_prior;

/* Time info */
static time_t now; 
static struct tm *time_info; 
#define print_status(M, ...) time(&now); time_info = localtime(&now); fprintf(stderr, M, ##__VA_ARGS__);

/* Mapping table */
static int seqnt_map[26];
static char NT[4] = "ATGC";

/* Fastq quality to probability table */
static double p_match[50], p_mismatch[50];

KHASH_MAP_INIT_STR(rsh, vector_t)   // hashmap: string key, vector value
static khash_t(rsh) *refseq_hash; // pointer to hashmap
static pthread_mutex_t refseq_lock; 

vector_t *bed_read(FILE *file) {
    vector_t *reg_list = vector_create(64, REGION_T);

    char *line = NULL;
    ssize_t read_file = 0;
    size_t line_length = 0;
    while ((read_file = getline(&line, &line_length, file)) != -1) {
        if (line_length <= 0 || line[strspn(line, " \t\v\r\n")] == '\0') continue; // blank line
        if (line[0] == '#') continue;

        int pos1, pos2;
        char chr[line_length];
        int t = sscanf(line, "%s %d %d", chr, &pos1, &pos2);
        if (t < 3) { exit_err("bad fields in BED file\n"); }

        region_t *g = region_create(chr, pos1, pos2);
        vector_add(reg_list, g);
    }
    free(line); line = NULL;
    fclose(file);
    qsort(reg_list->data, reg_list->len, sizeof (void *), nat_sort_region);
    return reg_list;
}

vector_t *bam_fetch(const char *bam_file, const char *chr, const int pos1, const int pos2) {
    /* Reads in region coordinates */
    vector_t *read_list = vector_create(64, READ_T);

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
            read_t *read = read_create((char *)aln->data, aln->core.tid, bam_header->target_name[aln->core.tid], aln->core.pos);

            char *flag = bam_flag2str(aln->core.flag);
            if (flag != NULL) read->flag = strdup(flag);
            else read->flag = NULL;
            free(flag); flag = NULL;

            int n;
            char *s, token[strlen(read->flag) + 1];
            for (s = read->flag; sscanf(s, "%[^,]%n", token, &n) == 1; s += n + 1) {
                if (strcmp("UNMAP", token) == 0) read->is_unmap = 1;
                else if (strcmp("DUP", token) == 0) read->is_dup = 1;
                else if (strcmp("REVERSE", token) == 0) read->is_reverse = 1;
                else if (strcmp("SECONDARY", token) == 0 || strcmp("SUPPLEMENTARY", token) == 0) read->is_secondary = 1;
                if (*(s + n) != ',') break;
            }
            if (read->is_unmap || (nodup && read->is_dup) || (pao && read->is_secondary)) {
                read_destroy(read);
                continue;
            }

            int start_align = 0;
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
            }
            read->cigar_opchr[read->n_cigar] = '\0';
            read->inferred_length = bam_cigar2qlen(read->n_cigar, cigar);
            read->n_splice = j;

            if (!isc) {
                read->pos -= s_offset; // compensate for soft clip in mapped position
                s_offset = 0;
                e_offset = 0;
            }
            read->length = aln->core.l_qseq - (s_offset + e_offset);
            read->qseq = malloc((read->length + 1) * sizeof read->qseq);
            read->qual = malloc(read->length  * sizeof read->qual);
            uint8_t *qual = bam_get_qual(aln);
            for (i = 0; i < read->length; i++) {
                read->qseq[i] = toupper(seq_nt16_str[bam_seqi(bam_get_seq(aln), i + s_offset)]); // get nucleotide id and convert into IUPAC id.
                read->qual[i] = (qual[i] > 41) ? qual[i] - 31 : qual[i]; // account for phred64
            }
            read->qseq[read->length] = '\0';

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

fasta_t *refseq_fetch(char *name, const char *fa_file) {
    pthread_mutex_lock(&refseq_lock);
    size_t i;
	khiter_t k = kh_get(rsh, refseq_hash, name);
    if (k != kh_end(refseq_hash)) {
        vector_t *node = &kh_val(refseq_hash, k);
        fasta_t **f = (fasta_t **)node->data;
        for (i = 0; i < node->len; i++) {
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

    fasta_t *f = fasta_create(name);
    f->seq = fai_fetch(fai, name, &f->seq_length);
    char *s;
    for (s = f->seq; *s != '\0'; s++) *s = toupper(*s);

    int absent;
    k = kh_put(rsh, refseq_hash, f->name, &absent);
    vector_t *node = &kh_val(refseq_hash, k);
    if (absent) vector_init(node, 8, FASTA_T);
    vector_add(node, f);
    fai_destroy(fai);
    pthread_mutex_unlock(&refseq_lock);
    return f;
}

void set_prob_matrix(double *matrix, const char *seq, int read_length, const double *is_match, const double *no_match) {
    int i, b; // array[width * row + col] = value
    for (b = 0; b < read_length; b++) {
        for (i = 0; i < NT_CODES; i++) matrix[NT_CODES * b + i] = no_match[b];
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

static inline double calc_read_prob(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int baseline) {
    int i, c; // array[width * row + col] = value
    double probability = 0;
    for (i = pos;  i < pos + read_length; i++) {
        if (i < 0) continue;
        if (i >= seq_length) break;

        c = seq[i] - 'A';
        if (c < 0 || c >= 26) { exit_err("Character %c at pos %d (%d) not in valid alphabet\n", seq[i], i, seq_length); }

        probability += matrix[NT_CODES * (i - pos) + seqnt_map[c]]; 
        if (probability < baseline - 10) break; // stop if less than 1% contribution to baseline (best/highest) probability mass
    }
    return probability;
}

static inline double calc_prob_region(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int start, int end) {
    double probability = 0;
    if (start < 0) start = 0;

    probability = calc_read_prob(matrix, read_length, seq, seq_length, pos, -1e6);
    double baseline = probability;

    int i;
    for (i = start; i < end; i++) {
        if (i >= seq_length) break;
        if (i != pos) {
            probability = log_add_exp(probability, calc_read_prob(matrix, read_length, seq, seq_length, i, baseline));
            if (probability > baseline) baseline = probability;
        }
    }
    return probability;
}

static inline double calc_prob(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int *splice_pos, int *splice_offset, int n_splice) {
    /* Get the sequence g in G and its neighborhood (half a read length flanking regions) */
    int n = read_length / 2;
    int start = pos - n;
    int end = pos + n;

    int i;
    double probability = 0;
    if (n_splice == 0) {
        probability = calc_prob_region(matrix, read_length, seq, seq_length, pos, start, end);
    }
    else { // calculate the probability for each splice section separately
        int r_pos = 0;
        int g_pos = pos;
        int r_len;
        double *submatrix;
        for (i = 0; i <= n_splice; i++) {
            if (i < n_splice) r_len = splice_pos[i] - r_pos + 1;
            else r_len = read_length -  1 - r_pos + 1;
            n = r_len / 2;
            start = g_pos - n;
            end = g_pos + n;

            submatrix = malloc(NT_CODES * r_len * sizeof(double));
            memcpy(submatrix, &matrix[NT_CODES * r_pos], NT_CODES * r_len * sizeof(double));
            probability += calc_prob_region(submatrix, r_len, seq, seq_length, pos, start, end);
            free(submatrix); submatrix = NULL;

            g_pos += r_len + splice_offset[i];
            r_pos = splice_pos[i] + 1;
        }
    }
    return probability;
}

static inline void calc_prob_snps_region(double *prgu, double *prgv, int g_pos, const double *matrix, int read_length, const char *seq, int seq_length, int pos, int start, int end) {
    int i, r_pos, k, x;

    double prgu_n = 0;
    double prgv_n = 0;
    double probability, p;

    if (start < 0) start = 0;

    vector_double_t *r = vector_double_create(abs(end - start)); // reference probability per position i
    vector_double_t *a = vector_double_create(abs(end - start)); // alternative probability per position i
    for (i = start; i < end; i++) {
        if (i >= seq_length || g_pos >= seq_length) break;
        probability = calc_read_prob(matrix, read_length, seq, seq_length, i, -1e6);
        vector_double_add(r, probability);
        vector_double_add(a, probability);

        r_pos = g_pos - pos;
        if (r_pos >= 0 && r_pos < read_length) {
            x = seq[g_pos] - 'A';
            if (x < 0 || x >= 26) { exit_err("Ref character %c at pos %d (%d) not in valid alphabet\n", seq[g_pos], g_pos, seq_length); }

            probability = 0;
            for (k = 0; k < 4; k++) {
                if (NT[k] != seq[g_pos]) {
                    p = matrix[NT_CODES * r_pos + seqnt_map[NT[k] - 'A']];
                    probability = (probability == 0) ? p : log_add_exp(probability, p);
                    //printf("%c %f\t", NT[k], p);
                }
            }
            a->data[i - start] = a->data[i - start] - matrix[NT_CODES * r_pos + seqnt_map[x]] + probability; // update alternative array
            //printf("%d\t%d\t%c\t%d\t%f\t%f\t%f\t%f\n", i, g_pos, seq[g_pos], r_pos, r->data[i - start], a->data[i - start], matrix[NT_CODES * r_pos + seqnt_map[x]], probability);
        }
        prgu_n = (prgu_n == 0) ? r->data[i - start] : log_add_exp(prgu_n, r->data[i - start]);
        prgv_n = (prgv_n == 0) ? a->data[i - start] : log_add_exp(prgv_n, a->data[i - start]);
    }
    *prgu += prgu_n;
    *prgv += prgv_n;
    vector_double_free(r);
    vector_double_free(a);
}

static inline void calc_prob_snps(double *prgu, double *prgv, int g_pos, const double *matrix, int read_length, const char *seq, int seq_length, int pos, int *splice_pos, int *splice_offset, int n_splice) {
    /* Get the sequence g in G and its neighborhood (half a read length flanking regions) */
    int start = pos; // - (read_length / 2);
    int end = pos + 1; //(read_length / 2);

    *prgu = 0;
    *prgv = 0;

    int i;
    if (n_splice == 0) {
        calc_prob_snps_region(prgu, prgv, g_pos, matrix, read_length, seq, seq_length, pos, start, end);
    }
    else { // calculate the probability for each splice section separately
        int r_pos = 0;
        int g_pos = pos;
        int r_len;
        double *submatrix;
        for (i = 0; i <= n_splice; i++) {
            if (i < n_splice) r_len = splice_pos[i] - r_pos + 1;
            else r_len = read_length -  1 - r_pos + 1;
            start = g_pos - (r_len / 2);
            end = g_pos + (r_len / 2);

            submatrix = malloc(NT_CODES * r_len * sizeof(double));
            memcpy(submatrix, &matrix[NT_CODES * r_pos], NT_CODES * r_len * sizeof(double));
            calc_prob_snps_region(prgu, prgv, g_pos, submatrix, r_len, seq, seq_length, pos, start, end);
            free(submatrix); submatrix = NULL;

            g_pos += r_len + splice_offset[i];
            r_pos = splice_pos[i] + 1;
        }
    }
}

static char *evaluate_nomutation(const region_t *g) {
    size_t i, readi;

    /* Reference sequence */
    fasta_t *f = refseq_fetch(g->chr, fa_file);
    if (f == NULL) return NULL;
    char *refseq = f->seq;
    int refseq_length = f->seq_length;

    /* Reads in region coordinates */
    vector_t *read_list = bam_fetch(bam_file, g->chr, g->pos1, g->pos2);
    if (read_list->len == 0) {
        free(read_list); read_list = NULL;
        return NULL;
    }
    read_t **read_data = (read_t **)read_list->data;

    int g_pos;
    double alt_probability = -1e6;
    double ref_probability = 0;
    int ref_count = 0;
    int alt_count = 0;
    for (g_pos = g->pos1; g_pos <= g->pos2; g_pos++) {
        double ref = 0;
        double alt = 0;
        int r_count = 0;
        int a_count = 0;

        /* Aligned reads */
        for (readi = 0; readi < read_list->len; readi++) {
            //if (read_data[readi]->pos < g->pos1 || read_data[readi]->pos > g->pos2) continue;

            double is_match[read_data[readi]->length], no_match[read_data[readi]->length];
            for (i = 0; i < read_data[readi]->length; i++) {
                is_match[i] = p_match[read_data[readi]->qual[i]];
                no_match[i] = p_mismatch[read_data[readi]->qual[i]];
            }
            double readprobmatrix[NT_CODES * read_data[readi]->length];
            set_prob_matrix(readprobmatrix, read_data[readi]->qseq, read_data[readi]->length, is_match, no_match);

            double prgu, prgv;
            //printf("\t%d\t%s\t%d\t%d\t%s\n", g_pos, read_data[readi]->name, read_data[readi]->pos, read_data[readi]->length, read_data[readi]->qseq);
            calc_prob_snps(&prgu, &prgv, g_pos, readprobmatrix, read_data[readi]->length, refseq, refseq_length, read_data[readi]->pos, read_data[readi]->splice_pos, read_data[readi]->splice_offset, read_data[readi]->n_splice);
            //printf("%f\t%f\n\n", prgu, prgv);

            /* Multi-map alignments from XA tags: chr8,+42860367,97M3S,3;chr9,-44165038,100M,4; */
            if (read_data[readi]->multimapXA != NULL) {
                int xa_pos, n;
                char *s, xa_chr[strlen(read_data[readi]->multimapXA) + 1];
                for (s = read_data[readi]->multimapXA; sscanf(s, "%[^,],%d,%*[^;]%n", xa_chr, &xa_pos, &n) == 2; s += n + 1) {
                    if (strcmp(xa_chr, read_data[readi]->chr) != 0 && abs(xa_pos - read_data[readi]->pos) < read_data[readi]->length) { // if secondary alignment does not overlap primary aligment
                        fasta_t *f = refseq_fetch(xa_chr, fa_file);
                        if (f == NULL) continue;
                        char *xa_refseq = f->seq;
                        int xa_refseq_length = f->seq_length;

                        double *p_readprobmatrix = readprobmatrix;
                        double *newreadprobmatrix = NULL;
                        if ((xa_pos < 0 && !read_data[readi]->is_reverse) || (xa_pos > 0 && read_data[readi]->is_reverse)) { // opposite of primary alignment strand
                            newreadprobmatrix = reverse(readprobmatrix, read_data[readi]->length * NT_CODES);
                            p_readprobmatrix = newreadprobmatrix;
                        }

                        xa_pos = abs(xa_pos);
                        double readprobability = calc_prob(p_readprobmatrix, read_data[readi]->length, xa_refseq, xa_refseq_length, xa_pos, read_data[readi]->splice_pos, read_data[readi]->splice_offset, read_data[readi]->n_splice);
                        prgu = log_add_exp(prgu, readprobability);
                        prgv = log_add_exp(prgv, readprobability);
                        free(newreadprobmatrix); newreadprobmatrix = NULL;
                    }
                    if (*(s + n) != ';') break;
                }
            }
            else if (read_data[readi]->multimapNH > 1) { // scale by the number of multimap positions
                double n = log(read_data[readi]->multimapNH - 1);
                double readprobability = prgu + n;
                prgu = log_add_exp(prgu, readprobability);
                prgv = log_add_exp(prgv, readprobability);
            }

            if (prgu == 0 && prgv == 0) continue;

            if (prgu > prgv && prgu - prgv > 0.69) r_count += 1;
            else if (prgv > prgu && prgv - prgu > 0.69) a_count += 1;

            double phet   = log_add_exp(LG50 + prgv, LG50 + prgu);
            double phet10 = log_add_exp(LG10 + prgv, LG90 + prgu);
            double phet90 = log_add_exp(LG90 + prgv, LG10 + prgu);
            if (phet10 > phet) phet = phet10;
            if (phet90 > phet) phet = phet90;

            /* Priors */
            ref += prgu + ref_prior;
            alt += log_add_exp(prgv + alt_prior, phet + het_prior);

            if (debug > 1) {
                fprintf(stderr, "::\t%d\t%s\t%d\t%f\t%f\t%f\t%d\t%d\t", g_pos, read_data[readi]->name, read_data[readi]->pos, prgu, prgv, prgu-prgv, r_count, a_count);
                for (i = 0; i < read_data[readi]->n_cigar; i++) fprintf(stderr, "%d%c ", read_data[readi]->cigar_oplen[i], read_data[readi]->cigar_opchr[i]);
                fprintf(stderr, "\n");
            }
        }
        ref_probability = ref;
        alt_probability = (alt_probability == 0) ? alt : log_add_exp(alt_probability, alt);
        if (a_count > alt_count) {
            ref_count = r_count;
            alt_count = a_count;
        }
        if (debug > 0) fprintf(stderr, "++\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\n", g->pos1, g->pos2, g_pos, (int)read_list->len, ref_count, alt_count, ref, alt, ref - alt, alt_probability);
    }

    ref_probability += nomut_prior;
    alt_probability += mut_prior;
    double odds = (ref_probability - alt_probability) * M_1_LN10;
    size_t n = snprintf(NULL, 0, "%s\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\n", g->chr, g->pos1, g->pos2, (int)read_list->len, ref_count, alt_count, ref_probability, alt_probability, odds) + 1;
    char *output = malloc(n * sizeof *output);
    snprintf(output, n, "%s\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\n", g->chr, g->pos1, g->pos2, (int)read_list->len, ref_count, alt_count, ref_probability, alt_probability, odds);

    vector_destroy(read_list); free(read_list); read_list = NULL;
    return output;
}

typedef struct {
    vector_t *queue, *results;
    pthread_mutex_t q_lock;
    pthread_mutex_t r_lock;
} work_t;

static void *pool(void *work) {
    work_t *w = (work_t *)work;
    vector_t *queue = (vector_t *)w->queue;
    vector_t *results = (vector_t *)w->results;

    while (1) { //pthread_t ptid = pthread_self(); uint64_t threadid = 0; memcpy(&threadid, &ptid, min(sizeof(threadid), sizeof(ptid)));
        pthread_mutex_lock(&w->q_lock);
        region_t *g = (region_t *)vector_pop(queue);
        pthread_mutex_unlock(&w->q_lock);
        if (g == NULL) break;
        
        char *outstr = evaluate_nomutation(g);
        if (outstr == NULL) continue;

        pthread_mutex_lock(&w->r_lock);
        vector_add(results, outstr);
        pthread_mutex_unlock(&w->r_lock);
    }
    return NULL;
}

static void process(const vector_t *reg_list, FILE *out_fh) {
    size_t i;

    region_t **reg_data = (region_t **)reg_list->data;
    size_t nregions = reg_list->len;

    print_status("# Options: pao=%d isc=%d nodup=%d splice=%d\n", pao, isc, nodup, splice);
    print_status("# Start: %d threads \t%s\t%s", nthread, bam_file, asctime(time_info));

    vector_t *queue = vector_create(nregions, VOID_T);
    vector_t *results = vector_create(nregions, VOID_T);
    for (i = 0; i < nregions; i++) {
        vector_add(queue, reg_data[i]);
    }
    work_t *w = malloc(sizeof (work_t));
    w->queue = queue;
    w->results = results;

    pthread_mutex_init(&w->q_lock, NULL);
    pthread_mutex_init(&w->r_lock, NULL);

    pthread_t tid[nthread];
    for (i = 0; i < nthread; i++) pthread_create(&tid[i], NULL, pool, w);
    for (i = 0; i < nthread; i++) pthread_join(tid[i], NULL);

    pthread_mutex_destroy(&w->q_lock);
    pthread_mutex_destroy(&w->r_lock);

    free(w); w = NULL;

    qsort(results->data, results->len, sizeof (void *), nat_sort_vector);
    fprintf(out_fh, "# CHR\tPOS1\tPOS2\tReads\tRefCount\tAltCount\tRefProb\tAltProb\tOdds\n");
    for (i = 0; i < results->len; i++) fprintf(out_fh, "%s", (char *)results->data[i]);
    vector_destroy(queue); free(queue); queue = NULL;
    vector_destroy(results); free(results); results = NULL;
    print_status("# Done:\t%s\t%s", bam_file, asctime(time_info));
}

static void print_usage() {
    printf("\nUsage: eagle [options] -v regions.bed -a alignment.bam -r reference.fasta\n\n");
    printf("Required:\n");
    printf("  -v --bed        FILE   Genome regions, BED file. [stdin]\n");
    printf("  -a --bam        FILE   Alignment data bam files, ref-coord sorted with bai index file.\n");
    printf("  -r --ref        FILE   Reference sequence, fasta file with fai index file.\n");
    printf("Options:\n");
    printf("  -o --out        FILE   Output file. [stdout]\n");
    printf("  -t --nthread    INT    Number of threads. [1]\n");
    printf("     --pao               Primary alignments only.\n");
    printf("     --isc               Ignore soft-clipped bases.\n");
    printf("     --nodup             Ignore marked duplicate reads (based on SAM flag).\n");
    printf("     --splice            Allow spliced reads.\n");
    printf("     --mut_prior  FLOAT  Prior probability for a mutation at any given reference position [0.001].\n");
    printf("     --verbose           Verbose mode, output likelihoods for each read seen for each hypothesis to stderr.\n");
}

int main(int argc, char **argv) {
    /* Command line parameters defaults */
    bed_file = NULL;
    bam_file = NULL;
    fa_file = NULL;
    out_file = NULL;
    nthread = 1;
    pao = 0;
    isc = 0;
    nodup = 0;
    verbose = 0;
    debug = 0;

    mut_prior = 0.001;

    static struct option long_options[] = {
        {"bed", required_argument, NULL, 'v'},
        {"bam", required_argument, NULL, 'a'},
        {"ref", required_argument, NULL, 'r'},
        {"out", optional_argument, NULL, 'o'},
        {"nthread", optional_argument, NULL, 't'},
        {"pao", no_argument, &pao, 1},
        {"isc", no_argument, &isc, 1},
        {"nodup", no_argument, &nodup, 1},
        {"splice", no_argument, &splice, 1},
        {"verbose", no_argument, &verbose, 1},
        {"debug", optional_argument, NULL, 'd'},
        {"mut_prior", optional_argument, NULL, 990},
        {0, 0, 0, 0}
    };

    int opt = 0;
    while ((opt = getopt_long(argc, argv, "v:a:r:o:t:d:", long_options, &opt)) != -1) {
        switch (opt) {
            case 0: 
                //if (long_options[option_index].flag != 0) break;
                break;
            case 'v': bed_file = optarg; break;
            case 'a': bam_file = optarg; break;
            case 'r': fa_file= optarg; break;
            case 'o': out_file = optarg; break;
            case 't': nthread = parse_int(optarg); break;
            case 'd': debug = parse_int(optarg); break;
            case 990: mut_prior = parse_float(optarg); break;
            default: exit_usage("Bad options");
        }
    }
    if (optind > argc) { exit_usage("Bad program call"); }

    FILE *bed_fh = stdin;
    if (bed_file != NULL) { // default bed file handle is stdin unless a bed file option is used
        bed_fh = fopen(bed_file, "r");
        if (bed_fh == NULL) { exit_err("failed to open BED file %s\n", bed_file); }
    }
    else {
        bed_file = "stdin";
    }
    if (bam_file == NULL) { exit_usage("Missing alignments given as BAM file!"); } 
    if (fa_file == NULL) { exit_usage("Missing reference genome given as Fasta file!"); }
    if (nthread < 1) nthread = 1;
    if (mut_prior < 0 || mut_prior > 1) mut_prior = 0.001;
    nomut_prior = 1 - mut_prior;

    ref_prior = log(0.5);
    alt_prior = log(0.25);
    het_prior = log(0.25);

    FILE *out_fh = stdout;
    if (out_file != NULL) out_fh = fopen(out_file, "w"); // default output file handle is stdout unless output file option is used

    init_seqnt_map(seqnt_map);
    mut_prior = log(mut_prior);
    nomut_prior = log(nomut_prior);
    init_q2p_table(p_match, p_mismatch, 50);
    
    /* Start processing data */
    clock_t tic = clock();
    vector_t *reg_list = bed_read(bed_fh);
    print_status("# Read BED: %s\t%i entries\t%s", bed_file, (int)reg_list->len, asctime(time_info));

    refseq_hash = kh_init(rsh);

    pthread_mutex_init(&refseq_lock, NULL);
    process(reg_list, out_fh);
    if (out_file != NULL) fclose(out_fh);
    else fflush(stdout);
    pthread_mutex_destroy(&refseq_lock);

    khiter_t k;
    for (k = kh_begin(refseq_hash); k != kh_end(refseq_hash); k++) {
        if (kh_exist(refseq_hash, k)) vector_destroy(&kh_val(refseq_hash, k));
    }
    kh_destroy(rsh, refseq_hash);
    vector_destroy(reg_list); free(reg_list); reg_list = NULL;

    clock_t toc = clock();
    print_status("# CPU time (hr):\t%f\n", (double)(toc - tic) / CLOCKS_PER_SEC / 3600);

    return 0;
}
