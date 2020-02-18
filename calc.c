/*
EAGLE: explicit alternative genome likelihood evaluator
Given the sequencing data and candidate variant, explicitly test 
the alternative hypothesis against the reference hypothesis

Copyright 2016 Tony Kuo
This program is distributed under the terms of the GNU General Public License
*/

#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include "calc.h"
//#include "calc_gpu.h"

#define M_1_LOG10E (1.0/M_LOG10E)
#define LG3 (log(3.0))

void init_q2p_table(double *p_match, double *p_mismatch, int size) {
    /* FastQ quality score to ln probability lookup table */
    int i;
    for (i = 0; i < size; i++) { 
        double a = (i == 0) ? -0.01 : (double)i / -10 * M_1_LOG10E; //convert to ln
        p_match[i] = log(1 - exp(a)); // log(1-err)
        p_mismatch[i] = a - LG3; // log(err/3)
     }
}

void init_seqnt_map(int *seqnt_map) {
    /* Mapping table, symmetrical according to complement */
    memset(seqnt_map, 0, sizeof (int) * 58);

    seqnt_map['A'-'A'] = 0;
    seqnt_map['C'-'A'] = 1;

    /* Ambiguous codes */
    seqnt_map['H'-'A'] = 2; // A, C, T
    seqnt_map['B'-'A'] = 3; // C, G, T
    seqnt_map['R'-'A'] = 4; // A, G
    seqnt_map['K'-'A'] = 5; // G, T
    seqnt_map['S'-'A'] = 6; // G, C
    seqnt_map['W'-'A'] = 7; // A, T

    seqnt_map['c'-'A'] = 8; // methylated C
    seqnt_map['t'-'A'] = 9; // unmethylated C

    seqnt_map['N'-'A'] = 10;
    seqnt_map['X'-'A'] = 10;

    seqnt_map['g'-'A'] = 11; // methylated C, other strand, i.e. G
    seqnt_map['a'-'A'] = 12; // unmethylated C, other strand, i.e. A

    // W also in 13, S also in 14
    seqnt_map['M'-'A'] = 15; // A, C
    seqnt_map['Y'-'A'] = 16; // C, T
    seqnt_map['V'-'A'] = 17; // A, C, G
    seqnt_map['D'-'A'] = 18; // A, G, T

    seqnt_map['G'-'A'] = 19;
    seqnt_map['T'-'A'] = 20;
    seqnt_map['U'-'A'] = 20;
}

void set_prob_matrix(double *matrix, const read_t *read, const double *is_match, const double *no_match, const int *seqnt_map, const int bisulfite) {
    int i, b; // array[row * width + col] = value
    for (b = 0; b < read->length; b++) {
        for (i = 0; i < NT_CODES; i++) matrix[read->length * i + b] = no_match[b];
        matrix[read->length * seqnt_map[read->qseq[b] - 'A'] + b] = is_match[b];
        switch (read->qseq[b]) {
        case 'A':
            matrix[read->length * seqnt_map['M' - 'A'] + b] = is_match[b];
            matrix[read->length * seqnt_map['R' - 'A'] + b] = is_match[b];
            matrix[read->length * seqnt_map['V' - 'A'] + b] = is_match[b];
            matrix[read->length * seqnt_map['H' - 'A'] + b] = is_match[b];
            matrix[read->length * seqnt_map['D' - 'A'] + b] = is_match[b];
            matrix[read->length * seqnt_map['W' - 'A'] + b] = is_match[b];
            matrix[read->length * 13 + b] = is_match[b]; // also W
            break;
        case 'T':
            matrix[read->length * seqnt_map['K' - 'A'] + b] = is_match[b];
            matrix[read->length * seqnt_map['Y' - 'A'] + b] = is_match[b];
            matrix[read->length * seqnt_map['B' - 'A'] + b] = is_match[b];
            matrix[read->length * seqnt_map['H' - 'A'] + b] = is_match[b];
            matrix[read->length * seqnt_map['D' - 'A'] + b] = is_match[b];
            matrix[read->length * seqnt_map['W' - 'A'] + b] = is_match[b];
            matrix[read->length * 13 + b] = is_match[b]; // also W
            break;
        case 'C':
            matrix[read->length * seqnt_map['M' - 'A'] + b] = is_match[b];
            matrix[read->length * seqnt_map['Y' - 'A'] + b] = is_match[b];
            matrix[read->length * seqnt_map['B' - 'A'] + b] = is_match[b];
            matrix[read->length * seqnt_map['V' - 'A'] + b] = is_match[b];
            matrix[read->length * seqnt_map['H' - 'A'] + b] = is_match[b];
            matrix[read->length * seqnt_map['S' - 'A'] + b] = is_match[b];
            matrix[read->length * 14 + b] = is_match[b]; // also S
            break;
        case 'G':
            matrix[read->length * seqnt_map['K' - 'A'] + b] = is_match[b];
            matrix[read->length * seqnt_map['R' - 'A'] + b] = is_match[b];
            matrix[read->length * seqnt_map['B' - 'A'] + b] = is_match[b];
            matrix[read->length * seqnt_map['V' - 'A'] + b] = is_match[b];
            matrix[read->length * seqnt_map['D' - 'A'] + b] = is_match[b];
            matrix[read->length * seqnt_map['S' - 'A'] + b] = is_match[b];
            matrix[read->length * 14 + b] = is_match[b]; // also S
            break;
        }
        if (bisulfite > 0) {
            switch (read->qseq[b]) {
            case 'A':
                matrix[read->length * seqnt_map['a' - 'A'] + b] = is_match[b]; // unmethylated reverse strand
                break;
            case 'T':
                matrix[read->length * seqnt_map['t' - 'A'] + b] = is_match[b]; // unmethylated forward strand
                break;
            case 'C':
                matrix[read->length * seqnt_map['c' - 'A'] + b] = is_match[b]; // methylated forward strand
                break;
            case 'G':
                matrix[read->length * seqnt_map['g' - 'A'] + b] = is_match[b]; // methylated reverse strand
                break;
            }
            if ((bisulfite == 1) && (read->qseq[b] == 'T') && ((!read->is_read2 && !read->is_reverse) || (read->is_read2 && read->is_reverse))) matrix[read->length * seqnt_map['C' - 'A'] + b] = is_match[b]; // unmethylated forward strand, top strand
            else if ((bisulfite == 2) && (read->qseq[b] == 'A') && ((!read->is_read2 && read->is_reverse) || (read->is_read2 && !read->is_reverse))) matrix[read->length * seqnt_map['G' - 'A'] + b] = is_match[b]; // unmethylated reverse strand, bottom strand
            else if ((bisulfite >= 3) && (read->qseq[b] == 'T') && ((!read->is_read2 && !read->is_reverse) || (read->is_read2 && read->is_reverse))) matrix[read->length * seqnt_map['C' - 'A'] + b] = is_match[b]; // unmethylated forward strand, top strand
            else if ((bisulfite >= 3) && (read->qseq[b] == 'A') && ((!read->is_read2 && read->is_reverse) || (read->is_read2 && !read->is_reverse))) matrix[read->length * seqnt_map['G' - 'A'] + b] = is_match[b]; // unmethylated reverse strand, bottom strand
        }
    }
}

double calc_read_prob(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int *seqnt_map) {
    int i; // array[width * row + col] = value
    int end = (pos + read_length < seq_length) ? pos + read_length : seq_length;

    double probability[end - pos];
    for (i = pos;  i < end; i++) {
        int c = seq[i] - 'A';
        if (c < 0 || c > 57 || (c > 25 && c < 32)) { exit_err("Character %c at pos %d (%d) not in valid alphabet\n", seq[i], i, seq_length); }

        probability[i - pos] = matrix[read_length * seqnt_map[c] + (i - pos)];
    }
    return sum_d(probability, end - pos);
}

double calc_prob_region(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int start, int end, int *seqnt_map) {
    if (start < 0) start = 0;
    else if (start >= seq_length) start = seq_length - 1;
    if (end < 0) end = 0;
    else if (end >= seq_length) end = seq_length - 1;

    int i;
    double p[end - start];
    for (i = start; i < end; i++) {
        p[i - start] = calc_read_prob(matrix, read_length, seq, seq_length, i, seqnt_map);
    }
    return log_sum_exp(p, end - start);
}

double calc_prob(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int *splice_pos, int *splice_offset, int n_splice, int *seqnt_map) {
    /* Get the sequence g in G and its neighborhood (half a read length flanking regions) */
    int start = pos - (read_length / 2);
    int end = pos + (read_length / 2);

    int i, j;
    double probability = 0;
    if (n_splice == 0) {
        probability = calc_prob_region(matrix, read_length, seq, seq_length, pos, start, end, seqnt_map);
    }
    else { // calculate the probability for each splice section separately
        int r_pos = 0;
        int g_pos = pos;
        for (i = 0; i <= n_splice; i++) {
            int r_len = (i < n_splice) ? splice_pos[i] - r_pos + 1 : read_length - r_pos;
            int n = r_len / 2;
            start = g_pos - n;
            end = g_pos + n;

            double *submatrix = malloc(NT_CODES * r_len * sizeof (double));
            for (j = 0; j < NT_CODES; j++) memcpy(&submatrix[r_len * j], &matrix[read_length * j + r_pos], r_len * sizeof (double));
            probability += calc_prob_region(submatrix, r_len, seq, seq_length, g_pos, start, end, seqnt_map);
            free(submatrix); submatrix = NULL;

            g_pos += r_len + splice_offset[i];
            r_pos = splice_pos[i] + 1;
        }
    }
    return probability;
}

double smith_waterman_gotoh(const double *matrix, int read_length, const char *seq, int seq_length, int start, int end, int gap_op, int gap_ex, int *seqnt_map) { /* short in long version */
    int i, j;

    double prev[read_length + 1], curr[read_length + 1];
    double a_gap_curr[read_length + 1];
    double b_gap_prev[read_length + 1], b_gap_curr[read_length + 1];

    for (j = 0; j < read_length + 1; j++) prev[j] = 0;
    for (j = 0; j < read_length + 1; j++) b_gap_prev[j] = 0;

    double max_score = 0;
    double x_drop = 20 + (8 * log(read_length)) - 1; // minimum alignment score - 1
    for (i = start; i < end; i++) {
        double row_max = 0;
        double upleft, open, extend;

        curr[0] = 0;
        a_gap_curr[0] = 0;
        b_gap_curr[0] = 0;
        for (j = 1; j <= read_length; j++) {
            int c = seq[i] - 'A';
            if (c < 0 || c > 57 || (c > 25 && c < 32)) { exit_err("Character %c at pos %d (%d) not in valid alphabet\n", seq[i], i, seq_length); }

            upleft = prev[j - 1] + matrix[read_length * seqnt_map[c] + (j - 1)];

            open = curr[j - 1] - gap_op;
            extend = a_gap_curr[j - 1] - gap_ex;
            a_gap_curr[j] = (open >= extend) ? open : extend;

            open = prev[j] - gap_op;
            extend = b_gap_prev[j] - gap_ex;
            b_gap_curr[j] = (open >= extend) ? open : extend;

            curr[j] = upleft;
            if (a_gap_curr[j] >= curr[j]) curr[j] = a_gap_curr[j];
            if (b_gap_curr[j] >= curr[j]) curr[j] = b_gap_curr[j];
            if (curr[j] < row_max - x_drop) break;
            else if (curr[j] > row_max) row_max = curr[j];
        }
        if (row_max < max_score - x_drop) break;
        else if (row_max > max_score) max_score = row_max;

        memcpy(prev, curr, sizeof (prev));
        memcpy(b_gap_prev, b_gap_curr, sizeof (b_gap_prev));
    }
    return max_score;
}

double calc_prob_region_dp(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int start, int end, int gap_op, int gap_ex, int *seqnt_map) {
    if (start < 0) start = 0;
    else if (start >= seq_length) start = seq_length - 1;
    end += read_length;
    if (end < 0) end = 0;
    else if (end >= seq_length) end = seq_length - 1;
    return smith_waterman_gotoh(matrix, read_length, seq, seq_length, start, end, gap_op, gap_ex, seqnt_map);
}

double calc_prob_dp(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int *splice_pos, int *splice_offset, int n_splice, int gap_op, int gap_ex, int *seqnt_map) {
    /* Get the sequence g in G and its neighborhood (half a read length flanking regions) */
    int start = pos - (read_length / 2);
    int end = pos + (read_length / 2);

    int i, j;
    double probability = 0;
    if (n_splice == 0) {
        probability = calc_prob_region_dp(matrix, read_length, seq, seq_length, pos, start, end, gap_op, gap_ex, seqnt_map);
    }
    else { // calculate the probability for each splice section separately
        int r_pos = 0;
        int g_pos = pos;
        for (i = 0; i <= n_splice; i++) {
            int r_len = (i < n_splice) ? splice_pos[i] - r_pos + 1 : read_length - r_pos;
            int n = r_len / 2;
            start = g_pos - n;
            end = g_pos + n;

            double *submatrix = malloc(NT_CODES * r_len * sizeof (double));
            for (j = 0; j < NT_CODES; j++) memcpy(&submatrix[r_len * j], &matrix[read_length * j + r_pos], r_len * sizeof (double));
            probability += calc_prob_region_dp(submatrix, r_len, seq, seq_length, g_pos, start, end, gap_op, gap_ex, seqnt_map);
            free(submatrix); submatrix = NULL;

            g_pos += r_len + splice_offset[i];
            r_pos = splice_pos[i] + 1;
        }
    }
    return probability;
}

void calc_prob_snps_region(double *prgu, double *prgv, vector_int_t *combo, variant_t **var_data, double *matrix, int read_length, const char *seq, int seq_length, int pos, int start, int end, int *seqnt_map) {
    if (start < 0) start = 0;
    else if (start >= seq_length) start = seq_length;
    if (end < 0) end = 0;
    else if (end >= seq_length) end = seq_length;

    int i, j, m;
    double prgu_i[end - start], prgv_i[end - start];
    //ALIGN_t *a = ALIGN_create(0, 0, matrix, read_length, seq, seq_length, pos, start, end, seqnt_map);
    //calc_read_prob_cpu(a, prgu_i);
    //ALIGN_destroy(a);
    for (i = start; i < end; i++) {
        int n = i - start;
        prgu_i[n] = calc_read_prob(matrix, read_length, seq, seq_length, i, seqnt_map); // reference probability per position i
        prgv_i[n] = prgu_i[n]; // alternative probability per position i

        int offset = 0;
        for (j = 0; j < combo->len; j++) {
            variant_t *v = var_data[combo->data[j]];
            int v_pos = v->pos - 1;
            int ref_len = strlen(v->ref);
            int alt_len = strlen(v->alt);
            if (v->ref[0] == '-') ref_len = 0;
            else if (v->alt[0] == '-') alt_len = 0;

            int l = (ref_len == alt_len) ? ref_len : read_length + ref_len + alt_len; // if snp(s), consider each change; if indel, consider the frameshift as a series of snps in the rest of the read
            for (m = 0; m < l; m++) {
                int g_pos = v_pos + m;
                int r_pos = g_pos - i + offset;
                if (r_pos < 0) continue;
                if (r_pos >= read_length || g_pos >= seq_length) break;

                int x;
                if (m >= ref_len) x = seq[g_pos] - 'A';
                else x = v->ref[m] - 'A';

                int y;
                if (m >= alt_len) {
                    if (g_pos + ref_len - alt_len >= seq_length) break;
                    y = seq[g_pos + ref_len - alt_len] - 'A';
                }
                else {
                    y = v->alt[m] - 'A';
                }

                if (x < 0 || x > 57 || (x > 25 && x < 32)) { exit_err("Ref character %c at gpos %d (%d) not in valid alphabet\n", seq[g_pos], g_pos, seq_length); }
                if (y < 0 || y > 57 || (y > 25 && y < 32)) { exit_err("Alt character %c at rpos %d for %s;%d;%s;%s not in valid alphabet\n", v->alt[m], m, v->chr, v->pos, v->ref, v->alt); }

                prgv_i[n] = prgv_i[n] - matrix[read_length * seqnt_map[x] + r_pos] + matrix[read_length * seqnt_map[y] + r_pos]; // update alternative array
            }
            offset += alt_len - ref_len;
        }
    }
    *prgu += log_sum_exp(prgu_i, end - start);
    *prgv += log_sum_exp(prgv_i, end - start);
}

void calc_prob_snps(double *prgu, double *prgv, vector_int_t *combo, variant_t **var_data, double *matrix, int read_length, const char *seq, int seq_length, int pos, int *splice_pos, int *splice_offset, int n_splice, int *seqnt_map) {
    /* Get the sequence g in G and its neighborhood (half a read length flanking regions) */
    int start = pos - (read_length / 2);
    int end = pos + (read_length / 2);

    *prgu = 0;
    *prgv = 0;

    int i, j;
    if (n_splice == 0) {
        calc_prob_snps_region(prgu, prgv, combo, var_data, matrix, read_length, seq, seq_length, pos, start, end, seqnt_map);
    }
    else { // calculate the probability for each splice section separately
        int r_pos = 0;
        int g_pos = pos;
        for (i = 0; i <= n_splice; i++) {
            int r_len = (i < n_splice) ? splice_pos[i] - r_pos + 1 : read_length - r_pos;
            start = g_pos - (r_len / 2);
            end = g_pos + (r_len / 2);

            double *submatrix = malloc(NT_CODES * r_len * sizeof (double));
            for (j = 0; j < NT_CODES; j++) memcpy(&submatrix[r_len * j], &matrix[read_length * j + r_pos], r_len * sizeof (double));
            calc_prob_snps_region(prgu, prgv, combo, var_data, submatrix, r_len, seq, seq_length, pos, start, end, seqnt_map);
            free(submatrix); submatrix = NULL;

            g_pos += r_len + splice_offset[i];
            r_pos = splice_pos[i] + 1;
        }
    }
}
