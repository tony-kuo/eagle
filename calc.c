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

#define M_1_LOG10E (1.0/M_LOG10E)
#define LG3 (log(3.0))

void init_q2p_table(double *p_match, double *p_mismatch, size_t size) {
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
    memset(seqnt_map, 0, sizeof (int) * 26);

    seqnt_map['A'-'A'] = 0;
    seqnt_map['C'-'A'] = 1;

    /* Ambiguous codes */
    seqnt_map['H'-'A'] = 2; // A, C, T
    seqnt_map['B'-'A'] = 3; // C, G, T
    seqnt_map['R'-'A'] = 4; // A, G
    seqnt_map['K'-'A'] = 5; // G, T
    seqnt_map['S'-'A'] = 6; // G, C
    seqnt_map['W'-'A'] = 7; // A, T

    seqnt_map['N'-'A'] = 8;
    seqnt_map['X'-'A'] = 8;

    // W also in 9, S also in 10
    seqnt_map['M'-'A'] = 11; // A, C
    seqnt_map['Y'-'A'] = 12; // C, T
    seqnt_map['V'-'A'] = 13; // A, C, G
    seqnt_map['D'-'A'] = 14; // A, G, T

    seqnt_map['G'-'A'] = 15;
    seqnt_map['T'-'A'] = 16;
    seqnt_map['U'-'A'] = 16;
}

void set_prob_matrix(double *matrix, const char *seq, int read_length, const double *is_match, const double *no_match, int *seqnt_map) {
    int i, b; // array[row * width + col] = value
    for (b = 0; b < read_length; b++) {
        for (i = 0; i < NT_CODES; i++) matrix[read_length * i + b] = no_match[b];
        matrix[read_length * seqnt_map[seq[b] - 'A'] + b] = is_match[b];
        switch (seq[b]) {
        case 'A':
            matrix[read_length * seqnt_map['M' - 'A'] + b] = is_match[b];
            matrix[read_length * seqnt_map['R' - 'A'] + b] = is_match[b];
            matrix[read_length * seqnt_map['V' - 'A'] + b] = is_match[b];
            matrix[read_length * seqnt_map['H' - 'A'] + b] = is_match[b];
            matrix[read_length * seqnt_map['D' - 'A'] + b] = is_match[b];
            matrix[read_length * seqnt_map['W' - 'A'] + b] = is_match[b];
            matrix[read_length * 9 + b] = is_match[b]; // also W
            break;
        case 'T':
            matrix[read_length * seqnt_map['K' - 'A'] + b] = is_match[b];
            matrix[read_length * seqnt_map['Y' - 'A'] + b] = is_match[b];
            matrix[read_length * seqnt_map['B' - 'A'] + b] = is_match[b];
            matrix[read_length * seqnt_map['H' - 'A'] + b] = is_match[b];
            matrix[read_length * seqnt_map['D' - 'A'] + b] = is_match[b];
            matrix[read_length * seqnt_map['W' - 'A'] + b] = is_match[b];
            matrix[read_length * 9 + b] = is_match[b]; // also W
            break;
        case 'C':
            matrix[read_length * seqnt_map['M' - 'A'] + b] = is_match[b];
            matrix[read_length * seqnt_map['Y' - 'A'] + b] = is_match[b];
            matrix[read_length * seqnt_map['B' - 'A'] + b] = is_match[b];
            matrix[read_length * seqnt_map['V' - 'A'] + b] = is_match[b];
            matrix[read_length * seqnt_map['H' - 'A'] + b] = is_match[b];
            matrix[read_length * seqnt_map['S' - 'A'] + b] = is_match[b];
            matrix[read_length * 10 + b] = is_match[b]; // also S
            break;
        case 'G':
            matrix[read_length * seqnt_map['K' - 'A'] + b] = is_match[b];
            matrix[read_length * seqnt_map['R' - 'A'] + b] = is_match[b];
            matrix[read_length * seqnt_map['B' - 'A'] + b] = is_match[b];
            matrix[read_length * seqnt_map['V' - 'A'] + b] = is_match[b];
            matrix[read_length * seqnt_map['D' - 'A'] + b] = is_match[b];
            matrix[read_length * seqnt_map['S' - 'A'] + b] = is_match[b];
            matrix[read_length * 10 + b] = is_match[b]; // also S
            break;
        }
    }
}

double calc_read_prob(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int *seqnt_map) {
    int i; // array[width * row + col] = value
    int end = (pos + read_length < seq_length) ? pos + read_length : seq_length;

    double probability[end - pos];
    for (i = pos;  i < end; i++) {
        int c = seq[i] - 'A';
        if (c < 0 || c >= 26) { exit_err("Character %c at pos %d (%d) not in valid alphabet\n", seq[i], i, seq_length); }

        probability[i - pos] = matrix[read_length * seqnt_map[c] + (i - pos)];
    }
    return sum_d(probability, end - pos);
}

double calc_prob_region(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int start, int end, int *seqnt_map) {
    if (start < 0) start = 0;
    if (end >= seq_length) end = seq_length - 1;

    int i;
    double p[end - start];
    for (i = start; i < end; i++) {
        p[i - start] = calc_read_prob(matrix, read_length, seq, seq_length, i, seqnt_map);
    }
    return sum_d(p, end - start);
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
            probability += calc_prob_region(submatrix, r_len, seq, seq_length, pos, start, end, seqnt_map);
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

    prev[0] = 0;
    prev[1] = 0 - gap_op;
    for (j = 2; j < read_length + 1; j++) prev[j] = prev[j - 1] - gap_ex;
    for (j = 0; j < read_length + 1; j++) b_gap_prev[j] = -DBL_MAX;

    double max_score = -DBL_MAX;
    for (i = start; i < end; i++) {
        double row_max = -DBL_MAX;
        double upleft, open, extend;

        curr[0] = 0;
        a_gap_curr[0] = -DBL_MAX;
        b_gap_curr[0] = -DBL_MAX;
        for (j = 1; j <= read_length; j++) {
            int c = seq[i] - 'A';
            if (c < 0 || c >= 26) { exit_err("Character %c at pos %d (%d) not in valid alphabet\n", seq[i], i, seq_length); }

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
            if (curr[j] > row_max) row_max = curr[j];
        }
        if (row_max > max_score) max_score = row_max;

        memcpy(prev, curr, sizeof (prev));
        memcpy(b_gap_prev, b_gap_curr, sizeof (b_gap_prev));
    }
    return max_score;
}

double calc_prob_region_dp(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int start, int end, int gap_op, int gap_ex, int *seqnt_map) {
    if (start < 0) start = 0;
    end += read_length;
    if (end >= seq_length) end = seq_length - 1;
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
            probability += calc_prob_region_dp(submatrix, r_len, seq, seq_length, pos, start, end, gap_op, gap_ex, seqnt_map);
            free(submatrix); submatrix = NULL;

            g_pos += r_len + splice_offset[i];
            r_pos = splice_pos[i] + 1;
        }
    }
    return probability;
}
