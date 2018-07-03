/*
EAGLE: explicit alternative genome likelihood evaluator
Given the sequencing data and candidate variant, explicitly test 
the alternative hypothesis against the reference hypothesis

Copyright 2016 Tony Kuo
This program is distributed under the terms of the GNU General Public License
*/

#ifndef _calc_h_
#define _calc_h_

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "util.h"

#define NT_CODES 17   // Size of nucleotide code table

void init_seqnt_map(int *seqnt_map);
void init_q2p_table(double *p_match, double *p_mismatch, size_t size);

void set_prob_matrix(double *matrix, const char *seq, int read_length, const double *is_match, const double *no_match, int *seqnt_map);
double calc_read_prob(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int *seqnt_map);
double calc_prob_region(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int start, int end, int *seqnt_map);
double calc_prob(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int *splice_pos, int *splice_offset, int n_splice, int *seqnt_mp);
double smith_waterman_gotoh(const double *matrix, int read_length, const char *seq, int seq_length, int start, int end, int gap_op, int gap_ex, int *seqnt_map);
double calc_prob_region_dp(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int start, int end, int gap_op, int gap_ex, int *seqnt_map);
double calc_prob_dp(const double *matrix, int read_length, const char *seq, int seq_length, int pos, int *splice_pos, int *splice_offset, int n_splice, int gap_op, int gap_ex, int *seqnt_map);
#endif
