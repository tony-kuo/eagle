/*
EAGLE: explicit alternative genome likelihood evaluator
Given the sequencing data and candidate variant, explicitly test 
the alternative hypothesis against the reference hypothesis

Copyright 2016 Tony Kuo
This program is distributed under the terms of the GNU General Public License
*/

#ifndef _util_h_
#define _util_h_

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "vector.h"

#ifdef NDEBUG
#define debug(M, ...)
#else
#define debug(M, ...) fprintf(stderr, "DEBUG %s:%d: " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)
#endif
#define clean_errno() (errno == 0 ? "None" : strerror(errno))
#define log_err(M, ...) fprintf(stderr, "ERROR: (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__)
#define log_warn(M, ...) fprintf(stderr, "WARN: (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__)
#define exit_err(M, ...) fprintf(stderr, "ERROR: (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__); exit(EXIT_FAILURE)
#define exit_usage(M, ...) print_usage(); fprintf(stderr, "\n" M "\n"); exit(EXIT_FAILURE)

char *strdup1(const char *src);
void str_resize(char **str, size_t size);

int has_numbers(const char *str);
int parse_int(const char *str);
float parse_float(const char *str);

double sum(const double *a, int size);
double *reverse(double *a, int size);

double log_add_exp(double a, double b);
double log_sum_exp(const double *a, size_t size);

void init_seqnt_map(int *seqnt_map);
void init_q2p_table(double *p_match, double *p_mismatch, size_t size);
void init_dp_q2p_table(double *p_match, double *p_mismatch, size_t size, int match, int mismatch);

void combinations(Vector *combo, int k, int n);
void derive_combo(Vector *combo, Vector_Int *prev, int n);
Vector *powerset(int n, int maxh);
int is_subset (int arr1[], int arr2[], int m, int n);

#endif
