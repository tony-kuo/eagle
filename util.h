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

#define FNV_32_PRIME ((u_int32_t)0x01000193)
#define FNV1_32_INIT ((u_int32_t)0x811c9dc5)
#define MASK_16 (((u_int32_t)1<<16)-1) /* i.e., (u_int32_t)0xffff */
#define MASK_24 (((u_int32_t)1<<24)-1) /* i.e., (u_int32_t)0xffff */

char *strdup1(const char *src);
void str_resize(char **str, int size);

int has_numbers(const char *str);
int parse_int(const char *str);
float parse_float(const char *str);
double parse_double(const char *str);

int sum_i(const int *a, int size);
double sum_d(const double *a, int size);
double *reverse(double *a, int size);

double log_add_exp(double a, double b);
double log_sum_exp(const double *a, int size);

void combinations(vector_t *combo, int k, int n);
void derive_combo(vector_t *combo, vector_int_t *prev, int n);
vector_t *powerset(int n, int maxh);
vector_t *all_and_singletons(int n);
int is_subset (int arr1[], int arr2[], int m, int n);

u_int32_t fnv_32a_str(char *str);

#endif
