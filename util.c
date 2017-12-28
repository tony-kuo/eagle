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
#include "vector.h"
#include "util.h"

#define M_1_LOG10E (1.0/M_LOG10E)
#define LG3 (log(3.0))

char *strdup1(const char *src) {
    size_t n = strlen(src) + 1;
    char *des = malloc(n * sizeof *des);
    des[n] = '\0';
    return des ? memcpy(des, src, n) : NULL;
}

void str_resize(char **str, size_t size) {
    char *p = realloc(*str, size * sizeof *str);
    if (p == NULL) { exit_err("failed to realloc in str_resize\n"); }
    else { *str = p; }
}

int has_numbers(const char *str) {
    while (*str != '\0') {
        if (isdigit(*str++) == 1) return 1;
    }
    return 0;
}

int parse_int(const char *str) {
    errno = 0;
    char *end;
    int num = strtol(str, &end, 0);
    if (end != str && *end != '\0') { exit_err("failed to convert '%s' to int with leftover string '%s'\n", str, end); }
    return num;
}

float parse_float(const char *str) {
    errno = 0;
    char *end;
    double num = strtof(str, &end);
    if (end != str && *end != '\0') { exit_err("failed to convert '%s' to float with leftover string '%s'\n", str, end); }
    return num;
}

double sum(const double *a, int size) {
    double s = 0;
    while (--size >= 0) s += a[size];
    return s;
}

double *reverse(double *a, int size) {
    int i = 0;
    double *b = malloc(size * sizeof *b);
    while (--size >= 0) b[i++] = a[size];
    return b;
}

double log_add_exp(double a, double b) {
    double max_exp = a > b ? a : b;
    return log(exp(a - max_exp) + exp(b - max_exp)) + max_exp;
}

double log_sum_exp(const double *a, size_t size) {
    int i;
    double max_exp = a[0]; 
    for (i = 1; i < size; ++i) { 
        if (a[i] > max_exp) max_exp = a[i]; 
    }
    double s = 0;
    for (i = 0; i < size; ++i) s += exp(a[i] - max_exp);
    return log(s) + max_exp;
}

void init_seqnt_map(int *seqnt_map) {
    /* Mapping table, symmetrical according to complement */
    memset(seqnt_map, 0, sizeof(int) * 26);

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

void init_q2p_table(double *p_match, double *p_mismatch, size_t size) {
    /* FastQ quality score to ln probability lookup table */
    int i;
    double a;
    for (i = 0; i < size; ++i) { 
        if (i == 0) a = -0.01;
        else a = (double)i / -10 * M_1_LOG10E; //convert to ln
        p_match[i] = log(1 - exp(a)); // log(1-err)
        p_mismatch[i] = a - LG3; // log(err/3)
     }
}

void combinations(vector_t *combo, int k, int n) {
    int i, c[k];
    for (i = 0; i < k; ++i) c[i] = i; // first combination
    while (1) { // while (next_comb(c, k, n)) {
        // record the combination
        vector_int_t *v = vector_int_create(k);
        for (i = 0; i < k; ++i) vector_int_add(v, c[i]);
        vector_add(combo, v);

        i = k - 1;
        c[i]++;
        while ((i >= 0 && i < k) && (c[i] >= n - k + 1 + i)) {
            i--;
            c[i]++;
        }
        /* Combination (n-k, n-k+1, ..., n) reached. No more combinations can be generated */
        if (c[0] > n - k) break; // return 0;
        /* c now looks like (..., x, n, n, n, ..., n), turn it into (..., x, x + 1, x + 2, ...) */
        for (i = i + 1; i < k; ++i) c[i] = c[i - 1] + 1;
        // return 1;
    }
}

void derive_combo(vector_t *combo, vector_int_t *prev, int n) { // Derive the combinations in k+1 that contain the previous elements
    if (prev->size + 1 >= n) return;

    int k = prev->size + 1;
    int i, c[k];

    for (i = 0; i < prev->size; ++i) c[i] = prev->data[i]; // first combination
    c[prev->size] = c[prev->size - 1] + 1;

    while (c[prev->size] < n) { // generate and record combinations
        //for (i = 0; i < k; ++i) { fprintf(stderr, "%d;", c[i]); } fprintf(stderr, "\n");
        vector_int_t *v = vector_int_create(k);
        for (i = 0; i < k; ++i) vector_int_add(v, c[i]);
        vector_add(combo, v);
        c[prev->size]++;
    }
    //int ii, jj; for (ii = 0; ii < combo->size; ++ii) { vector_int_t **c = (vector_int_t **)combo->data; fprintf(stderr, "%d\t", (int)ii); for (jj = 0; jj < c[ii]->size; ++jj) { fprintf(stderr, "%d;", c[ii]->data[jj]); } fprintf(stderr, "\n"); } fprintf(stderr, "\n");
}

vector_t *powerset(int n, int maxh) {
    vector_t *combo = vector_create(n + 1, VOID_T);
    if (n == 1) {
        combinations(combo, 1, n);
    }
    else if (n > 1) {
        combinations(combo, n, n);
        combinations(combo, 1, n);
        /*
        int k; 
        for (k = 2; k <= n - 1 && (int)combo->size - n - 1 < maxh; ++k) combinations(combo, k, n);
        */
        /*
        int i = 0;
        while (++i < combo->size) {
            vector_int_t **c = (vector_int_t **)combo->data; 
            //fprintf(stderr, "%d\t%d\t", i, n); int jj; for (jj = 0; jj < c[i]->size; ++jj) { fprintf(stderr, "%d;", c[i]->data[jj]); } fprintf(stderr, "\n"); 
            derive_combo(combo, c[i], n);
        }
        */
    }
    return combo;
}

int is_subset (int *arr1, int *arr2, int m, int n) { // Check if arr2 is a subset of arr1.  Requires sorted arrays
    int i = 0;
    int j = 0;

    if (m < n) return 0;
    while (i < n && j < m) {
        if (arr1[j] < arr2[i]) j++;
        else if (arr1[j] == arr2[i]) {
            j++;
            i++;
        }
        else if (arr1[j] > arr2[i]) return 0;
    }
    return (i < n) ? 0 : 1;
} 
