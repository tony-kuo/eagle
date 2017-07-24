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
#include "util.h"

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

double log_sum_exp(const double *a, int size) {
    int i;
    double max_exp = a[0]; 
    for (i = 1; i < size; ++i) { 
        if (a[i] > max_exp) max_exp = a[i]; 
    }
    double s = 0;
    for (i = 0; i < size; ++i) s += exp(a[i] - max_exp);
    return log(s) + max_exp;
}
