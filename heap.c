/*
EAGLE: explicit alternative genome likelihood evaluator
Given the sequencing data and candidate variant, explicitly test 
the alternative hypothesis against the reference hypothesis

Copyright 2016 Tony Kuo
This program is distributed under the terms of the GNU General Public License
*/

#include <stdlib.h>
#include "heap.h"
#include "util.h"
#include "vector.h"

void heap_init(heap_t *a, enum type var_type) {
    a->len = 0;
    a->size = 4;
    a->type = var_type;
    a->node = malloc(4 * sizeof (node_t));
}

heap_t *heap_create(enum type var_type) {
    heap_t *a = malloc(sizeof (heap_t));
    heap_init(a, var_type);
    return a;
}

void heap_free(heap_t *a) {
    a->len = a->size = 0;
    free(a->node); a->node = NULL;
    free(a); a = NULL;
}

void heap_destroy(heap_t *a) {
    size_t i;
    enum type var_type = a->type;
    for (i = 1; i <= a->len; i++) {
        switch (var_type) {
            case STATS_T:
                stats_destroy((stats_t *)a->node[i].data);
                break;
            case VARIANT_T:
                variant_destroy((variant_t *)a->node[i].data);
                break;
            case READ_T:
                read_destroy((read_t *)a->node[i].data);
                break;
            case FASTA_T:
                fasta_destroy((fasta_t *)a->node[i].data);
                break;
            default:
                break;
        }
        free(a->node[i].data); a->node[i].data = NULL;
    }
    a->len = a->size = 0;
    free(a->node); a->node = NULL;
}

void heap_push(heap_t *a, double priority, void *data) {
    if (a->len + 1 >= a->size) {
        a->size *= 2;
        node_t *p = realloc(a->node, a->size * sizeof (node_t));
        if (p == NULL) { exit_err("failed to realloc in heap_push\n"); }
        else { a->node = p; }
    }
    size_t i = a->len + 1;
    size_t j = i / 2;
    while (i > 1 && a->node[j].priority < priority) {
        a->node[i] = a->node[j];
        i = j;
        j = j / 2;
    }
    a->node[i].priority = priority;
    a->node[i].data = data;
    a->len++;
}

void *heap_pop(heap_t *a) {
    size_t i, j, k;
    if (a->len == 0) return NULL;

    void *data = a->node[1].data;

    a->node[1] = a->node[a->len];
    double priority = a->node[1].priority;

    a->len--;

    i = 1;
    while (1) {
        k = i;
        j = 2 * i;
        if (j <= a->len && a->node[j].priority < priority) k = j;
        if (j + 1 <= a->len && a->node[j + 1].priority < a->node[k].priority) k = j + 1;
        if (k == i) break;
        a->node[i] = a->node[k];
        i = k;
    }
    a->node[i] = a->node[a->len + 1];
    return data;
}
