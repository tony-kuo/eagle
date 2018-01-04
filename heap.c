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
    a->nodes = malloc(4 * sizeof (node_t));
}

heap_t *heap_create(enum type var_type) {
    heap_t *a = malloc(sizeof (heap_t));
    heap_init(a, var_type);
    return a;
}

void heap_free(heap_t *a) {
    a->len = a->size = 0;
    free(a->nodes); a->nodes = NULL;
    free(a); a = NULL;
}

void heap_push(heap_t *a, double priority, void *data) {
    if (a->len + 1 >= a->size) {
        a->size *= 2;
        node_t *p = realloc(a->nodes, a->size * sizeof (node_t));
        if (p == NULL) { exit_err("failed to realloc in heap_push\n"); }
        else { a->nodes = p; }
    }
    int i = a->len + 1;
    int j = i / 2;
    while (i > 1 && a->nodes[j].priority < priority) {
        a->nodes[i] = a->nodes[j];
        i = j;
        j = j / 2;
    }
    a->nodes[i].priority = priority;
    a->nodes[i].data = data;
    a->len++;
}

void *heap_pop(heap_t *a) {
    int i, j, k;
    if (a->len == 0) return NULL;

    void *data = a->nodes[1].data;

    a->nodes[1] = a->nodes[a->len];
    double priority = a->nodes[1].priority;

    a->len--;

    i = 1;
    while (1) {
        k = i;
        j = 2 * i;
        if (j <= a->len && a->nodes[j].priority < priority) k = j;
        if (j + 1 <= a->len && a->nodes[j + 1].priority < a->nodes[k].priority) k = j + 1;
        if (k == i) break;
        a->nodes[i] = a->nodes[k];
        i = k;
    }
    a->nodes[i] = a->nodes[a->len + 1];
    return data;
}
