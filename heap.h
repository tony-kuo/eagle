/*
EAGLE: explicit alternative genome likelihood evaluator
Given the sequencing data and candidate variant, explicitly test 
the alternative hypothesis against the reference hypothesis

Copyright 2016 Tony Kuo
This program is distributed under the terms of the GNU General Public License
*/

#ifndef _heap_h_
#define _heap_h_

#include "vector.h"

typedef struct {
    double priority;
    void *data;
} node_t;

typedef struct {
    node_t *node;
    size_t len, size;
    enum type type;
} heap_t;

void heap_init(heap_t *a, enum type var_type);
heap_t *heap_create(enum type var_type);
void heap_free(heap_t *a);
void heap_destroy(heap_t *a);
void heap_push(heap_t *a, double priority, void *data);
void *heap_pop(heap_t *a);

#endif
