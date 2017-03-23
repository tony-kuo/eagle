/*
EAGLE: explicit alternative genome likelihood evaluator
Given the sequencing data and candidate variant, explicitly test 
the alternative hypothesis against the reference hypothesis

Copyright 2016 Tony Kuo
This program is distributed under the terms of the GNU General Public License
*/

#include <stdlib.h>
#include <ctype.h>

#include "util.h"
#include "vector.h"

void vector_init(Vector *a, size_t initial_size, enum type var_type) {
    a->size = 0;
    a->type = var_type;
    a->capacity = initial_size;
    a->data = malloc(initial_size * sizeof (void *));
}

void vector_add(Vector *a, void *entry) {
    if (a->size >= a->capacity) {
        a->capacity *= 2;
        void **p = realloc(a->data, a->capacity * sizeof (void *));
        if (p == NULL) { exit_err("failed to realloc in vector_add\n"); }
        else { a->data = p; }
    }
    a->data[a->size++] = entry;
}

void vector_del(Vector *a, int i) {
    a->data[i] = NULL;
    if (i == --a->size) return;
    memmove(&(a->data[i]), &(a->data[i + 1]), (a->size - i) * sizeof (void *));
    a->data[a->size] = NULL;
}

void *vector_pop(Vector *a) {
    if (a->size <= 0) return NULL;
    void *entry = a->data[--a->size];
    a->data[a->size] = NULL;
    return entry;
}

Vector *vector_create(size_t initial_size, enum type var_type) {
    Vector *a = malloc(sizeof (Vector));
    vector_init(a, initial_size, var_type);
    return a;
}

Vector *vector_dup(Vector *a) {
    Vector *v = vector_create(a->capacity, a->type);
    v->size = a->size;
    memcpy(&(v->data[0]), &(a->data[0]), a->size * sizeof (void *));
    return v;
}

void vector_free(Vector *a) {
    a->size = 0;
    free(a->data); a->data = NULL;
}

void vector_destroy(Vector *a) {
    size_t i;
    enum type var_type = a->type;
    for (i = 0; i < a->size; ++i) {
        switch (var_type) {
            case VARIANT_T:
                variant_destroy((Variant *)a->data[i]);
                break;
            case READ_T:
                read_destroy((Read *)a->data[i]);
                break;
            case FASTA_T:
                fasta_destroy((Fasta *)a->data[i]);
                break;
            default:
                break;
        }
        free(a->data[i]); a->data[i] = NULL;
    }
    a->size = a->capacity = 0;
    free(a->data); a->data = NULL;
}

Variant *variant_create(char *chr, int pos, char *ref, char *alt) {
    Variant *v = malloc(sizeof (Variant));
    v->chr = strdup(chr);
    v->pos = pos;
    v->ref = strdup(ref);
    v->alt = strdup(alt);
    return v;
}

void variant_destroy(Variant *v) {
    if (v == NULL) return;
    v->pos = 0;
    free(v->chr); v->chr = NULL;      
    free(v->ref); v->ref = NULL;
    free(v->alt); v->alt = NULL;
}

Read *read_create(char *name, int tid, char *chr, int pos) {
    Read *r = malloc(sizeof (Read));
    r->name = strdup(name);
    r->tid = tid;
    r->chr = strdup(chr);
    r->pos = pos;
    r->prgu = -10000;
    r->prgv = -10000;
    r->pout = -10000;
    r->maxseti = 0;
    r->var_list = vector_create(8, VOID_T);

    r->length = r->n_cigar = r->inferred_length = r->multimapNH = 0;
    r->qseq = NULL;
    r->qual = NULL;
    r->flag = NULL;
    r->cigar_opchr = NULL;
    r->cigar_oplen = NULL;
    r->multimapXA = NULL;
    return r;
}

void read_destroy(Read *r) {
    if (r == NULL) return;
    r->tid = r->pos = r->length = r->n_cigar = r->inferred_length = r->multimapNH = 0;
    r->prgu = r->prgv = r->pout = 0;
    r->maxseti = 0;
    free(r->name); r->name = NULL;
    free(r->chr); r->chr = NULL;
    free(r->qseq); r->qseq = NULL;
    free(r->qual); r->qual = NULL;
    free(r->flag); r->flag = NULL;
    free(r->cigar_opchr); r->cigar_opchr = NULL;
    free(r->cigar_oplen); r->cigar_oplen = NULL;
    free(r->multimapXA); r->multimapXA = NULL;
    vector_destroy(r->var_list); free(r->var_list); r->var_list = NULL;
}

Fasta *fasta_create(char *name) {
    Fasta *f = malloc(sizeof (Fasta));
    f->name = strdup(name);
    return f;
}

void fasta_destroy(Fasta *f) {
    if (f == NULL) return;
    f->seq_length = 0;
    free(f->seq); f->seq = NULL;      
    free(f->name); f->name = NULL;
}

int nat_sort_cmp(const void *a, const void *b, enum type var_type) {
    char *str1, *str2;
    switch (var_type) {
        case VARIANT_T: {
            Variant *c1 = *(Variant **)a;
            Variant *c2 = *(Variant **)b;
            if (strcasecmp(c1->chr, c2->chr) == 0) return (c1->pos > c2->pos) - (c1->pos < c2->pos);
            str1 = strdup(c1->chr);
            str2 = strdup(c2->chr);
            c1 = NULL;
            c2 = NULL;
            break;
        }
        case STR_T: {
            str1 = strdup((char *)a);
            str2 = strdup((char *)b);
            break;
        }
        default:
            str1 = strdup(*(char **)a);
            str2 = strdup(*(char **)b);
            break;
    }
    char *s1 = str1;
    char *s2 = str2;
    int cmp = 0;
    while (cmp == 0 && *s1 != '\0' && *s2 != '\0') {
        if (isspace(*s1) && isspace(*s2)) { // ignore whitespace
            s1 += 1;
            s2 += 1;
        }
        else if ((isalpha(*s1) && isalpha(*s2)) || (ispunct(*s1) && ispunct(*s2))) { // compare alphabet and punctuation
            *s1 = tolower(*s1);
            *s2 = tolower(*s2);
            cmp = (*s1 > *s2) - (*s1 < *s2);
            s1 += 1;
            s2 += 1;
        }
        else { // compare digits
            int i1, i2, n1, n2, t1, t2;
            t1 = sscanf(s1, "%d%n", &i1, &n1);
            if (t1 == 0) t1 = sscanf(s1, "%*[^0123456789]%d%n", &i1, &n1);
            t2 = sscanf(s2, "%d%n", &i2, &n2);
            if (t2 == 0) t2 = sscanf(s2, "%*[^0123456789]%d%n", &i2, &n2);

            if (t1 < 1 || t2 < 1) { // one string has no digits
                cmp = strcmp(s1, s2);
            }
            else {
                cmp = (i1 > i2) - (i1 < i2);
                if (cmp == 0) { // first set of digits are equal, check further
                    s1 += n1;
                    s2 += n2;
                }
            }
        }
    }
    free(str1); str1 = NULL;
    free(str2); str2 = NULL;
    return cmp;
}

int nat_sort_vector(const void *a, const void *b) {
    return nat_sort_cmp(a, b, VOID_T);
}

int nat_sort_variant(const void *a, const void *b) {
    return nat_sort_cmp(a, b, VARIANT_T);
}
