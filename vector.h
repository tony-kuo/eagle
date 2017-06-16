/*
EAGLE: explicit alternative genome likelihood evaluator
Given the sequencing data and candidate variant, explicitly test 
the alternative hypothesis against the reference hypothesis

Copyright 2016 Tony Kuo
This program is distributed under the terms of the GNU General Public License
*/

#ifndef _vector_h_
#define _vector_h_

enum type {VOID_T, STR_T, VARIANT_T, READ_T, FASTA_T};

typedef struct {
    size_t size, capacity;
    enum type type;
    void **data;
} Vector;

void vector_init(Vector *a, size_t initial_size, enum type var_type);
void vector_add(Vector *a, void *entry);
void vector_del(Vector *a, int i);
void *vector_pop(Vector *a);
Vector *vector_create(size_t initial_size, enum type var_type);
Vector *vector_dup(Vector *a);
void vector_free(Vector *a);
void vector_destroy(Vector *a);

typedef struct {
    int pos;
    char *chr, *ref, *alt;
} Variant;

Variant *variant_create(char *chr, int pos, char *ref, char *alt);
void variant_destroy(Variant *v);

typedef struct {
    int *cigar_oplen, *splice_pos, *splice_offset;
    int tid, pos, length, inferred_length, n_cigar, n_splice, multimapNH;
    char *qseq, *chr, *name, *flag, *cigar_opchr, *multimapXA;
    double *qual;
    double prgu, prgv, pout;
    int index;
    Vector *var_list;
} Read;

Read *read_create(char *name, int tid, char *chr, int pos);
void read_destroy(Read *r);

typedef struct {
    int seq_length;
    char *name, *seq;
} Fasta;

Fasta *fasta_create(char *name);
void fasta_destroy(Fasta *f);

int nat_sort_cmp(const void *a, const void *b, enum type var_type);
int nat_sort_vector(const void *a, const void *b);
int nat_sort_variant(const void *a, const void *b);

#endif
