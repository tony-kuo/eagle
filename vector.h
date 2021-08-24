/*
EAGLE: explicit alternative genome likelihood evaluator
Given the sequencing data and candidate variant, explicitly test 
the alternative hypothesis against the reference hypothesis

Copyright 2016 Tony Kuo
This program is distributed under the terms of the GNU General Public License
*/

#ifndef _vector_h_
#define _vector_h_

#include <stdint.h>
#include "htslib/sam.h"

enum type {VOID_T, STATS_T, VARIANT_T, READ_T, FASTA_T, REGION_T};

typedef struct {
    size_t len, size;
    void **data;
    enum type type;
} vector_t;

void vector_init(vector_t *a, size_t initial_size, enum type var_type);
vector_t *vector_create(size_t initial_size, enum type var_type);
void vector_free(vector_t *a);
void vector_destroy(vector_t *a);
void vector_add(vector_t *a, void *entry);
void vector_del(vector_t *a, size_t i);
void *vector_pop(vector_t *a);
vector_t *vector_dup_shallow(vector_t *a);

typedef struct {
    size_t len, size;
    int *data;
} vector_int_t;

void vector_int_init(vector_int_t *a, size_t initial_size);
vector_int_t *vector_int_create(size_t initial_size);
void vector_int_free(vector_int_t *a);
void vector_int_add(vector_int_t *a, int entry);
void vector_int_del(vector_int_t *a, size_t i);
vector_int_t *vector_int_dup(vector_int_t *a);

typedef struct {
    double *data;
    size_t len, size;
} vector_double_t;

void vector_double_init(vector_double_t *a, size_t initial_size);
vector_double_t *vector_double_create(size_t initial_size);
void vector_double_free(vector_double_t *a);
void vector_double_add(vector_double_t *a, double entry);
void vector_double_del(vector_double_t *a, size_t i);
vector_double_t *vector_double_dup(vector_double_t *a);

typedef struct {
    int pos;
    char *chr, *ref, *alt;
} variant_t;

variant_t *variant_create(char *chr, int pos, char *ref, char *alt);
void variant_destroy(variant_t *v);

typedef struct {
    vector_t *var_list;
    float prgu, prgv, pout;
    int32_t *qual, *cigar_oplen, *splice_pos, *splice_offset;
    int32_t tid, pos, end, length, inferred_length, n_cigar, n_splice, multimapNH;
    int16_t index;
    int8_t is_dup, is_reverse, is_secondary, is_read2;
    char *qseq, *chr, *name, *flag, *cigar_opchr, *multimapXA;
} read_t;

read_t *read_create(char *name, int tid, char *chr, int pos);
void read_destroy(read_t *r);

typedef struct {
    int seq_length;
    char *name, *seq;
} fasta_t;

fasta_t *fasta_create(char *name);
void fasta_destroy(fasta_t *f);

typedef struct {
    int pos1, pos2;
    char *chr;
} region_t;

region_t *region_create(char *chr, int pos1, int pos2);
void region_destroy(region_t *g);

typedef struct {
    vector_int_t *combo;
    vector_double_t *read_prgv;
    double ref, alt, het, mut;
    int ref_count, alt_count, seen;
} stats_t;

stats_t *stats_create(vector_int_t *combo, int nreads);
void stats_destroy(stats_t *s);

int nat_sort_cmp(const void *a, const void *b, enum type var_type);
int nat_sort_vector(const void *a, const void *b);
int nat_sort_variant(const void *a, const void *b);
int nat_sort_region(const void *a, const void *b);

read_t *read_fetch(bam_hdr_t *bam_header, bam1_t *aln, int pao, int isc, int nodup, int splice, int phred64, int const_qual);

#endif
