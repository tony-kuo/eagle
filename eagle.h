/*
EAGLE: explicit alternative genome likelihood evaluator
Given the sequencing data and candidate variant, explicitly test 
the alternative hypothesis against the reference hypothesis

Copyright 2016 Tony Kuo
This program is distributed under the terms of the GNU General Public License
*/

#ifndef _eagle_h_
#define _eagle_h_

#include <stdio.h>
#include <errno.h>
#include <string.h>

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

enum type {VOID_T, VARIANT_T, READ_T, FASTA_T};

typedef struct {
    int pos;
    char *chr, *ref, *alt;
} Variant;

typedef struct {
    int *cigar_oplen;
    int tid, pos, length, inferred_length, n_cigar, multimapNH;
    char *qseq, *chr, *name, *flag, *cigar_opchr, *multimapXA;
    double *qual;
    double prgu, prgv, pout;
} Read;

typedef struct {
    int seq_length;
    char *name, *seq;
} Fasta;

void freeVariant(Variant *v);
void freeRead(Read *r);
void freeFasta(Fasta *f);

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

Vector *vcf_read(FILE *file);
void fasta_read(const char *fa_file);
void process(const Vector *var_vector, char *bam_file, char *fa_file, FILE *out_fh);
#endif
