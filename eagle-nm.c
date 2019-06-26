/*
EAGLE: explicit alternative genome likelihood evaluator
Utility program that calculates the likelihood that there is no mutation in the given genome regions

Given the sequencing data and genomic regions in the reference, 
explicitly test the reference hypothesis against total hamming 1

Copyright 2016 Tony Kuo
This program is distributed under the terms of the GNU General Public License
*/

#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <pthread.h>

#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/khash.h"
#include "vector.h"
#include "util.h"
#include "calc.h"

/* Constants */
#define ALPHA 1.3     // Factor to account for longer read lengths lowering the probability a sequence matching an outside paralogous source

/* Precalculated log values */
#define M_1_LOG10E (1.0/M_LOG10E)
#define M_1_LN10 (1.0/M_LN10)
#define LG3 (log(3.0))
#define LG50 (log(0.5))
#define LG10 (log(0.1))
#define LG90 (log(0.9))
#define LGALPHA (log(ALPHA))

/* Command line arguments */
static char *bed_file;
static char *bam_file;
static char *fa_file;
static char *out_file;
static int nthread;
static int pao;
static int isc;
static int nodup;
static int splice;
static int verbose;
static int phred64;
static int bisulfite;
static int const_qual;
static int debug;
static double ref_prior, alt_prior, het_prior;
static double mut_prior, nomut_prior;

/* Time info */
static time_t now; 
static struct tm *time_info; 
#define print_status(M, ...) time(&now); time_info = localtime(&now); fprintf(stderr, M, ##__VA_ARGS__);

static char NT[4] = "ATGC";

KHASH_MAP_INIT_STR(rsh, vector_t)   // hashmap: string key, vector value
static khash_t(rsh) *refseq_hash; // pointer to hashmap
static pthread_mutex_t refseq_lock; 

vector_t *bed_read(FILE *file) {
    vector_t *reg_list = vector_create(64, REGION_T);

    char *line = NULL;
    ssize_t read_file = 0;
    size_t line_length = 0;
    while ((read_file = getline(&line, &line_length, file)) != -1) {
        if (line_length <= 0 || line[strspn(line, " \t\v\r\n")] == '\0') continue; // blank line
        if (line[0] == '#') continue;

        int pos1, pos2;
        char chr[line_length];
        int t = sscanf(line, "%s %d %d", chr, &pos1, &pos2);
        if (t < 3) { exit_err("bad fields in BED file\n"); }

        region_t *g = region_create(chr, pos1, pos2);
        vector_add(reg_list, g);
    }
    free(line); line = NULL;
    fclose(file);
    qsort(reg_list->data, reg_list->len, sizeof (void *), nat_sort_region);
    return reg_list;
}

static vector_t *bam_fetch(const char *bam_file, const char *chr, const int pos1, const int pos2) {
    /* Reads in region coordinates */
    vector_t *read_list = vector_create(64, READ_T);

    samFile *sam_in = sam_open(bam_file, "r"); // open bam file
    if (sam_in == NULL) { exit_err("failed to open BAM file %s\n", bam_file); }
    bam_hdr_t *bam_header = sam_hdr_read(sam_in); // bam header
    if (bam_header == 0) { exit_err("bad header %s\n", bam_file); }
    hts_idx_t *bam_idx = sam_index_load(sam_in, bam_file); // bam index
    if (bam_idx == NULL) { exit_err("failed to open BAM index %s\n", bam_file); }

    int tid = bam_name2id(bam_header, chr);
    hts_itr_t *iter = sam_itr_queryi(bam_idx, tid, pos1-1, pos2); // read iterator
    if (iter != NULL) {
        bam1_t *aln = bam_init1(); // initialize an alignment
        while (sam_itr_next(sam_in, iter, aln) >= 0) {
            read_t *read = read_fetch(bam_header, aln, pao, isc, nodup, splice, phred64, const_qual);
            if (read != NULL) vector_add(read_list, read);
        }
        bam_destroy1(aln);
    }
    hts_itr_destroy(iter);
    hts_idx_destroy(bam_idx);
    bam_hdr_destroy(bam_header);
    sam_close(sam_in);
    return read_list;
}

static fasta_t *refseq_fetch(char *name, const char *fa_file) {
    pthread_mutex_lock(&refseq_lock);
    int i;
	khiter_t k = kh_get(rsh, refseq_hash, name);
    if (k != kh_end(refseq_hash)) {
        vector_t *node = &kh_val(refseq_hash, k);
        fasta_t **f = (fasta_t **)node->data;
        for (i = 0; i < node->len; i++) {
            if (strcmp(f[i]->name, name) == 0) {
                pthread_mutex_unlock(&refseq_lock);
                return f[i];
            }
        }
        exit_err("failed to find %s in hash key %d\n", name, k);
    }

    faidx_t *fai = fai_load(fa_file);
    if (fai == NULL) { 
        errno = fai_build(fa_file);
        if (errno == 0) { fai = fai_load(fa_file); }
        else { exit_err("failed to build and open FA index %s\n", fa_file); }
    }
    if (!faidx_has_seq(fai, name)) { exit_err("failed to find %s in reference %s\n", name, fa_file); }

    fasta_t *f = fasta_create(name);
    //f->seq = fai_fetch(fai, name, &f->seq_length);
    f->seq = faidx_fetch_seq(fai, f->name, 0, faidx_seq_len(fai, f->name) - 1, &f->seq_length);
    char *s;
    for (s = f->seq; *s != '\0'; s++) *s = toupper(*s);

    int absent;
    k = kh_put(rsh, refseq_hash, f->name, &absent);
    vector_t *node = &kh_val(refseq_hash, k);
    if (absent) vector_init(node, 8, FASTA_T);
    vector_add(node, f);
    fai_destroy(fai);
    pthread_mutex_unlock(&refseq_lock);
    return f;
}

static inline void calc_prob_snps_mut_region(double *prgu, double *prgv, int g_pos, const double *matrix, int read_length, const char *seq, int seq_length, int pos, int start, int end, int *seqnt_map) {
    if (start < 0) start = 0;
    if (end >= seq_length) end = seq_length;

    int i, k;
    double prgu_i[end - start], prgv_i[end - start];
    for (i = start; i < end; i++) {
        int n = i - start;
        prgu_i[n] = calc_read_prob(matrix, read_length, seq, seq_length, i, seqnt_map); // reference probability per position i
        prgv_i[n] = prgu_i[n]; // alternative probability per position i

        int r_pos = g_pos - pos;
        if (r_pos >= 0 && r_pos < read_length) {
            int x = seq[g_pos] - 'A';
            if (x < 0 || x >= 26) { exit_err("Ref character %c at pos %d (%d) not in valid alphabet\n", seq[g_pos], g_pos, seq_length); }

            double probability = 0;
            for (k = 0; k < 4; k++) {
                if (NT[k] != seq[g_pos]) {
                    double p = matrix[read_length * seqnt_map[NT[k] - 'A'] + r_pos];
                    probability = (probability == 0) ? p : log_add_exp(probability, p);
                    //printf("%c %f\t", NT[k], p);
                }
            }
            prgv_i[n] = prgv_i[n] - matrix[read_length * seqnt_map[x] + r_pos] + probability; // update alternative array
            //printf("%d\t%d\t%c\t%d\t%f\t%f\t%f\t%f\n", i, g_pos, seq[g_pos], r_pos, prgu_i[n], prgv_i[n], (double)matrix[read_length * seqnt_map[x] + r_pos], probability);
        }
    }
    *prgu += log_sum_exp(prgu_i, end - start);
    *prgv += log_sum_exp(prgv_i, end - start);
}

static inline void calc_prob_snps_mut(double *prgu, double *prgv, int g_pos, const double *matrix, int read_length, const char *seq, int seq_length, int pos, int *splice_pos, int *splice_offset, int n_splice, int *seqnt_map) {
    /* Get the sequence g in G and its neighborhood (half a read length flanking regions) */
    int start = pos; // - (read_length / 2);
    int end = pos + 1; //(read_length / 2);

    *prgu = 0;
    *prgv = 0;

    int i, j;
    if (n_splice == 0) {
        calc_prob_snps_mut_region(prgu, prgv, g_pos, matrix, read_length, seq, seq_length, pos, start, end, seqnt_map);
    }
    else { // calculate the probability for each splice section separately
        int r_pos = 0;
        int g_pos = pos;
        for (i = 0; i <= n_splice; i++) {
            int r_len = (i < n_splice) ? splice_pos[i] - r_pos + 1 : read_length - r_pos;
            start = g_pos - (r_len / 2);
            end = g_pos + (r_len / 2);

            double *submatrix = malloc(NT_CODES * r_len * sizeof (double));
            for (j = 0; j < NT_CODES; j++) memcpy(&submatrix[r_len * j], &matrix[read_length * j + r_pos], r_len * sizeof (double));
            calc_prob_snps_mut_region(prgu, prgv, g_pos, submatrix, r_len, seq, seq_length, pos, start, end, seqnt_map);
            free(submatrix); submatrix = NULL;

            g_pos += r_len + splice_offset[i];
            r_pos = splice_pos[i] + 1;
        }
    }
}

static char *evaluate_nomutation(const region_t *g) {
    int i, readi;

    /* Reference sequence */
    fasta_t *f = refseq_fetch(g->chr, fa_file);
    if (f == NULL) return NULL;
    char *refseq = f->seq;
    int refseq_length = f->seq_length;

    /* Reads in region coordinates */
    vector_t *read_list = bam_fetch(bam_file, g->chr, g->pos1, g->pos2);
    if (read_list->len == 0) {
        free(read_list); read_list = NULL;
        return NULL;
    }
    read_t **read_data = (read_t **)read_list->data;

    int g_pos;
    double alt_probability = -DBL_MAX;
    double ref_probability = 0;
    int ref_count = 0;
    int alt_count = 0;
    int start = g->pos1;
    int end = g->pos2;
    if (start < 0) start = 0;
    if (end >= refseq_length) end = refseq_length;
    for (g_pos = start; g_pos <= end; g_pos++) {
        double ref = 0;
        double alt = 0;
        int r_count = 0;
        int a_count = 0;

        /* Aligned reads */
        for (readi = 0; readi < read_list->len; readi++) {
            //if (read_data[readi]->pos < g->pos1 || read_data[readi]->pos > g->pos2) continue;

            double is_match[read_data[readi]->length], no_match[read_data[readi]->length];
            for (i = 0; i < read_data[readi]->length; i++) {
                is_match[i] = p_match[read_data[readi]->qual[i]];
                no_match[i] = p_mismatch[read_data[readi]->qual[i]];
            }
            double readprobmatrix[NT_CODES * read_data[readi]->length];
            set_prob_matrix(readprobmatrix, read_data[readi], is_match, no_match, seqnt_map, bisulfite);

            double prgu, prgv;
            //printf("\t%d\t%s\t%d\t%d\t%s\n", g_pos, read_data[readi]->name, read_data[readi]->pos, read_data[readi]->length, read_data[readi]->qseq);
            calc_prob_snps_mut(&prgu, &prgv, g_pos, readprobmatrix, read_data[readi]->length, refseq, refseq_length, read_data[readi]->pos, read_data[readi]->splice_pos, read_data[readi]->splice_offset, read_data[readi]->n_splice, seqnt_map);
            //printf("%f\t%f\n\n", prgu, prgv);

            /* Multi-map alignments from XA tags: chr8,+42860367,97M3S,3;chr9,-44165038,100M,4; */
            if (read_data[readi]->multimapXA != NULL) {
                int xa_pos, n;
                char *s, xa_chr[strlen(read_data[readi]->multimapXA) + 1];
                for (s = read_data[readi]->multimapXA; sscanf(s, "%[^,],%d,%*[^;]%n", xa_chr, &xa_pos, &n) == 2; s += n + 1) {
                    if (strcmp(xa_chr, read_data[readi]->chr) != 0 && abs(xa_pos - read_data[readi]->pos) < read_data[readi]->length) { // if secondary alignment does not overlap primary aligment
                        fasta_t *f = refseq_fetch(xa_chr, fa_file);
                        if (f == NULL) continue;
                        char *xa_refseq = f->seq;
                        int xa_refseq_length = f->seq_length;

                        double *p_readprobmatrix = readprobmatrix;
                        double *newreadprobmatrix = NULL;
                        if ((xa_pos < 0 && !read_data[readi]->is_reverse) || (xa_pos > 0 && read_data[readi]->is_reverse)) { // opposite of primary alignment strand
                            newreadprobmatrix = reverse(readprobmatrix, read_data[readi]->length * NT_CODES);
                            p_readprobmatrix = newreadprobmatrix;
                        }

                        xa_pos = abs(xa_pos);
                        double readprobability = calc_prob(p_readprobmatrix, read_data[readi]->length, xa_refseq, xa_refseq_length, xa_pos, read_data[readi]->splice_pos, read_data[readi]->splice_offset, read_data[readi]->n_splice, seqnt_map);
                        prgu = log_add_exp(prgu, readprobability);
                        prgv = log_add_exp(prgv, readprobability);
                        free(newreadprobmatrix); newreadprobmatrix = NULL;
                    }
                    if (*(s + n) != ';') break;
                }
            }
            else if (read_data[readi]->multimapNH > 1) { // scale by the number of multimap positions
                double n = log(read_data[readi]->multimapNH - 1);
                double readprobability = prgu + n;
                prgu = log_add_exp(prgu, readprobability);
                prgv = log_add_exp(prgv, readprobability);
            }

            if (prgu == 0 && prgv == 0) continue;

            if (prgu > prgv && prgu - prgv > 0.69) r_count += 1;
            else if (prgv > prgu && prgv - prgu > 0.69) a_count += 1;

            double phet   = log_add_exp(LG50 + prgv, LG50 + prgu);
            double phet10 = log_add_exp(LG10 + prgv, LG90 + prgu);
            double phet90 = log_add_exp(LG90 + prgv, LG10 + prgu);
            if (phet10 > phet) phet = phet10;
            if (phet90 > phet) phet = phet90;

            /* Priors */
            ref += prgu + ref_prior;
            alt += log_add_exp(prgv + alt_prior, phet + het_prior);

            if (debug > 1) {
                fprintf(stderr, "::\t%d\t%s\t%d\t%f\t%f\t%f\t%d\t%d\t", g_pos, read_data[readi]->name, read_data[readi]->pos, prgu, prgv, prgu-prgv, r_count, a_count);
                for (i = 0; i < read_data[readi]->n_cigar; i++) fprintf(stderr, "%d%c ", read_data[readi]->cigar_oplen[i], read_data[readi]->cigar_opchr[i]);
                fprintf(stderr, "\n");
            }
        }
        ref_probability = ref;
        alt_probability = (alt_probability == 0) ? alt : log_add_exp(alt_probability, alt);
        if (a_count > alt_count) {
            ref_count = r_count;
            alt_count = a_count;
        }
        if (debug > 0) fprintf(stderr, "++\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\n", g->pos1, g->pos2, g_pos, (int)read_list->len, ref_count, alt_count, ref, alt, ref - alt, alt_probability);
    }

    ref_probability += nomut_prior;
    alt_probability += mut_prior;
    double odds = (ref_probability - alt_probability) * M_1_LN10;
    int n = snprintf(NULL, 0, "%s\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\n", g->chr, g->pos1, g->pos2, (int)read_list->len, ref_count, alt_count, ref_probability, alt_probability, odds) + 1;
    char *output = malloc(n * sizeof (*output));
    snprintf(output, n, "%s\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\n", g->chr, g->pos1, g->pos2, (int)read_list->len, ref_count, alt_count, ref_probability, alt_probability, odds);

    vector_destroy(read_list); free(read_list); read_list = NULL;
    return output;
}

typedef struct {
    vector_t *queue, *results;
    pthread_mutex_t q_lock;
    pthread_mutex_t r_lock;
} work_t;

static void *pool(void *work) {
    work_t *w = (work_t *)work;
    vector_t *queue = (vector_t *)w->queue;
    vector_t *results = (vector_t *)w->results;

    while (1) { //pthread_t ptid = pthread_self(); uint64_t threadid = 0; memcpy(&threadid, &ptid, min(sizeof (threadid), sizeof (ptid)));
        pthread_mutex_lock(&w->q_lock);
        region_t *g = (region_t *)vector_pop(queue);
        pthread_mutex_unlock(&w->q_lock);
        if (g == NULL) break;
        
        char *outstr = evaluate_nomutation(g);
        if (outstr == NULL) continue;

        pthread_mutex_lock(&w->r_lock);
        vector_add(results, outstr);
        pthread_mutex_unlock(&w->r_lock);
    }
    return NULL;
}

static void process(const vector_t *reg_list, FILE *out_fh) {
    int i;

    region_t **reg_data = (region_t **)reg_list->data;
    int nregions = reg_list->len;

    print_status("# Options: pao=%d isc=%d nodup=%d splice=%d bs=%d phred64=%d cq=%d\n", pao, isc, nodup, splice, bisulfite, phred64, const_qual);
    print_status("# Start: %d threads \t%s\t%s", nthread, bam_file, asctime(time_info));

    vector_t *queue = vector_create(nregions, VOID_T);
    vector_t *results = vector_create(nregions, VOID_T);
    for (i = 0; i < nregions; i++) {
        vector_add(queue, reg_data[i]);
    }
    work_t *w = malloc(sizeof (work_t));
    w->queue = queue;
    w->results = results;

    pthread_mutex_init(&w->q_lock, NULL);
    pthread_mutex_init(&w->r_lock, NULL);

    pthread_t tid[nthread];
    for (i = 0; i < nthread; i++) pthread_create(&tid[i], NULL, pool, w);
    for (i = 0; i < nthread; i++) pthread_join(tid[i], NULL);

    pthread_mutex_destroy(&w->q_lock);
    pthread_mutex_destroy(&w->r_lock);

    free(w); w = NULL;

    qsort(results->data, results->len, sizeof (void *), nat_sort_vector);
    fprintf(out_fh, "# CHR\tPOS1\tPOS2\tReads\tRefCount\tAltCount\tRefProb\tAltProb\tOdds\n");
    for (i = 0; i < results->len; i++) fprintf(out_fh, "%s", (char *)results->data[i]);
    vector_destroy(queue); free(queue); queue = NULL;
    vector_destroy(results); free(results); results = NULL;
    print_status("# Done:\t%s\t%s", bam_file, asctime(time_info));
}

static void print_usage() {
    printf("\nUsage: eagle [options] -v regions.bed -a alignment.bam -r reference.fasta\n\n");
    printf("Required:\n");
    printf("  -v --bed        FILE   Genome regions, BED file. [stdin]\n");
    printf("  -a --bam        FILE   Alignment data bam files, ref-coord sorted with bai index file.\n");
    printf("  -r --ref        FILE   Reference sequence, fasta file with fai index file.\n");
    printf("Options:\n");
    printf("  -o --out        FILE   Output file. [stdout]\n");
    printf("  -t --nthread    INT    Number of threads. [1]\n");
    printf("     --pao               Primary alignments only.\n");
    printf("     --isc               Ignore soft-clipped bases.\n");
    printf("     --nodup             Ignore marked duplicate reads (based on SAM flag).\n");
    printf("     --splice            RNA-seq spliced reads.\n");
    printf("     --bs         INT    Bisulfite treated reads. 0: off, 1: top/forward strand, 2: bottom/reverse strand, 3: both. [0]\n");
    printf("     --mut_prior  FLOAT  Prior probability for a mutation at any given reference position [0.001].\n");
    printf("     --verbose           Verbose mode, output likelihoods for each read seen for each hypothesis to stderr.\n");
    printf("     --phred64           Read quality scores are in phred64.\n");
    printf("     --cq         INT    Constant quality as a phred score, ignoring the quality field in SAM. [0 is off]\n");
}

int main(int argc, char **argv) {
    /* Command line parameters defaults */
    bed_file = NULL;
    bam_file = NULL;
    fa_file = NULL;
    out_file = NULL;
    nthread = 1;
    pao = 0;
    isc = 0;
    nodup = 0;
    verbose = 0;
    phred64 = 0;
    bisulfite = 0;
    const_qual = 0;
    debug = 0;

    mut_prior = 0.001;

    static struct option long_options[] = {
        {"bed", required_argument, NULL, 'v'},
        {"bam", required_argument, NULL, 'a'},
        {"ref", required_argument, NULL, 'r'},
        {"out", optional_argument, NULL, 'o'},
        {"nthread", optional_argument, NULL, 't'},
        {"pao", no_argument, &pao, 1},
        {"isc", no_argument, &isc, 1},
        {"nodup", no_argument, &nodup, 1},
        {"splice", no_argument, &splice, 1},
        {"verbose", no_argument, &verbose, 1},
        {"phred64", no_argument, &phred64, 1},
        {"debug", optional_argument, NULL, 'd'},
        {"mut_prior", optional_argument, NULL, 990},
        {"bs", optional_argument, NULL, 992},
        {"cq", optional_argument, NULL, 993},
        {0, 0, 0, 0}
    };

    int opt = 0;
    while ((opt = getopt_long(argc, argv, "v:a:r:o:t:d:", long_options, &opt)) != -1) {
        switch (opt) {
            case 0: 
                //if (long_options[option_index].flag != 0) break;
                break;
            case 'v': bed_file = optarg; break;
            case 'a': bam_file = optarg; break;
            case 'r': fa_file= optarg; break;
            case 'o': out_file = optarg; break;
            case 't': nthread = parse_int(optarg); break;
            case 'd': debug = parse_int(optarg); break;
            case 990: mut_prior = parse_double(optarg); break;
            case 992: bisulfite = parse_int(optarg); break;
            case 993: const_qual = parse_int(optarg); break;
            default: exit_usage("Bad options");
        }
    }
    if (optind > argc) { exit_usage("Bad program call"); }

    FILE *bed_fh = stdin;
    if (bed_file != NULL) { // default bed file handle is stdin unless a bed file option is used
        bed_fh = fopen(bed_file, "r");
        if (bed_fh == NULL) { exit_err("failed to open BED file %s\n", bed_file); }
    }
    else {
        bed_file = "stdin";
    }
    if (bam_file == NULL) { exit_usage("Missing alignments given as BAM file!"); } 
    if (fa_file == NULL) { exit_usage("Missing reference genome given as Fasta file!"); }
    if (nthread < 1) nthread = 1;
    if (mut_prior < 0 || mut_prior > 1) mut_prior = 0.001;
    nomut_prior = 1 - mut_prior;

    ref_prior = log(0.5);
    alt_prior = log(0.25);
    het_prior = log(0.25);

    FILE *out_fh = stdout;
    if (out_file != NULL) out_fh = fopen(out_file, "w"); // default output file handle is stdout unless output file option is used

    init_seqnt_map(seqnt_map);
    mut_prior = log(mut_prior);
    nomut_prior = log(nomut_prior);
    init_q2p_table(p_match, p_mismatch, 50);
    
    /* Start processing data */
    clock_t tic = clock();
    vector_t *reg_list = bed_read(bed_fh);
    print_status("# Read BED: %s\t%i entries\t%s", bed_file, (int)reg_list->len, asctime(time_info));

    refseq_hash = kh_init(rsh);

    pthread_mutex_init(&refseq_lock, NULL);
    process(reg_list, out_fh);
    if (out_file != NULL) fclose(out_fh);
    else fflush(stdout);
    pthread_mutex_destroy(&refseq_lock);

    khiter_t k;
    for (k = kh_begin(refseq_hash); k != kh_end(refseq_hash); k++) {
        if (kh_exist(refseq_hash, k)) vector_destroy(&kh_val(refseq_hash, k));
    }
    kh_destroy(rsh, refseq_hash);
    vector_destroy(reg_list); free(reg_list); reg_list = NULL;

    clock_t toc = clock();
    print_status("# CPU time (hr):\t%f\n", (double)(toc - tic) / CLOCKS_PER_SEC / 3600);

    return 0;
}
