/*
EAGLE: explicit alternative genome likelihood evaluator
Given the sequencing data and candidate variant, explicitly test 
the alternative hypothesis against the reference hypothesis

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
#include "heap.h"

/* Constants */
#define VERSION "1.1.1"
#define ALPHA 1.3     // Factor to account for longer read lengths lowering the probability a sequence matching an outside paralogous source

/* Precalculated log values */
#define M_1_LOG10E (1.0/M_LOG10E)
#define M_1_LN10 (1.0/M_LN10)
#define LOG50 (log(0.5))
#define LOG10 (log(0.1))
#define LOG90 (log(0.9))
#define LGALPHA (log(ALPHA))

/* Command line arguments */
static int debug;
static char *vcf_file;
static char *bam_file;
static char *fa_file;
static char *out_file;
static int nthread;
static int sharedr;
static int distlim;
static int maxdist;
static int maxh;
static int mvh;
static int pao;
static int isc;
static int nodup;
static int splice;
static int verbose;
static int lowmem;
static int phred64;
static int bisulfite;
static int const_qual;
static double hetbias;
static double omega, lgomega;
static int dp, gap_op, gap_ex;
static int rc;
static double ref_prior, alt_prior, het_prior;

/* Time info */
static time_t now; 
static struct tm *time_info; 
#define print_status(M, ...) time(&now); time_info = localtime(&now); fprintf(stderr, M, ##__VA_ARGS__);

KHASH_MAP_INIT_STR(rsh, fasta_t *)   // hashmap: string key, vector value
static khash_t(rsh) *refseq_hash; // pointer to hashmap
static pthread_mutex_t refseq_lock; 

static vector_t *vcf_read(FILE *file) {
    vector_t *var_list = vector_create(8, VARIANT_T);

    char *line = NULL;
    ssize_t read_file = 0;
    size_t line_length = 0;
    while ((read_file = getline(&line, &line_length, file)) != -1) {
        if (line_length <= 0 || line[strspn(line, " \t\v\r\n")] == '\0') continue; // blank line
        if (line[0] == '#') continue;

        int pos;
        char chr[line_length], ref[line_length], alt[line_length];
        int t = sscanf(line, "%s %d %*[^\t] %s %s", chr, &pos, ref, alt);
        if (t < 4 || has_numbers(ref) || has_numbers(alt)) { exit_err("bad fields in VCF file\n%s\n", line); }

        int n1, n2;
        char *s1, *s2, ref_token[strlen(ref) + 1], alt_token[strlen(alt) + 1];
        for (s1 = ref; sscanf(s1, "%[^, ]%n", ref_token, &n1) == 1 || sscanf(s1, "%[-]%n", ref_token, &n1) == 1; s1 += n1 + 1) { // heterogenenous non-reference (comma-delimited) as separate entries
            for (s2 = alt; sscanf(s2, "%[^, ]%n", alt_token, &n2) == 1 || sscanf(s2, "%[-]%n", alt_token, &n2) == 1; s2 += n2 + 1) {
                if (alt_token[0] != '.' && alt_token[0] != '*' && strcmp("<*:DEL>", alt_token) != 0) {
                    variant_t *v = variant_create(chr, pos, ref_token, alt_token);
                    vector_add(var_list, v);
                }
                if (*(s2 + n2) != ',') break;
            }
            if (*(s1 + n1) != ',') break;
        }
    }
    free(line); line = NULL;
    fclose(file);
    qsort(var_list->data, var_list->len, sizeof (void *), nat_sort_variant);
    return var_list;
}

static int bam_fetch_last(const char *bam_file, const char *chr, const int pos1, const int pos2) {
    /* Reads in variant j = i + 1 region coordinates */
    samFile *sam_in = sam_open(bam_file, "r"); // open bam file
    if (sam_in == NULL) { exit_err("failed to open BAM file %s\n", bam_file); }
    bam_hdr_t *bam_header = sam_hdr_read(sam_in); // bam header
    if (bam_header == 0) { exit_err("bad header %s\n", bam_file); }
    hts_idx_t *bam_idx = sam_index_load(sam_in, bam_file); // bam index
    if (bam_idx == NULL) { exit_err("failed to open BAM index %s\n", bam_file); }

    int last = -1;
    int tid = bam_name2id(bam_header, chr);
    hts_itr_t *iter = sam_itr_queryi(bam_idx, tid, pos1-1, pos2); // read iterator
    if (iter != NULL) {
        bam1_t *aln = bam_init1(); // initialize an alignment
        while (sam_itr_next(sam_in, iter, aln) >= 0) {
            if (aln->core.tid < 0) continue; // not mapped
            last = aln->core.pos + aln->core.l_qseq;
        }
        bam_destroy1(aln);
    }
    hts_itr_destroy(iter);
    hts_idx_destroy(bam_idx);
    bam_hdr_destroy(bam_header);
    sam_close(sam_in);
    return(last);
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
	khiter_t k = kh_get(rsh, refseq_hash, name);
    if (k != kh_end(refseq_hash)) {
        pthread_mutex_unlock(&refseq_lock);
        return kh_val(refseq_hash, k);
    }

    faidx_t *fai = fai_load(fa_file);
    if (fai == NULL) { 
        errno = fai_build(fa_file);
        if (errno == 0) { fai = fai_load(fa_file); }
        else { exit_err("failed to build and open FA index %s\n", fa_file); }
    }
    if (!faidx_has_seq(fai, name)) { exit_err("failed to find %s in reference %s\n", name, fa_file); }

    fasta_t *f = fasta_create(name);
    //f->seq = fai_fetch(fai, f->name, &f->seq_length);
    f->seq = faidx_fetch_seq(fai, f->name, 0, faidx_seq_len(fai, f->name) - 1, &f->seq_length);
    char *s;
    for (s = f->seq; *s != '\0'; s++) *s = toupper(*s);

    int absent;
    k = kh_put(rsh, refseq_hash, f->name, &absent);
    if (absent) { kh_val(refseq_hash, k) = f; }
    else { exit_err("# refseq_hash collision: %s", asctime(time_info)); }
    fai_destroy(fai);
    pthread_mutex_unlock(&refseq_lock);
    return f;
}

static char *construct_altseq(const char *refseq, int refseq_length, const vector_int_t *combo, variant_t **var_data, int *altseq_length) {
    int i;
    int offset = 0;
    char *altseq = strdup(refseq);
    *altseq_length = refseq_length;
    for (i = 0; i < combo->len; i++) {
        variant_t *v = var_data[combo->data[i]];
        int pos = v->pos - 1 + offset;
        if (pos < 0 || pos > *altseq_length) { exit_err("Variant at %s:%d is out of bounds in reference\n", v->chr, v->pos); }

        char *var_ref, *var_alt;
        if (v->ref[0] == '-') { // account for "-" variant representations 
            var_ref = "";
            var_alt = v->alt;
        }
        else if (v->alt[0] == '-') { // account for "-" variant representations
            var_ref = v->ref;
            var_alt = "";
        }
        else {
            char *s1 = v->ref;
            char *s2 = v->alt;
            while (*s1 == *s2) { // account for and disregard common prefix in variant
                s1++;
                s2++;
                pos++;
            }
            var_ref = s1;
            var_alt = s2;
        }
        size_t var_ref_length = strlen(var_ref);
        size_t var_alt_length = strlen(var_alt);
        int delta = var_alt_length - var_ref_length;
        offset += delta;
        if (delta == 0) { // snps, equal length haplotypes
            memcpy(altseq + pos, var_alt, var_alt_length * sizeof (*var_alt));
        }
        else { // indels
            char *newalt = malloc((*altseq_length + delta + 1) * sizeof (*newalt));
            memcpy(newalt, altseq, pos * sizeof (*newalt));
            memcpy(newalt + pos, var_alt, var_alt_length * sizeof (*newalt));
            memcpy(newalt + pos + var_alt_length, altseq + pos + var_ref_length, (*altseq_length - pos - var_ref_length) * sizeof (*newalt));
            *altseq_length += delta;
            newalt[*altseq_length] = '\0';
            free(altseq); altseq = NULL;
            altseq = newalt;
        }
    }
    return altseq;
}

static inline int variant_find(const vector_int_t *a, int v) {
    int i = 0;
    int j = a->len - 1;
    int n = (i + j) / 2;
    while (i <= j) {
        if (v == a->data[n]) return n;
        if (v > a->data[n]) i = n + 1;
        else j = n - 1;
        n = (i + j) / 2;
    }
    return -1;
}

static inline void variant_print(char **output, const vector_t *var_set, int i, int nreads, int not_alt_count, int has_alt_count, double total, double has_alt, double not_alt) {
    variant_t **var_data = (variant_t **)var_set->data;

    double prob = (has_alt - total) * M_1_LN10;
    double odds = (has_alt - not_alt) * M_1_LN10;

    int n = snprintf(NULL, 0, "%s\t%d\t%s\t%s\t%d\t%d\t%d\t%e\t%f\t", var_data[i]->chr, var_data[i]->pos, var_data[i]->ref, var_data[i]->alt, nreads, not_alt_count, has_alt_count, prob, odds) + 1;
    char token[n];
    snprintf(token, n, "%s\t%d\t%s\t%s\t%d\t%d\t%d\t%e\t%f\t", var_data[i]->chr, var_data[i]->pos, var_data[i]->ref, var_data[i]->alt, nreads, not_alt_count, has_alt_count, prob, odds);
    str_resize(output, strlen(*output) + n);
    strcat(*output, token);

    str_resize(output, strlen(*output) + 2);
    strcat(*output, "[");
    if (var_set->len > 1) {
        for (i = 0; i < var_set->len; i++) {
            n = snprintf(NULL, 0, "%s,%d,%s,%s;", var_data[i]->chr, var_data[i]->pos, var_data[i]->ref, var_data[i]->alt) + 1;
            char token[n];
            snprintf(token, n, "%s,%d,%s,%s;", var_data[i]->chr, var_data[i]->pos, var_data[i]->ref, var_data[i]->alt);
            str_resize(output, strlen(*output) + n);
            strcat(*output, token);
        }
    }
    str_resize(output, strlen(*output) + 3);
    strcat(*output, "]\n");
}

static void calc_likelihood(stats_t *stat, vector_t *var_set, const char *refseq, const int refseq_length, read_t **read_data, const int nreads, int seti, int *seqnt_map) {
    size_t i, readi;
    stat->ref = 0;
    stat->alt = 0;
    stat->het = 0;
    stat->ref_count = 0;
    stat->alt_count = 0;
    stat->seen = 0;

    variant_t **var_data = (variant_t **)var_set->data;
    double log_nv = log((double)var_set->len);

    int has_indel = 0;
    if (!lowmem) {
        for (i = 0; i < stat->combo->len; i++) {
            variant_t *v = var_data[stat->combo->data[i]];
            if (v->ref[0] == '-' || v->alt[0] == '-' || strlen(v->ref) != strlen(v->alt)) {
                has_indel = 1;
                break;
            }
        }
    }

    /* Alternative sequence */
    int altseq_length = 0;
    char *altseq = NULL;
    if (has_indel || dp) altseq = construct_altseq(refseq, refseq_length, stat->combo, var_data, &altseq_length); 

    /* Aligned reads */
    for (readi = 0; readi < nreads; readi++) {
        if (read_data[readi]->pos > var_data[stat->combo->data[0]]->pos || read_data[readi]->end < var_data[stat->combo->data[stat->combo->len - 1]]->pos) { // read must cross all variants in current combo
            vector_double_add(stat->read_prgv, -DBL_MAX);
            continue; // read must cross all variants in current combo
        }
        stat->seen++;

        double is_match[read_data[readi]->length], no_match[read_data[readi]->length];
        for (i = 0; i < read_data[readi]->length; i++) {
            is_match[i] = p_match[read_data[readi]->qual[i]];
            no_match[i] = p_mismatch[read_data[readi]->qual[i]];
            if (dp) {
                double n = is_match[i] - 1;
                is_match[i] += 1 - n;
                no_match[i] += 1 - n;
            }
        }
        /* Read probability matrix */
        double readprobmatrix[NT_CODES * read_data[readi]->length];
        set_prob_matrix(readprobmatrix, read_data[readi], is_match, no_match, seqnt_map, bisulfite);

        /* Outside Paralog Exact Formuation: Probability that read is from an outside the reference paralogous "elsewhere", f in F.  Approximate the bulk of probability distribution P(r|f):
           a) perfect match = prod[ (1-e) ]
           b) hamming/edit distance 1 = prod[ (1-e) ] * sum[ (e/3) / (1-e) ]
           c) lengthfactor = alpha ^ (read length - expected read length). Length distribution, for reads with different lengths (hard clipped), where longer reads should have a relatively lower P(r|f):
        P(r|f) = (perfect + hamming_1) / lengthfactor */
        double delta[read_data[readi]->length];
        for (i = 0; i < read_data[readi]->length; i++) delta[i] = no_match[i] - is_match[i];
        double a = sum_d(is_match, read_data[readi]->length);
        double elsewhere = log_add_exp(a, a + log_sum_exp(delta, read_data[readi]->length)) - (LGALPHA * (read_data[readi]->length - read_data[readi]->inferred_length));

        double prgu, prgv;
        //for (i =0; i < stat->combo->len; i++) { variant_t *v = var_data[stat->combo->data[i]]; printf("%d;%s;%s;", v->pos, v->ref, v->alt); }
        //printf("\t%s\t%d\t%d\t%s\n", read_data[readi]->name, read_data[readi]->pos, read_data[readi]->length, read_data[readi]->qseq);
        if (dp) {
            prgu = calc_prob_dp(readprobmatrix, read_data[readi]->length, refseq, refseq_length, read_data[readi]->pos, read_data[readi]->splice_pos, read_data[readi]->splice_offset, read_data[readi]->n_splice, gap_op, gap_ex, seqnt_map);
            prgv = calc_prob_dp(readprobmatrix, read_data[readi]->length, altseq, altseq_length, read_data[readi]->pos, read_data[readi]->splice_pos, read_data[readi]->splice_offset, read_data[readi]->n_splice, gap_op, gap_ex, seqnt_map);
        }
        else if (has_indel) {
            prgu = calc_prob(readprobmatrix, read_data[readi]->length, refseq, refseq_length, read_data[readi]->pos, read_data[readi]->splice_pos, read_data[readi]->splice_offset, read_data[readi]->n_splice, seqnt_map);
            prgv = calc_prob(readprobmatrix, read_data[readi]->length, altseq, altseq_length, read_data[readi]->pos, read_data[readi]->splice_pos, read_data[readi]->splice_offset, read_data[readi]->n_splice, seqnt_map);
        }
        else {
            calc_prob_snps(&prgu, &prgv, stat->combo, var_data, readprobmatrix, read_data[readi]->length, refseq, refseq_length, read_data[readi]->pos, read_data[readi]->splice_pos, read_data[readi]->splice_offset, read_data[readi]->n_splice, seqnt_map);
        }
        //printf("%f\t%f\n\n", prgv, prgu);
        double pout = elsewhere;

        /* Multi-map alignments from XA tags: chr8,+42860367,97M3S,3;chr9,-44165038,100M,4; */
        if (read_data[readi]->multimapXA != NULL) {
            int xa_pos, n;
            char *s, xa_chr[strlen(read_data[readi]->multimapXA) + 1];
            for (s = read_data[readi]->multimapXA; sscanf(s, "%[^,],%d,%*[^;]%n", xa_chr, &xa_pos, &n) == 2; s += n + 1) {
                pout = log_add_exp(pout, elsewhere); // the more multi-mapped, the more likely it is the read is from elsewhere (paralogous), hence it scales (multiplied) with the number of multi-mapped locations
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
            pout = log_add_exp(pout, elsewhere + n);
            prgu = log_add_exp(prgu, readprobability);
            prgv = log_add_exp(prgv, readprobability);
        }

        /* Mixture model: probability that the read is from elsewhere, outside paralogous source */
        pout += lgomega;
        prgu = log_add_exp(pout, prgu);
        prgv = log_add_exp(pout, prgv);

        /* Track combination with highest variant likelihood */
        if (prgv > read_data[readi]->prgv) {
            read_data[readi]->index = seti;
            read_data[readi]->prgu = (float)prgu;
            read_data[readi]->prgv = (float)prgv;
            read_data[readi]->pout = (float)pout;
        }

        /* Mixture model: heterozygosity or heterogeneity as explicit allele frequency mu such that P(r|GuGv) = (mu)(P(r|Gv)) + (1-mu)(P(r|Gu)) */
        double phet   = log_add_exp(LOG50 + prgv, LOG50 + prgu);
        double phet10 = log_add_exp(LOG10 + prgv, LOG90 + prgu);
        double phet90 = log_add_exp(LOG90 + prgv, LOG10 + prgu);
        if (phet10 > phet) phet = phet10;
        if (phet90 > phet) phet = phet90;

        /* Priors */
        prgu += ref_prior;
        prgv += alt_prior - log_nv;
        phet += het_prior - log_nv;
        stat->ref += prgu;
        stat->alt += prgv;
        stat->het += phet;

        vector_double_add(stat->read_prgv, log_add_exp(prgv, phet));

        /* Read count incremented only when the difference in probability is not ambiguous, > ~log(2) difference and more likely than pout */
        if (prgv > prgu && prgv - prgu > 0.69 && prgv - pout > 0.69) stat->alt_count += 1;
        else if (prgu > prgv && prgu - prgv > 0.69 && prgu - pout > 0.69) stat->ref_count += 1;

        if (debug >= 2) {
            fprintf(stderr, "%f\t%f\t%f\t%f\t%d\t%d\t", prgu, phet, prgv, pout, stat->ref_count, stat->alt_count);
            fprintf(stderr, "%s\t%s\t%d\t%d\t", read_data[readi]->name, read_data[readi]->chr, read_data[readi]->pos, read_data[readi]->end);
            for (i = 0; i < read_data[readi]->n_cigar; i++) fprintf(stderr, "%d%c ", read_data[readi]->cigar_oplen[i], read_data[readi]->cigar_opchr[i]);
            fprintf(stderr, "\t");
            for (i = 0; i < stat->combo->len; i++) { variant_t *v = var_data[stat->combo->data[i]]; fprintf(stderr, "%s,%d,%s,%s;", v->chr, v->pos, v->ref, v->alt); }
            fprintf(stderr, "\t");
            if (read_data[readi]->multimapXA != NULL) fprintf(stderr, "%s\t", read_data[readi]->multimapXA);
            else fprintf(stderr, "%d\t", read_data[readi]->multimapNH);
            if (read_data[readi]->flag != NULL) fprintf(stderr, "%s\t", read_data[readi]->flag);
            fprintf(stderr, "%s\t", read_data[readi]->qseq);
            for (i = 0; i < read_data[readi]->length; i++) fprintf(stderr, "%d ", read_data[readi]->qual[i]);
            fprintf(stderr, "\n");
        }
    }
    stat->mut = log_add_exp(stat->alt, stat->het);
    free(altseq); altseq = NULL;
    if (debug >= 1) {
        fprintf(stderr, "==\t%f\t%f\t%f\t%d\t%d\t%d\t", stat->ref, stat->het, stat->alt, stat->ref_count, stat->alt_count, (int)nreads);
        for (i = 0; i < stat->combo->len; i++) { variant_t *v = var_data[stat->combo->data[i]]; fprintf(stderr, "%s,%d,%s,%s;", v->chr, v->pos, v->ref, v->alt); } fprintf(stderr, "\n");
    }
}

static char *evaluate(vector_t *var_set) {
    size_t i, readi, seti;

    variant_t **var_data = (variant_t **)var_set->data;

    /* Reference sequence */
    fasta_t *f = refseq_fetch(var_data[0]->chr, fa_file);
    if (f == NULL) return NULL;
    char *refseq = f->seq;
    int refseq_length = f->seq_length;

    /* Reads in variant region coordinates */
    vector_t *read_list = bam_fetch(bam_file, var_data[0]->chr, var_data[0]->pos, var_data[var_set->len - 1]->pos);
    if (read_list->len == 0) {
        vector_destroy(read_list); free(read_list); read_list = NULL;
        return NULL;
    }
    read_t **read_data = (read_t **)read_list->data;

    /* Variant combinations as a vector of vectors */
    vector_t *combo = powerset(var_set->len, maxh, 0);
    /*
    for (seti = 0; seti < combo->len; seti++) { // Print combinations
        fprintf(stderr, "%d\t", (int)seti); 
        for (i = 0; i < ((vector_int_t *)combo->data[seti])->len; i++) { fprintf(stderr, "%d;", ((vector_int_t *)combo->data[seti])->data[i]); } fprintf(stderr, "\t"); 
        for (i = 0; i < ((vector_int_t *)combo->data[seti])->len; i++) { variant_t *v = var_data[((vector_int_t *)combo->data[seti])->data[i]]; fprintf(stderr, "%s,%d,%s,%s;", v->chr, v->pos, v->ref, v->alt); } fprintf(stderr, "\n"); 
    }
    */

    vector_t *stats = vector_create(var_set->len + 1, STATS_T);

    for (seti = 0; seti < combo->len; seti++) { // all, singles
        stats_t *s = stats_create((vector_int_t *)combo->data[seti], read_list->len);
        calc_likelihood(s, var_set, refseq, refseq_length, read_data, read_list->len, seti, seqnt_map);
        vector_add(stats, s);
    }
    if (var_set->len > 1) { // doubles and beyond
        heap_t *h = heap_create(STATS_T);
        for (seti = 1; seti < combo->len; seti++) heap_push(h, ((stats_t *)stats->data[seti])->mut, stats->data[seti]);

        stats_t *s;
        while (s = heap_pop(h), s != NULL) {
            if ((int)stats->len - var_set->len - 1 >= maxh) break;

            vector_t *c = vector_create(8, VOID_T);
            derive_combo(c, s->combo, var_set->len);
            for (i = 0; i < c->len; i++) {
                stats_t *s = stats_create((vector_int_t *)c->data[i], read_list->len);
                calc_likelihood(s, var_set, refseq, refseq_length, read_data, read_list->len, stats->len, seqnt_map);
                vector_add(stats, s);
                heap_push(h, s->mut, s);
            }
            vector_free(c); //combos in stat so don't destroy
        }
        heap_free(h);
    }
    vector_free(combo); //combos in stat so don't destroy

    stats_t **stat = (stats_t **)stats->data;

    /* Heterozygous non-reference haplotypes as mixture model hypotheses */
    int c[stats->len];
    memset(c, 0, sizeof (c));
    for (readi = 0; readi < read_list->len; readi++) c[read_data[readi]->index]++; // combinations, based on best combination in each read

    vector_int_t *haplotypes = vector_int_create(stats->len);
    for (i = 0; i < stats->len; i++) {
        if ((double)c[i] / (double)read_list->len >= 0.1) vector_int_add(haplotypes, i); // relevant combination if read count >= 10% of reads seen
    }
    combo = vector_create(haplotypes->len, VOID_T);
    if (haplotypes->len > 1) combinations(combo, 2, haplotypes->len); // combination pairs

    vector_double_t *prhap = vector_double_create(combo->len);
    for (seti = 0; seti < combo->len; seti++) { // mixture model probabilities of combination pairs
        int x = haplotypes->data[((vector_int_t *)combo->data[seti])->data[0]];
        int y = haplotypes->data[((vector_int_t *)combo->data[seti])->data[1]];
        vector_double_add(prhap, 0);
        for (readi = 0; readi < read_list->len; readi++) {
            if (stat[x]->read_prgv->data[readi] == -DBL_MAX && stat[y]->read_prgv->data[readi] == -DBL_MAX) continue;
            double phet   = log_add_exp(LOG50 + stat[x]->read_prgv->data[readi], LOG50 + stat[y]->read_prgv->data[readi]);
            double phet10 = log_add_exp(LOG10 + stat[x]->read_prgv->data[readi], LOG90 + stat[y]->read_prgv->data[readi]);
            double phet90 = log_add_exp(LOG90 + stat[x]->read_prgv->data[readi], LOG10 + stat[y]->read_prgv->data[readi]);
            if (phet10 > phet) phet = phet10;
            if (phet90 > phet) phet = phet90;
            prhap->data[seti] += phet; // equal prior probability to ref since this assumes heterozygous non-reference variant
        }
    }
    if (debug >= 1) {
        for (seti = 0; seti < combo->len; seti++) {
            int x = haplotypes->data[((vector_int_t *)combo->data[seti])->data[0]];
            int y = haplotypes->data[((vector_int_t *)combo->data[seti])->data[1]];
            fprintf(stderr, "==\t%d, %d, %f\n", x, y, prhap->data[seti]);
        }
    }

    double total = log_add_exp(stat[0]->mut, stat[0]->ref);
    for (seti = 1; seti < stats->len; seti++) {
        total = log_add_exp(total, stat[seti]->mut);
        total = log_add_exp(total, stat[seti]->ref);
    }
    for (seti = 0; seti < combo->len; seti++) total = log_add_exp(total, prhap->data[seti]);

    char *output = malloc(sizeof (*output));
    output[0] = '\0';
    if (mvh) { /* Max likelihood variant hypothesis */
        size_t max_seti = 0;
        double r = stat[0]->mut - stat[0]->ref;
        double has_alt = stat[0]->mut;
        for (seti = 1; seti < stats->len; seti++) { 
            if (stat[seti]->mut - stat[seti]->ref > r) {
                r = stat[seti]->mut - stat[seti]->ref;
                has_alt = stat[seti]->mut;
                max_seti = seti;
            }
        }
        vector_t *v = vector_create(var_set->len, VARIANT_T);
        for (i = 0; i < stat[max_seti]->combo->len; i++) vector_add(v, var_data[stat[max_seti]->combo->data[i]]);
        variant_print(&output, v, 0, stat[max_seti]->seen, stat[max_seti]->ref_count, stat[max_seti]->alt_count, log_add_exp(total, stat[max_seti]->ref), has_alt, stat[max_seti]->ref);
        vector_free(v); //variants in var_list so don't destroy
    }
    else { /* Marginal probabilities & likelihood ratios*/
        for (i = 0; i < var_set->len; i++) {
            double has_alt = 0;
            double not_alt = 0;
            int acount = -1;
            int rcount = -1;
            int seen = -1;
            for (seti = 0; seti < stats->len; seti++) {
                if (variant_find(stat[seti]->combo, i) != -1) { // if variant is in this combination
                    has_alt = (has_alt == 0) ? stat[seti]->mut : log_add_exp(has_alt, stat[seti]->mut);
                    not_alt = (not_alt == 0) ? stat[seti]->ref : log_add_exp(not_alt, stat[seti]->ref);
                    if (stat[seti]->seen > seen) seen = stat[seti]->seen;
                    if (stat[seti]->alt_count > acount) {
                        acount = stat[seti]->alt_count;
                        rcount = stat[seti]->ref_count;
                    }
                }
                else {
                    not_alt = (not_alt == 0) ? stat[seti]->mut : log_add_exp(not_alt, stat[seti]->mut);
                }
            }
            for (seti = 0; seti < combo->len; seti++) {
                int x = haplotypes->data[((vector_int_t *)combo->data[seti])->data[0]];
                int y = haplotypes->data[((vector_int_t *)combo->data[seti])->data[1]];
                if (variant_find(stat[x]->combo, i) != -1 || variant_find(stat[y]->combo, i) != -1) has_alt = log_add_exp(has_alt, prhap->data[seti]);
                else not_alt = log_add_exp(not_alt, prhap->data[seti]);
            }
            variant_print(&output, var_set, i, seen, rcount, acount, total, has_alt, not_alt);
        }
    }

    if (verbose) {
        for (readi = 0; readi < read_list->len; readi++) {
            if (read_data[readi]->prgu == 0 && read_data[readi]->prgv == 0 && read_data[readi]->pout == 0) continue; // unprocessed read
            flockfile(stderr);
            fprintf(stderr, "%s\t%s\t%d\t", read_data[readi]->name, read_data[readi]->chr, read_data[readi]->pos);
            fprintf(stderr, "%f\t%f\t%f\t", read_data[readi]->prgu, read_data[readi]->prgv, read_data[readi]->pout);
            for (i = 0; i < read_data[readi]->n_cigar; i++) fprintf(stderr, "%d%c", read_data[readi]->cigar_oplen[i], read_data[readi]->cigar_opchr[i]);
            fprintf(stderr, "\t");
            if (read_data[readi]->multimapXA != NULL) fprintf(stderr, "%s\t", read_data[readi]->multimapXA);
            else fprintf(stderr, "%d\t", read_data[readi]->multimapNH);
            if (read_data[readi]->flag != NULL) fprintf(stderr, "%s\t", read_data[readi]->flag);
            else fprintf(stderr, "NONE\t");
            fprintf(stderr, "[");
            for (i = 0; i < stat[read_data[readi]->index]->combo->len; i++) {
                variant_t *v = var_data[stat[read_data[readi]->index]->combo->data[i]];
                fprintf(stderr, "%s,%d,%s,%s;", v->chr, v->pos, v->ref, v->alt);
            }
            fprintf(stderr, "]\n");
            funlockfile(stderr);
        }
    }

    for (i = 0; i < combo->len; i++) vector_int_free(combo->data[i]);
    vector_free(combo); //not destroyed because previously vector_int_free all elements
    vector_int_free(haplotypes);
    vector_double_free(prhap);
    vector_destroy(read_list); free(read_list); read_list = NULL;
    vector_destroy(stats); free(stats); stats = NULL;
    return output;
}

typedef struct {
    vector_t *queue, *results;
    pthread_mutex_t q_lock;
    pthread_mutex_t r_lock;
    size_t len;
} work_t;

static void *pool(void *work) {
    work_t *w = (work_t *)work;

    size_t n = w->len / 10;
    while (1) { //pthread_t ptid = pthread_self(); uint64_t threadid = 0; memcpy(&threadid, &ptid, min(sizeof (threadid), sizeof (ptid)));
        pthread_mutex_lock(&w->q_lock);
        vector_t *var_set = (vector_t *)vector_pop(w->queue);
        pthread_mutex_unlock(&w->q_lock);
        if (var_set == NULL) break;
        
        char *outstr = evaluate(var_set);
        if (outstr != NULL) {
            pthread_mutex_lock(&w->r_lock);
            if (!verbose && n > 10 && w->results->len > 10 && w->results->len % n == 0) {
                print_status("# Progress: %zd%%: %zd / %zd\t%s", 10 * w->results->len / n, w->results->len, w->queue->len, asctime(time_info));
            }
            vector_add(w->results, outstr);
            pthread_mutex_unlock(&w->r_lock);
        }
        vector_free(var_set); //variants in var_list so don't destroy
    }
    return NULL;
}

static void process(const vector_t *var_list, FILE *out_fh) {
    size_t i, j;

    variant_t **var_data = (variant_t **)var_list->data;

    i = 0;
    vector_t *var_set = vector_create(var_list->len, VOID_T);
    if (sharedr == 1) { /* Variants that share a read: shared with a given first variant */
        while (i < var_list->len) {
            vector_t *curr = vector_create(8, VARIANT_T);
            vector_add(curr, var_data[i]);

            /* Reads in variant i region coordinates */
            int i_last = bam_fetch_last(bam_file, var_data[i]->chr, var_data[i]->pos, var_data[i]->pos);

            j = i + 1;
            while (j < var_list->len && strcmp(var_data[i]->chr, var_data[j]->chr) == 0) { // while last read in i will reach j
                if (var_data[j]->pos > i_last) break;
                vector_add(curr, var_data[j++]);
            }
            i = j;
            vector_add(var_set, curr);
        }
    }
    else if (sharedr == 2) { /* Variants that share a read: shared with any neighboring variant */
        while (i < var_list->len) {
            vector_t *curr = vector_create(8, VARIANT_T);
            vector_add(curr, var_data[i]);

            j = i + 1;
            while (j < var_list->len && strcmp(var_data[i]->chr, var_data[j]->chr) == 0) { // while last read in i will reach j
                /* Reads in variant i region coordinates */
                int i_last = bam_fetch_last(bam_file, var_data[i]->chr, var_data[i]->pos, var_data[i]->pos);
                if (var_data[j]->pos > i_last) break;
                vector_add(curr, var_data[j]);
                i++;
                j++;
            }
            i = j;
            vector_add(var_set, curr);
        }
    }
    else { /* Variants that are close together as sets */
        while (i < var_list->len) {
            vector_t *curr = vector_create(8, VARIANT_T);
            vector_add(curr, var_data[i]);
            j = i + 1;
            while (distlim > 0 && j < var_list->len && strcmp(var_data[j]->chr, var_data[j - 1]->chr) == 0 && abs(var_data[j]->pos - var_data[j - 1]->pos) <= distlim) {
                if (maxdist > 0 && abs(var_data[j]->pos - var_data[i]->pos) > maxdist) break;
                vector_add(curr, var_data[j++]);
            }
            i = j;
            vector_add(var_set, curr);
        }
    }
    /* Heterozygous non-reference variants as separate entries */
    int flag_add = 1;
    while (flag_add) {
        flag_add = 0;
        for (i = 0; i < var_set->len; i++) {
            vector_t *curr_set = (vector_t *)var_set->data[i];
            if (curr_set->len == 1) continue;

            int flag_nonset = 1;
            for (j = 0; j < curr_set->len - 1; j++) { // check if all entries have the same position
                variant_t *curr = (variant_t *)curr_set->data[j];
                variant_t *next = (variant_t *)curr_set->data[j + 1];
                if (curr->pos == next->pos && strcmp(curr->chr, next->chr) == 0 && strcmp(curr->ref, next->ref) == 0 && strcmp(curr->alt, next->alt) == 0) vector_del(curr_set, j + 1); // delete duplicate entries
                else if (curr->pos != next->pos) flag_nonset = 0;
            }
            if (flag_nonset) { // only 1 entry, with multiple heterozygous non-reference variants
                while (curr_set->len > 1) {
                    variant_t *curr = (variant_t *)vector_pop(curr_set);
                    vector_t *dup = vector_create(8, VARIANT_T);
                    vector_add(dup, curr);
                    vector_add(var_set, dup);
                }
            }
            else { // multiple entries comprising a set
                for (j = 0; j < curr_set->len - 1; j++) {
                    variant_t *curr = (variant_t *)curr_set->data[j];
                    variant_t *next = (variant_t *)curr_set->data[j + 1];
                    if (curr->pos == next->pos) {
                        flag_add = 1;
                        vector_t *dup = vector_dup(curr_set);
                        vector_del(curr_set, j);
                        vector_del(dup, j + 1);
                        vector_add(var_set, dup);
                    }
                }
            }
        }
    } 
    if (sharedr == 1) { print_status("# Variants with shared reads to first in set: %i entries\t%s", (int)var_set->len, asctime(time_info)); }
    else if (sharedr == 2) { print_status("# Variants with shared reads to any in set: %i entries\t%s", (int)var_set->len, asctime(time_info)); }
    else { print_status("# Variants within %d (max window: %d) bp: %i entries\t%s", distlim, maxdist, (int)var_set->len, asctime(time_info)); }

    print_status("# Options: maxh=%d mvh=%d pao=%d isc=%d nodup=%d splice=%d bs=%d lowmem=%d phred64=%d\n", maxh, mvh, pao, isc, nodup, splice, bisulfite, lowmem, phred64);
    print_status("#          dp=%d gap_op=%d gap_ex=%d\n", dp, gap_op, gap_ex);
    print_status("#          hetbias=%g omega=%g cq=%d\n", hetbias, omega, const_qual);
    print_status("#          verbose=%d\n", verbose);
    print_status("# Start: %d threads \t%s\t%s", nthread, bam_file, asctime(time_info));

    vector_t *queue = vector_create(var_set->len, VOID_T);
    vector_t *results = vector_create(var_set->len, VOID_T);
    for (i = 0; i < var_set->len; i++) vector_add(queue, var_set->data[i]);

    work_t *w = malloc(sizeof (work_t));
    w->queue = queue;
    w->results = results;
    w->len = var_set->len;

    pthread_mutex_init(&w->q_lock, NULL);
    pthread_mutex_init(&w->r_lock, NULL);

    pthread_t tid[nthread];
    for (i = 0; i < nthread; i++) pthread_create(&tid[i], NULL, pool, w);
    for (i = 0; i < nthread; i++) pthread_join(tid[i], NULL);

    pthread_mutex_destroy(&w->q_lock);
    pthread_mutex_destroy(&w->r_lock);

    free(w); w = NULL;
    vector_free(var_set); //variants in var_list so don't destroy

    qsort(results->data, results->len, sizeof (void *), nat_sort_vector);
    fprintf(out_fh, "# SEQ\tPOS\tREF\tALT\tReads\tRefReads\tAltReads\tProb\tOdds\tSet\n");
    for (i = 0; i < results->len; i++) fprintf(out_fh, "%s", (char *)results->data[i]);
    vector_destroy(queue); free(queue); queue = NULL;
    vector_destroy(results); free(results); results = NULL;
    print_status("# Done:\t%s\t%s", bam_file, asctime(time_info));
}

static void print_usage() {
    printf("\nUsage: eagle [options] -v variants.vcf -a alignment.bam -r reference.fasta\n\n");
    printf("Required:\n");
    printf("  -v --vcf      FILE   Variants VCF file. [stdin]\n");
    printf("  -a --bam      FILE   Alignment data bam files, ref-coord sorted with bai index file.\n");
    printf("  -r --ref      FILE   Reference sequence, fasta file with fai index file.\n");
    printf("Options:\n");
    printf("  -o --out      FILE   Output file. [stdout]\n");
    printf("  -t --nthread  INT    Number of threads. [1]\n");
    printf("  -s --sharedr  INT    Group nearby variants that share a read, 0:distance based/off, 1:shared with first, 2:shared with any. [0]\n");
    printf("  -n --distlim  INT    Group nearby variants within n bases, 0:off. [10]\n");
    printf("  -w --maxdist  INT    Maximum number of bases between any two variants in a set of hypotheses, 0:off. [0]\n");
    printf("  -m --maxh     INT    Maximum number of combinations in the set of hypotheses, instead of all 2^n. [1024]\n");
    printf("     --mvh             Output the maximum likelihood hypothesis in the set instead of marginal probabilities.\n");
    printf("     --pao             Primary alignments only.\n");
    printf("     --isc             Ignore soft-clipped bases.\n");
    printf("     --nodup           Ignore marked duplicate reads (based on SAM flag).\n");
    printf("     --splice          RNA-seq spliced reads.\n");
    printf("     --bs       INT    Bisulfite treated reads. 0: off, 1: top/forward strand, 2: bottom/reverse strand, 3: both. [0]\n");
    printf("     --dp              Use dynamic programming to calculate likelihood instead of the basic model.\n");
    printf("     --gap_op   INT    DP gap open penalty. [6]. Recommend 2 for long reads with indel errors.\n");
    printf("     --gap_ex   INT    DP gap extend penalty. [1].\n");
    printf("     --verbose         Verbose mode, output likelihoods for each read seen for each hypothesis to stderr.\n");
    printf("     --lowmem          Low memory usage mode, the default mode for snps, this may be slightly slower for indels but uses less memory.\n");
    printf("     --phred64         Read quality scores are in phred64.\n");
    printf("     --hetbias  FLOAT  Prior probability bias towards non-homozygous mutations, between [0,1]. [0.5]\n");
    printf("     --omega    FLOAT  Prior probability of originating from outside paralogous source, between [0,1]. [1e-6]\n");
    printf("     --cq       INT    Constant quality as a phred score, ignoring the quality field in SAM. [0 is off]\n");
    printf("     --rc              Wrapper for read classification settings: --omega=1.0e-40 --isc --mvh --verbose --lowmem.\n");
    printf("     --version         Display version.\n");
}

int main(int argc, char **argv) {
    /* Command line parameters defaults */
    debug = 0;
    vcf_file = NULL;
    bam_file = NULL;
    fa_file = NULL;
    out_file = NULL;
    nthread = 1;
    sharedr = 0;
    distlim = 10;
    maxdist = 0;
    maxh = 1024;
    mvh = 0;
    pao = 0;
    isc = 0;
    nodup = 0;
    splice = 0;
    bisulfite = 0;
    verbose = 0;
    lowmem = 0;
    phred64 = 0;
    dp = 0;
    gap_op = 6;
    gap_ex = 1;
    hetbias = 0.5;
    omega = 1.0e-6;
    const_qual = 0;
    rc = 0;

    static struct option long_options[] = {
        {"debug", optional_argument, NULL, 'd'},
        {"vcf", required_argument, NULL, 'v'},
        {"bam", required_argument, NULL, 'a'},
        {"ref", required_argument, NULL, 'r'},
        {"out", optional_argument, NULL, 'o'},
        {"nthread", optional_argument, NULL, 't'},
        {"sharedr", optional_argument, NULL, 's'},
        {"distlim", optional_argument, NULL, 'n'},
        {"maxdist", optional_argument, NULL, 'w'},
        {"maxh", optional_argument, NULL, 'm'},
        {"maxh", optional_argument, NULL, 'm'},
        {"mvh", no_argument, &mvh, 1},
        {"pao", no_argument, &pao, 1},
        {"isc", no_argument, &isc, 1},
        {"nodup", no_argument, &nodup, 1},
        {"splice", no_argument, &splice, 1},
        {"verbose", no_argument, &verbose, 1},
        {"phred64", no_argument, &phred64, 1},
        {"lowmem", no_argument, &lowmem, 1},
        {"dp", no_argument, &dp, 1},
        {"gap_op", optional_argument, NULL, 981},
        {"gap_ex", optional_argument, NULL, 982},
        {"hetbias", optional_argument, NULL, 990},
        {"omega", optional_argument, NULL, 991},
        {"bs", optional_argument, NULL, 992},
        {"cq", optional_argument, NULL, 993},
        {"rc", no_argument, &rc, 1},
        {"version", optional_argument, NULL, 999},
        {0, 0, 0, 0}
    };

    int opt = 0;
    while ((opt = getopt_long(argc, argv, "d:v:a:r:o:t:s:n:w:m:", long_options, &opt)) != -1) {
        switch (opt) {
            case 0: 
                //if (long_options[option_index].flag != 0) break;
                break;
            case 'd': debug = parse_int(optarg); break;
            case 'v': vcf_file = optarg; break;
            case 'a': bam_file = optarg; break;
            case 'r': fa_file= optarg; break;
            case 'o': out_file = optarg; break;
            case 't': nthread = parse_int(optarg); break;
            case 's': sharedr = parse_int(optarg); break;
            case 'n': distlim = parse_int(optarg); break;
            case 'w': maxdist = parse_int(optarg); break;
            case 'm': maxh = parse_int(optarg); break;
            case 981: gap_op = parse_int(optarg); break;
            case 982: gap_ex = parse_int(optarg); break;
            case 990: hetbias = parse_double(optarg); break;
            case 991: omega = parse_double(optarg); break;
            case 992: bisulfite = parse_int(optarg); break;
            case 993: const_qual = parse_int(optarg); break;
            case 999: printf("EAGLE %s\n", VERSION); exit(0);
            default: exit_usage("Bad options");
        }
    }
    if (optind > argc) { exit_usage("Bad program call"); }

    FILE *vcf_fh = stdin;
    if (vcf_file != NULL) { // default vcf file handle is stdin unless a vcf file option is used
        vcf_fh = fopen(vcf_file, "r");
        if (vcf_fh == NULL) { exit_err("failed to open VCF file %s\n", vcf_file); }
    }
    else {
        vcf_file = "stdin";
    }
    if (bam_file == NULL) { exit_usage("Missing alignments given as BAM file!"); } 
    if (fa_file == NULL) { exit_usage("Missing reference genome given as Fasta file!"); }
    if (nthread < 1) nthread = 1;
    if (sharedr < 0 || sharedr > 2) sharedr = 0;
    if (distlim < 0) distlim = 10;
    if (maxdist < 0) maxdist = 0;
    if (maxh < 0) maxh = 0;
    if (gap_op <= 0) gap_op = 6;
    if (gap_ex <= 0) gap_ex = 1;
    if (hetbias < 0 || hetbias > 1) hetbias = 0.5;
    if (omega < 0 || omega > 1) omega = 1e-6;
    if (rc) {
        omega = 1e-40;
        isc = 1;
        mvh = 1;
        verbose = 1;
        lowmem = 1;
    }
    lgomega = (log(omega) - log(1.0-omega));

    ref_prior = log(0.5);
    alt_prior = log(0.5 * (1 - hetbias));
    het_prior = log(0.5 * hetbias);

    FILE *out_fh = stdout;
    if (out_file != NULL) out_fh = fopen(out_file, "w"); // default output file handle is stdout unless output file option is used

    init_seqnt_map(seqnt_map);
    init_q2p_table(p_match, p_mismatch, 50);

    /* Start processing data */
    clock_t tic = clock();
    vector_t *var_list = vcf_read(vcf_fh);
    print_status("# Read VCF: %s\t%i entries\t%s", vcf_file, (int)var_list->len, asctime(time_info));

    refseq_hash = kh_init(rsh);

    pthread_mutex_init(&refseq_lock, NULL);
    process(var_list, out_fh);
    if (out_file != NULL) fclose(out_fh);
    else fflush(stdout);
    pthread_mutex_destroy(&refseq_lock);

    khiter_t k;
    for (k = kh_begin(refseq_hash); k != kh_end(refseq_hash); k++) {
        if (kh_exist(refseq_hash, k)) {
            fasta_destroy(kh_val(refseq_hash, k)); free(kh_val(refseq_hash, k)); kh_val(refseq_hash, k) = NULL;
        }
    }
    kh_destroy(rsh, refseq_hash);
    vector_destroy(var_list); free(var_list); var_list = NULL;

    clock_t toc = clock();
    print_status("# CPU time (hr):\t%f\n", (double)(toc - tic) / CLOCKS_PER_SEC / 3600);

    return 0;
}
