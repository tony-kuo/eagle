#!/usr/bin/python

# Copyright 2016 Tony Kuo
# This program is distributed under the terms of the GNU General Public License

# EAGLE: explicit alternative genome likelihood evaluator
# Given the sequencing data, explicitly test the alternative genome hypothesis
# and calculate its probability against the reference genome hypothesis

from __future__ import print_function
import argparse, re, sys, os;
import pysam;
import numpy as np;
from cffi import FFI;
from itertools import chain, combinations;
from multiprocessing import Pool, Process, Manager;
from datetime import datetime;
from signal import signal, SIGPIPE, SIG_DFL;
signal(SIGPIPE,SIG_DFL);

# Constants
omega = 1E-4; # Prior probability of read originating from an outside paralogous source
alpha = 1.3; # Factor to account for longer read lengths lowering the probability a sequence matching an outside paralogous source

# Precalculated log values 
l3 = np.log10(3);
l50 = np.log10(0.5);
l10 = np.log10(0.1);
l90 = np.log10(0.9);
lomega = np.log10(omega);
l1_omega = np.log10(1-omega);
lalpha = np.log10(alpha);

complement = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a' };

ffi = FFI();
ffi.cdef ("""
void initAlphaMap(void);
double log10addexp(double, double);
double log10sumexp(double*, size_t);
void setReadProbMatrix(char*, int, double*, double*, double*);
double calcReadProb(char*, int, int, int, double*, double);
double calcReadProbability(char*, int, int, int, double*);
""")
C = ffi.verify ("""
static int alphabetval[26];
static void initAlphaMap(void) {
    memset(alphabetval, 4, sizeof(alphabetval)); // Default value 4 to array elements
    alphabetval['A'-'A'] = 0;
    alphabetval['T'-'A'] = 1;
    alphabetval['G'-'A'] = 2;
    alphabetval['C'-'A'] = 3;
    alphabetval['N'-'A'] = 4;
}
double log10addexp(double a, double b) {
    double max_exp;
    if ( a > b ) { max_exp = a; }
    else { max_exp = b; }
    return (log10(pow(10, a - max_exp) + pow(10, b - max_exp)) + max_exp);
}
double log10sumexp(double* n, size_t count) {
    size_t i;
    double max_exp = n[0]; 
    for ( i = 1; i < count; i++ ) { if ( n[i] > max_exp ) max_exp = n[i]; }
    double sum = 0.0;
    for ( i = 0; i < count; i++ ) { sum += pow(10, n[i] - max_exp); }
    return (log10(sum) + max_exp);
}
void setReadProbMatrix(char* seq, int readlength, double* isbase, double* notbase, double* matrix) {
    int i, b; // array[width * row + col] = value
    memset(matrix, 0, readlength*5*sizeof(double)); // Zero out array
    for ( b = 0; b < readlength; b++ ) {
        for ( i = 0; i < 5; i++ ) matrix[5*b+i] = notbase[b];
        matrix[5*b+alphabetval[toupper(seq[b])-'A']] = isbase[b];
    }
}
double calcReadProb(char* seq, int seqlength, int pos, int readlength, double* matrix, double baseline) {
    int b; // array[width * row + col] = value
    double probability = 0.0;
    //printf("%d\\t", pos);
    for ( b = pos;  b < pos+readlength; b++ ) {
        if ( b < 0) continue; // Skip if position is before the start of reference seq
        if ( b >= seqlength ) break; // Stop if it reaches the end of reference seq
        //printf("%c", seq[b]);
        probability += matrix[5*(b-pos)+alphabetval[toupper(seq[b])-'A']]; 
        if ( probability < baseline - 10 ) break; // Stop sum if less than 1% contribution to baseline (best, highest) probability mass
    }
    //printf("\\t%f\\n", probability);
    return (probability);
}
double calcReadProbability(char* seq, int seqlength, int pos, int readlength, double* matrix) {
    int i;
    double probability = 0.0;
    double baseline = calcReadProb(seq, seqlength, pos, readlength, matrix, -1000); // First probability at refpos, likely the highest, to be used as first best
    for ( i = pos-readlength; i <= pos+readlength; i++ ) {
        if ( i + readlength < 0 ) continue;
        if ( i >= seqlength ) break;

        if ( probability == 0.0 ) { probability = calcReadProb(seq, seqlength, i, readlength, matrix, baseline); }
        else { probability = log10addexp(probability, calcReadProb(seq, seqlength, i, readlength, matrix, baseline)); }

        if ( probability > baseline ) baseline = probability; // Update the best probability so far
    }
    return (probability);
}
""")
C.initAlphaMap(); # Initialize alphabet to int mapping table

def naturalSort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower();
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)];
    return sorted(l, key=alphanum_key);

def readVCF(filename):
    entry = {};
    with open(filename, 'r') as fh:
        for line in fh:
            if line[0] == '#': continue;
            var = line.strip().split('\t');
            pos = int(var[1]);
            ref = var[3].upper().split(',');
            alt = var[4].upper().split(',');
            # Account for double heterozygous non-reference or entries with the same position
            for i in ref:
                for j in alt:
                    entry[(var[0], pos, i, j)] = [(pos, i, j)]; # key = (chr, pos, ref, alt)
    fh.close;
    print("Read VCF:\t{0}\t{1} entries\t{2}".format(filename, len(entry.keys()), datetime.now()), file=sys.stderr);
    return(entry);

def groupNearbyVariants(entry):
    # Check if current variant is within distance limit of the previous on the same chromosome, merge if so
    varid = sorted(entry.keys());
    skip2ind = -1;
    for i in range(0, len(varid)-1):
        if i <= skip2ind: continue;
        for j in range(i+1, len(varid)):
            if varid[j][0] != varid[j-1][0] or varid[j][1] - varid[j-1][1] > distancethreshold: break;
            entry[varid[i]].extend(entry[varid[j]]);
            del entry[varid[j]];
            skip2ind = j;
    # Split heterozygous non-reference variants into separate entries
    while True:
        misc = {};
        varid = sorted(entry.keys());
        for i in range(0, len(varid)):
            if len(entry[varid[i]]) <= 1: continue; # Skip singletons
            for j in range(0, len(entry[varid[i]])-1):
                if entry[varid[i]][j][0] == entry[varid[i]][j+1][0]:
                    newid = (varid[i][0], varid[i][1], varid[i][2], varid[i][3], j);
                    misc[newid] = list(entry[varid[i]]);
                    del entry[varid[i]][j];
                    del misc[newid][j+1];
                    break;
        if misc: entry.update(misc);
        else: break;
    # Combine adjacent SNPs
    #misc = {};
    #varid = sorted(entry.keys());
    #for i in range(0, len(varid)):
        #if len(entry[varid[i]]) <= 1: continue; # Ignore singletons
        #x = max([entry[varid[i]][j+1][0] - entry[varid[i]][j][0] for j in range(0, len(entry[varid[i]])-1)]); # If all variants in set are adjacent to each other
        #y = max([max(len(j[1]), len(j[2])) for j in entry[varid[i]]]); # If all variants are length 1
        #z = [j[1] for j in entry[varid[i]]] + [j[2] for j in entry[varid[i]]]; # Account for '-' representation
        #if x == 1 and y == 1 and '-' not in z:
            #ref = ''.join([j[1] for j in entry[varid[i]]]);
            #alt = ''.join([j[2] for j in entry[varid[i]]]);
            #misc[(varid[i][0], varid[i][1], ref, alt)] = [(varid[i][1], ref, alt)]; # Merge and delete extraneous entries
            #del entry[varid[i]];
    #if misc: entry.update(misc);
    #for i in entry:
        #if len(entry[i]) > maxk: print(i, len(entry[i]))
    #sys.exit(0);
    print("Group within:\t{0} bp\t{1} entries \t{2}".format(distancethreshold, len(entry.keys()), datetime.now()), file=sys.stderr);
    return(entry);

def readFasta(filename):
    print("Read genome:\t{0}\t{1}".format(filename, datetime.now()), file=sys.stderr);
    seqid = "";
    seq = {};
    seqlength = {};
    with open(filename, 'r') as fh:
        for line in fh:
            line = line.strip();
            if line[0] == '>': #re.match('^>.+', line):
                seqid = re.split('>| ', line)[1]; # >chr1 1 -> ['', 'chr1', '1']
                seq[seqid] = [];
            else:
                seq[seqid].append(line.upper());
    fh.close;
    for seqid in seq: 
        seq[seqid] = ''.join(seq[seqid]).encode('utf-8'); # Explicit byte string for python3, needed for cffi char*
        seqlength[seqid] = len(seq[seqid]);
    return(seq, seqlength);

def processWork(inqueue, results):
    while not inqueue.empty(): 
        args = inqueue.get();
        results.extend(evaluateVariant(args));

def readPYSAM(files, var_list, outfile):
    results = [];

    #manager = Manager();
    #work = manager.Queue(numprocesses);
    #results = manager.list();
    for fn in files: 
        print("Start:\t{0}\t{1}".format(fn, datetime.now()), file=sys.stderr);
        varid = sorted(list(var_list.keys()));

        # Testing speed with Process, seems slower than Pool
        #pool = [];
        #for i in range(0,numprocesses):
            #p = Process(target=processWork, args=(work, results));
            #p.start()
            #pool.append(p);
        #for n in range(0, len(varid)): 
            #if varid[n][0] in refseq: work.put((fn, varid[n], var_list[varid[n]])); 
        #for p in pool: p.join();

        args = [];
        for n in range(0, len(varid)): 
            if varid[n][0] in refseq: 
                args.append((fn, varid[n], var_list[varid[n]]));
                #results.extend(evaluateVariant((fn, varid[n], var_list[varid[n]]))); # single process, for testing
        #continue;
        try:
            pool = Pool(processes=numprocesses);
            poolresults = pool.map_async(evaluateVariant, args);
            for p in poolresults.get(): results.extend(p);
        finally:
            pool.close();
            pool.join();

    if len(outfile) > 0: fh = open(outfile, 'w');
    else: fh = sys.stdout;
    print('#SEQ\tPOS\tREF\tALT\tReads\tAltReads\tlog10_Prob\tlog10_Odds\tVarSet', file=fh);
    for i in naturalSort(results): print(i, file=fh);
    fh.close();
    print("Done:\t{0}\t{1}".format(fn, datetime.now()), file=sys.stderr);

def evaluateVariant(args):
    (fn, varid, var_set) = args;
    refentry = {};
    altentry = {};
    readentry = {};
    varstart = min([a[0] for a in var_set]);
    varend = max([a[0] for a in var_set]);

    # Fetch reads that cross the variant location
    samfile = pysam.AlignmentFile(fn, "rb");
    readset = list(samfile.fetch(reference=varid[0], start=max(0,varstart-1), end=min(varend+1,reflength[varid[0]])));

    if len(var_set) == 1: v = [1];
    else: v = [1, len(var_set)];
    if not multivariant: v.extend(range(2, len(var_set)));

    setid = 0;
    for currentset in chain(*map(lambda x: combinations(var_set, x), v)): # For each set combination in powerset of variants in set excluding empty set
        if setid > 0 and setid >= maxh+len(var_set)+1 and len(readentry[setid-1]['_VARSET_']) != len(currentset): break; # Stop if combinations limit reached and finished current k in n choose k
        if setid not in readentry: 
            refentry[setid] = {};
            altentry[setid] = {};
            readentry[setid] = {};
            readentry[setid]['_VARSET_'] = currentset;

        # Construct the variant sequence
        offset = 0;
        altseq = refseq[varid[0]];
        for v in sorted(currentset):
            pos = v[0] - 1 + offset;
            if v[1] == '-': # Account for '-' representation of variants, insertion
                pos += 1;
                ref = '';
                alt = i[2];
            elif v[2] == '-': # Account for '-' representation of variants, deletion
                alt = '';
                ref = v[1];
            else:
                ref = v[1];
                alt = v[2];
            offset += len(alt) - len(ref);
            altseq = altseq[:pos] + alt.encode('utf-8') + altseq[(pos+len(ref)):]; # Explicit byte string for python3, needed for cffi char*
        altseqlength = len(altseq);

        for read in readset: 
            if read.is_unmapped: continue;
            if primaryonly and read.is_secondary: continue;

            if read.is_read1: readid = '{0}/1'.format(read.query_name);
            else: readid = '{0}/2'.format(read.query_name);

            # Convert read sequence into a probability matrix based on base-call error
            callerror = np.array(read.query_qualities) / -10.0; # Already converted to number by pysam
            callerror[callerror == 0] = -0.01; # Ensure there are no zeros, defaulting to a very high base-call error probability
            isbase = np.log10(1 - np.power(10, callerror)); # log10(1-e);
            notbase = callerror - l3; #log10(e/3)

            readlength = len(read.query_sequence);
            readprobmatrix = np.zeros(readlength*5);
            p_readprobmatrix = ffi.cast("double *", readprobmatrix.ctypes.data); # Pointer to read probability matrix
            C.setReadProbMatrix(read.query_sequence.encode('utf-8'), readlength, list(isbase), list(notbase), p_readprobmatrix);

            # Reference genome probability and "elsewhere" probability only needs to be calculated once per readid
            if setid == 0:
                # Calculate the probability for "elsewhere", assuming the read is from somewhere paralogous. Approximate probability distribution by accounting for the bulk with:
                #   perfect match = prod[ (1-e) ]
                #   hamming/edit distance 1 = prod[ (1-e) ] * sum[ (e/3) / (1-e) ]
                # We also account for if reads have different lengths (hard clipped), where longer reads should have a relatively lower probability of originating from some paralogous elsewhere 
                #   lengthfactor = alpha ^ (readlength - expected readlength)
                # P(elsewhere) = (perfect + hamming) / lengthfactor
                elsewhereprobability = C.log10addexp(sum(isbase), sum(isbase) + C.log10sumexp(list(notbase - isbase), len(isbase))) - (lalpha * (readlength - read.infer_query_length())); 
                if readid not in readentry[setid]: readentry[setid][readid] = elsewhereprobability;
                else: readentry[setid][readid] = C.log10addexp(readentry[setid][readid], elsewhereprobability);

                # Calculate the probability given reference genome, once per varid for setid 0
                readprobability = C.calcReadProbability(refseq[samfile.getrname(read.reference_id)], reflength[samfile.getrname(read.reference_id)], read.reference_start, readlength, p_readprobmatrix);
                if readid not in refentry[setid]: refentry[setid][readid] = readprobability;
                else: refentry[setid][readid] = C.log10addexp(refentry[setid][readid], readprobability);

            # Calculate the probability given alternate genome
            readprobability = C.calcReadProbability(altseq, altseqlength, read.reference_start, readlength, p_readprobmatrix);
            if readid not in altentry[setid]: altentry[setid][readid] = readprobability;
            else: altentry[setid][readid] = C.log10addexp(altentry[setid][readid], readprobability);
            if debug: print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(setid, refentry[0][readid], altentry[setid][readid], readentry[0][readid], read, currentset));

            # Multi-mapped alignments
            if not primaryonly and read.has_tag('XA'):
                for j in read.get_tag('XA').split(';'):
                    if len(j) <= 0: break; # Because last element is empty
                    t = j.split(',');
                    if t[0] not in refseq: continue;
                    xa_pos = int(t[1]);
                    if (read.is_reverse == False and xa_pos < 0) or (read.is_reverse == True and xa_pos > 0): # If aligned strand is opposite of that from primary alignment
                        newreadseq = ''.join(complement.get(base,base) for base in reversed(read.query_sequence)).encode('utf-8'); # Reverse complement read sequence
                        newreadprobmatrix = np.zeros(readlength*5);
                        p_readprobmatrix = ffi.cast("double *", newreadprobmatrix.ctypes.data); # Pointer to read probability matrix
                        C.setReadProbMatrix(newreadseq, readlength, list(isbase[::-1]), list(notbase[::-1]), p_readprobmatrix);
                    else: 
                        p_readprobmatrix = ffi.cast("double *", readprobmatrix.ctypes.data); # Pointer to read probability matrix
                    xa_pos = abs(xa_pos);

                    # Probability given reference genome
                    readprobability = C.calcReadProbability(refseq[t[0]], reflength[t[0]], xa_pos, readlength, p_readprobmatrix);
                    if setid == 0: 
                        refentry[setid][readid] = C.log10addexp(refentry[setid][readid], readprobability);
                        # The more multi-mapped, the more likely it is the read is from elsewhere (paralogous), hence it scales (multiplied) with the number of multi-mapped locations
                        readentry[setid][readid] = C.log10addexp(readentry[setid][readid], elsewhereprobability);
                    if debug: print(readprobability, end='\t');
                    # Probability given alternate genome
                    if t[0] == samfile.getrname(read.reference_id): # If secondary alignments are in same chromosome (ie. contains variant thus has modified coordinates), in case it also crosses the variant position, otherwise is the same as probability given reference
                        if xa_pos > varid[1]-1: xa_pos += offset;
                        readprobability = C.calcReadProbability(altseq, altseqlength, xa_pos, readlength, p_readprobmatrix);
                    altentry[setid][readid] = C.log10addexp(altentry[setid][readid], readprobability);
                    if debug: print('{0}\t{1}\t{2}'.format(readprobability, refentry[0][readid], altentry[setid][readid]));
        setid += 1;
    if not readentry: return([]); # Return empty list if no read data

    outlist = [];
    ref = 0.0;
    alt = {};
    het = {};
    refcount = {};
    altcount = {};

    refprior = np.log10(0.5);
    if len(var_set) == 1 or multivariant: # Eeither one variant or multiple variants as a haplotype, homozygous & non-homozygous
        altprior = np.log10(0.5 * (1-hetbias));
        hetprior = np.log10(0.5 * hetbias);
    else:  # Variant set: divided evenly among the variant hypotheses, homozygous & non-homozygous
        altprior = np.log10(0.5 * (1-hetbias) / len(readentry));
        hetprior = np.log10(0.5 * hetbias / len(readentry));

    for setid in readentry:
        currentset = readentry[setid]['_VARSET_'];
        if currentset not in alt: 
            alt[currentset] = 0.0;
            het[currentset] = 0.0;
            refcount[currentset] = 0;
            altcount[currentset] = 0;

        for readid in refentry[0]:
            prgu = refentry[0][readid]; # log of P(r|Gu), where Gu is the unchanged from the reference genome
            prgv = altentry[setid][readid]; # log of P(r|Gv), where Gv is the variant genome
            pelsewhere = readentry[0][readid];

            # Mixture model: probability that the read is from elsewhere, outside paralogous source
            prgu = C.log10addexp(lomega - l1_omega + pelsewhere, prgu);
            prgv = C.log10addexp(lomega - l1_omega + pelsewhere, prgv);

            # Mixture model: heterozygosity or heterogeneity as explicit allele frequency mu such that P(r|GuGv) = (mu)(P(r|Gv)) + (1-mu)(P(r|Gu))
            phet = C.log10addexp(l50 + prgv, l50 + prgu);
            phet10 = C.log10addexp(l10 + prgv, l90 + prgu);
            phet90 = C.log10addexp(l90 + prgv, l10 + prgu);
            phet = max([phet, phet10, phet90]); # Use the best allele frequency probability

            # Read count is only incremented when the difference in probability is not ambiguous, ~log(2) difference
            if prgv > prgu and prgv - prgu > 0.3: altcount[currentset] += 1;
            elif prgu > prgv and prgu - prgv > 0.3: refcount[currentset] += 1;
            
            if setid == 0: ref += prgu + refprior; # Only one reference hypothesis
            alt[currentset] += prgv + altprior;
            het[currentset] += phet + hetprior;
            if debug: print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}'.format(prgu, phet, prgv, pelsewhere, varid[0], currentset, readid, altcount[currentset])); # log likelihoods
        if debug: print('-=-\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(ref, het[currentset], alt[currentset], varid[0], currentset, altcount[currentset])); # log likelihoods

    total = [ref] + list(alt.values()) + list(het.values());
    total = C.log10sumexp(total, len(total));
    for v in var_set:
        marginal_alt = [];
        not_alt = [ref];
        marginal_count = [];
        for currentset in alt:
            if v in currentset: 
                marginal_alt.extend([ alt[currentset], het[currentset] ]);
                marginal_count.extend([ altcount[currentset] ]);
            else: 
                not_alt.extend([ alt[currentset], het[currentset] ]);
        outstr = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t'.format(varid[0], v[0], v[1], v[2], max(refcount.values())+max(altcount.values()), max(marginal_count));
        # Probability and odds in log10
        m = C.log10sumexp(marginal_alt, len(marginal_alt));
        outstr += '{0}\t{1}\t'.format(m - total, m - C.log10sumexp(not_alt, len(not_alt))); 
        # Print the variant set entries if exists    
        if len(var_set) > 1: outstr += '{0}'.format(var_set);
        else: outstr += '[]';
        outlist.append(outstr.strip());
    return(outlist);

# Global Variables
debug = False;
primaryonly = False;
numprocesses = 1;
hetbias = 0.5;
distancethreshold = 10;
maxh = 1024;
multivariant = False;
refseq = {};
reflength = {};
def main():
    parser = argparse.ArgumentParser(description='Explicit alternative genome likelihood evaluator.');
    parser.add_argument('-v', help='variants VCF file');
    parser.add_argument('-a', nargs='+', help='alignment data bam files');
    parser.add_argument('-r', help='reference sequence fasta file');
    parser.add_argument('-o', type=str, default='', help='output file (default: stdout)');
    parser.add_argument('-hetbias', type=float, default=0.5, help='prior probability bias towards non-homozygous mutations (value between [0,1], default: 0.5 unbiased)');
    parser.add_argument('-n', type=int, default=10, help='consider nearby variants within n bases in the set of hypotheses (off: 0, default: 10)');
    parser.add_argument('-maxh', type=int, default=1024, help='the maximum number of combinations in the set of hypotheses, instead of all 2^n (default: 2^10 = 1024)');
    parser.add_argument('-mvh', action='store_true', help='consider nearby variants as *one* multi-variant hypothesis');
    parser.add_argument('-p', action='store_true', help='consider only primary alignments');
    parser.add_argument('-t', type=int, default=1, help='number of processes to use (default: 1)');
    parser.add_argument('-debug', action='store_true', help='debug mode, printing information on every read for every variant');
    if len(sys.argv) == 1:
        parser.print_help();
        sys.exit(1);
    args = parser.parse_args();

    global debug, primaryonly, numprocesses, hetbias, distancethreshold, maxh, multivariant;
    debug = args.debug;
    primaryonly = args.p;
    numprocesses = args.t;
    hetbias = args.hetbias;
    distancethreshold = args.n;
    maxh = args.maxh;
    multivariant = args.mvh;

    refvar_list = readVCF(args.v);
    if distancethreshold > 0: refvar_list = groupNearbyVariants(refvar_list);

    global refseq, reflength;
    (refseq, reflength) = readFasta(args.r);
    readPYSAM(args.a, refvar_list, args.o); # reads crossing variant in genome

if __name__ == '__main__':
    try:
        main();
    except KeyboardInterrupt: 
        try:
            sys.exit(0);
        except SystemExit:
            os._exit(0);
