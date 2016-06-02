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
from scipy.misc import logsumexp;
from itertools import chain, combinations;
from multiprocessing import Pool, Process, Manager;
from datetime import datetime;
from signal import signal, SIGPIPE, SIG_DFL;
signal(SIGPIPE,SIG_DFL);

# Constants
omega = 1E-4; # Prior probability of read originating from an outside paralogous source
alpha = 1.3; # Factor to account for longer read lengths lowering the probability a sequence matching an outside paralogous source

# Precalculated log values 
logln = np.log10(np.e);
log3 = np.log10(3);
ln50 = np.log(0.5);
ln10 = np.log(0.1);
ln90 = np.log(0.9);
lnomega = np.log(omega);
ln1_omega = np.log(1-omega);
lnalpha = np.log(alpha);

complement = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a' };

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
            ref = var[3].split(',');
            alt = var[4].split(',');
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
                seq[seqid].append(line);
    fh.close;
    for seqid in seq: 
        seq[seqid] = ''.join(seq[seqid]).encode('utf-8'); # Explicit unicode for python3, needed for cffi char*
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

        # Testing speed with Process, seems slower than pool
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
            if varid[n][0] in refseq: args.append((fn, varid[n], var_list[varid[n]]));
        try:
            pool = Pool(processes=numprocesses);
            poolresults = pool.map_async(evaluateVariant, args);
            for p in poolresults.get(): results.extend(p);
        finally:
            pool.close();
            pool.join();

    if len(outfile) > 0: fh = open(outfile, 'w');
    else: fh = sys.stdout;
    print('#SEQ\tPOS\tREF\tALT\tReads\tAltReads\tProb(log10)\tOdds(log10)\tVarSet', file=fh);
    for i in naturalSort(results): print(i, file=fh);
    fh.close();
    print("Done:\t{0}\t{1}".format(fn, datetime.now()), file=sys.stderr);

def evaluateVariant(args):
    (fn, varid, var_set) = args;
    refentry = {};
    altentry = {};
    readentry = {};
    readentry[varid] = {};
    refentry[varid] = {};
    altentry[varid] = {};
    varstart = min([a[0] for a in var_set]);
    varend = max([a[0] for a in var_set]);

    if len(var_set) == 1: i = [1];
    else: i = [1, len(var_set)];
    if not multivariant: i.extend(range(2, len(var_set)));
    hypotheses = chain(*map(lambda x: combinations(var_set, x), i)); # Powerset of variants in set excluding empty set

    setid = 0;
    for currentset in hypotheses:
        if setid >= maxh+len(var_set)+1 and setid > 0 and len(readentry[varid][setid-1]['_VARSET_']) != len(currentset): break; # Stop if max hypotheses limit reached and finished current n choose k
        if setid not in readentry[varid]: 
            refentry[varid][setid] = {};
            altentry[varid][setid] = {};
            readentry[varid][setid] = {};
            readentry[varid][setid]['_VARSET_'] = currentset;

        # Construct the variant sequence
        offset = 0;
        altseq = refseq[varid[0]];
        for i in sorted(currentset):
            pos = i[0] - 1 + offset;
            if i[1] == '-': # Account for '-' representation of variants, insertion
                pos += 1;
                ref = '';
                alt = i[2];
            elif i[2] == '-': # Account for '-' representation of variants, deletion
                alt = '';
                ref = i[1];
            else:
                ref = i[1];
                alt = i[2];
            offset += len(alt) - len(ref);
            altseq = altseq[:pos] + alt.encode('utf-8') + altseq[(pos+len(ref)):]; # Explicit unicode for python3, needed for cffi char*
        altseqlength = len(altseq);

        # Fetch reads that cross the variant location
        samfile = pysam.AlignmentFile(fn, "rb");
        readset = samfile.fetch(reference=varid[0], start=max(0,varstart-1), end=min(varend+1,reflength[varid[0]]));
        for read in readset: 
            if read.is_unmapped: continue;
            if primaryonly and read.is_secondary: continue;

            readid = read.query_name;
            if read.is_read1: readid += '/1';
            elif read.is_read2: readid += '/2';

            readlength = len(read.query_sequence.encode('utf-8'));
            callerror = np.array(read.query_qualities) / float(-10); # Already converted to ord by pysam
            callerror[callerror == 0] = -0.01; # Ensure there are no zeros, defaulting to a very high base-call error probability
            isbase = np.log10(1 - np.power(10, callerror)); # log10(1-e);
            notbase = callerror - log3; #log10(e/3)
            # Convert read sequence into a probability matrix based on base-call error
            readprobmatrix = ffi.new("double[]", readlength*5);
            C.setReadProbMatrix(read.query_sequence.encode('utf-8'), readlength, list(isbase), list(notbase), readprobmatrix);

            # Calculate the probability for "elsewhere", assuming the read is from somewhere paralogous
            #   perfect & edit distance 1: to approximate probability distribution of a read that describes a paralog elsewhere, this account for the bulk of the probability distribution
            #   we account for if reads have different lengths, where longer reads should have a lower probability of originating from some paralogous elsewhere 
            elsewhereprobability = np.logaddexp(sum(isbase) / logln, (sum(isbase) / logln) + logsumexp((notbase - isbase) / logln)) - (lnalpha * (readlength - read.infer_query_length())); 

            if readid not in readentry[varid][setid]: readentry[varid][setid][readid] = elsewhereprobability;
            else: readentry[varid][setid][readid] = np.logaddexp(readentry[varid][setid][readid], elsewhereprobability);

            # Calculate the probability given reference genome
            readprobability = calcReadProbability(refseq[samfile.getrname(read.reference_id)], reflength[samfile.getrname(read.reference_id)], read.reference_start, readlength, readprobmatrix);
            if readid not in refentry[varid][setid]: refentry[varid][setid][readid] = readprobability;
            else: refentry[varid][setid][readid] = np.logaddexp(refentry[varid][setid][readid], readprobability);

            # Calculate the probability given alternate genome
            readprobability = calcReadProbability(altseq, altseqlength, read.reference_start, readlength, readprobmatrix);
            if readid not in altentry[varid][setid]: altentry[varid][setid][readid] = readprobability;
            else: altentry[varid][setid][readid] = np.logaddexp(altentry[varid][setid][readid], readprobability);
            if debug: print('{0}\t{1}\t{2}\t{3}\t{4}'.format(refentry[varid][setid][readid], altentry[varid][setid][readid], readentry[varid][setid][readid], read, currentset));

            # Multi-mapped alignments
            if not primaryonly and read.has_tag('XA'):
                for j in read.get_tag('XA').split(';'):
                    if len(j) <= 0: break; # Because last element is empty
                    t = j.split(',');
                    if t[0] not in refseq: continue;
                    xa_pos = int(t[1]);
                    if (read.is_reverse == False and xa_pos < 0) or (read.is_reverse == True and xa_pos > 0): # If strand is opposite of that from primary alignment
                        newreadseq = ''.join(complement.get(base,base) for base in reversed(read.query_sequence)).encode('utf-8');
                        newreadprobmatrix = ffi.new("double[]", readlength*5);
                        C.setReadProbMatrix(newreadseq, readlength, list(isbase[::-1]), list(notbase[::-1]), newreadprobmatrix);
                    else: 
                        newreadprobmatrix = readprobmatrix;
                    xa_pos = abs(xa_pos);

                    # The more multi-mapped, the more likely it is the read is from elsewhere (paralogous), hence it scales (multiplied) with the number of multi-mapped locations
                    readentry[varid][setid][readid] = np.logaddexp(readentry[varid][setid][readid], elsewhereprobability);
                    # Probability given reference genome
                    readprobability = calcReadProbability(refseq[t[0]], reflength[t[0]], xa_pos, readlength, newreadprobmatrix);
                    refentry[varid][setid][readid] = np.logaddexp(refentry[varid][setid][readid], readprobability);
                    if debug: print(readprobability, end='\t');
                    # Probability given alternate genome
                    if t[0] == samfile.getrname(read.reference_id): # If secondary alignments are in same chromosome (ie. contains variant thus has modified coordinates), in case it also crosses the variant position, otherwise is the same as probability given reference
                        if xa_pos > varid[1]-1: xa_pos += offset;
                        readprobability = calcReadProbability(altseq, altseqlength, xa_pos, readlength, newreadprobmatrix);
                    altentry[varid][setid][readid] = np.logaddexp(altentry[varid][setid][readid], readprobability);
                    if debug: print('{0}\t{1}\t{2}'.format(readprobability, refentry[varid][setid][readid], altentry[varid][setid][readid]));
        setid += 1;

    outlist = [];
    if readentry[varid]: # Compile probabilities from read data if exists
        ref = float(0);
        alt = {};
        het = {};
        refcount = {};
        altcount = {};
        currentset = ();

        for setid in readentry[varid]:
            currentset = readentry[varid][setid]['_VARSET_'];
            if currentset not in alt: 
                alt[currentset] = float(0);
                het[currentset] = float(0);
                refcount[currentset] = 0;
                altcount[currentset] = 0;

            refprior = 0.5;
            if multivariant or len(var_set) == 1: altprior = np.log((1-refprior) / float(2)); # one hypothesis, either one variant or multiple variants as a haplotype, homozygous & non-homozygous
            else: altprior = np.log((1-refprior) / float(len(var_set)*2)); # remainder divided evenly among the variant hypotheses, homozygous & non-homozygous
            refprior = np.log(refprior);

            for readid in refentry[varid][setid]:
                prgu = refentry[varid][setid][readid]; # ln of P(r|Gu), where Gu is the unchanged from the reference genome
                prgv = altentry[varid][setid][readid]; # ln of P(r|Gv), where Gv is the variant genome
                pelsewhere = readentry[varid][setid][readid];

                # Mixture model: probability that the read is from elsewhere, outside paralogous source
                prgu = np.logaddexp(lnomega - ln1_omega + pelsewhere, prgu);
                prgv = np.logaddexp(lnomega - ln1_omega + pelsewhere, prgv);

                # Mixture model: heterozygosity or heterogeneity as explicit allele frequency mu such that P(r|GuGv) = (mu)(P(r|Gv)) + (1-mu)(P(r|Gu))
                phet = np.logaddexp(ln50 + prgv, ln50 + prgu);
                phet10 = np.logaddexp(ln10 + prgv, ln90 + prgu);
                phet90 = np.logaddexp(ln90 + prgv, ln10 + prgu);
                phet = max([phet, phet10, phet90]); # Use the best allele frequency probability

                # Read count is only incremented when the difference in probability is not ambiguous, ln(2) difference
                if prgv > prgu and prgv - prgu > 0.69: altcount[currentset] += 1;
                elif prgu > prgv and prgu - prgv > 0.69: refcount[currentset] += 1;
                
                if setid == 0: ref += prgu + refprior; # Only one reference hypothesis
                alt[currentset] += prgv + altprior;
                het[currentset] += phet + altprior;
                if debug: print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}'.format(prgu, prgv, pelsewhere, varid[0], currentset, readid, altcount[currentset])); # ln likelihoods
            if debug: print('-=-\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(ref, het[currentset], alt[currentset], varid[0], currentset, altcount[currentset])); # ln likelihoods

        total = logsumexp( [ref] + list(alt.values()) + list(het.values()) );
        for i in var_set:
            marginal_alt = [];
            not_alt = [ref];
            marginal_count = [];
            for v in alt:
                if i in v: 
                    marginal_alt.extend([ alt[v], het[v] ]);
                    marginal_count.extend([ altcount[v] ]);
                else: 
                    not_alt.extend([ alt[v], het[v] ]);
            outstr = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t'.format(varid[0], i[0], i[1], i[2], max(refcount.values())+max(altcount.values()), max(marginal_count));
            # Probability and odds in log10
            outstr += '{0}\t{1}\t'.format((logsumexp(marginal_alt) - total) / np.log(10), (logsumexp(marginal_alt) - logsumexp(not_alt)) / np.log(10)); 
            # Print the variant set entries if exists    
            if len(var_set) > 1: outstr += '{0}'.format(var_set[0]);
            else: outstr += '[]';
            outlist.append(outstr.strip());
    return(outlist);

ffi = FFI();
ffi.cdef ("""
void initAlphaMap(void);
void setReadProbMatrix(char*, int, double*, double*, double*);
double getReadProb(char*, int, int, int, double*, double);
void getReadProbList(char*, int, int, int, double*, double*);
""")
C = ffi.verify ("""
static int alphabetval[26];
static void initAlphaMap(void) {
    memset(alphabetval, 4, sizeof(alphabetval)); // Zero out array
    alphabetval['A'-'A'] = 0;
    alphabetval['T'-'A'] = 1;
    alphabetval['G'-'A'] = 2;
    alphabetval['C'-'A'] = 3;
    alphabetval['N'-'A'] = 4;
}
void setReadProbMatrix(char* seq, int readlength, double *isbase, double *notbase, double *matrix) {
    int i, b; //array[width * row + col] = value
    memset(matrix, 0, readlength*5*sizeof(double)); // Zero out array
    for ( b = 0; b < readlength; b++ ) {
        for ( i = 0; i < 5; i++ ) matrix[5*b+i] = notbase[b];
        matrix[5*b+alphabetval[toupper(seq[b])-'A']] = isbase[b];
    }
}
double getReadProb(char *seq, int seqlength, int refpos, int readlength, double *matrix, double baseline) {
    int b; //array[width * row + col] = value
    double probability = 0;
    //printf("%d\\t", refpos);
    for ( b = refpos;  b < refpos+readlength; b++ ) {
        if ( b < 0) continue; // Skip if position is before the start of reference seq
        if ( b >= seqlength ) break; // Stop if it reaches the end of reference seq
        //printf("%c", seq[b]);
        probability += matrix[5*(b-refpos)+alphabetval[toupper(seq[b])-'A']]; 
        if ( probability < baseline - 10 ) break; // Stop sum if less than 1% contribution to baseline (best, highest) probability mass
    }
    //printf("\\t%f\\n", probability);
    return (probability);
}
void getReadProbList(char *seq, int seqlength, int refpos, int readlength, double *matrix, double *probabilityarray) {
    memset(probabilityarray, 0, readlength*2*sizeof(double)); // Zero out array
    int i;
    int n = 0;
    int slidelength = readlength;
    double baseline = getReadProb(seq, seqlength, refpos, readlength, matrix, -1000); // First probability at refpos, likely the highest, to be used as first best
    for ( i = refpos-slidelength; i <= refpos+slidelength; i++ ) {
        if ( i + readlength < 0 ) continue;
        if ( i >= seqlength ) break;
        probabilityarray[n] = getReadProb(seq, seqlength, i, readlength, matrix, baseline);
        if ( probabilityarray[n] > baseline ) baseline = probabilityarray[n]; // Update the best probability so far
        n++;
    }
}
""")
C.initAlphaMap(); # Initialize alphabet to int mapping table
def calcReadProbability(refseq, reflength, refpos, readlength, readprobmatrix):
    probabilityarray = np.zeros(readlength*2);
    p_probabilityarray = ffi.cast("double *", probabilityarray.ctypes.data); # Pointer to probability array
    C.getReadProbList(refseq, reflength, refpos, readlength, readprobmatrix, p_probabilityarray);
    return(logsumexp(probabilityarray[np.nonzero(probabilityarray)] / logln)); # Return sum of probabilities as ln exponential (converted from log10)

# Global Variables
debug = False;
primaryonly = False;
numprocesses = 1;
distancethreshold = 10;
maxh = 1024;
multivariant = False;
refseq = {};
reflength = {};
def main():
    parser = argparse.ArgumentParser(description='Evaluate the significance of alternative genome using an explicit probabilistic model');
    parser.add_argument('-v', help='variants VCF file');
    parser.add_argument('-a', nargs='+', help='alignment data bam files');
    parser.add_argument('-r', help='reference sequence fasta file');
    parser.add_argument('-o', type=str, default='', help='output file (default: stdout)');
    parser.add_argument('-n', type=int, default=10, help='consider nearby variants within n bases in the set of hypotheses (off: 0, default: 10)');
    parser.add_argument('-maxh', type=int, default=1024, help='the maximum number of hypotheses, instead of all 2^n (default: 2^10 = 1024)');
    parser.add_argument('-mvh', action='store_true', help='consider nearby variants as *one* multi-variant hypothesis');
    parser.add_argument('-p', action='store_true', help='consider only primary alignments');
    parser.add_argument('-t', type=int, default=1, help='number of processes to use (default: 1)');
    parser.add_argument('-debug', action='store_true', help='debug mode, printing information on every read for every variant');
    if len(sys.argv) == 1:
        parser.print_help();
        sys.exit(1);
    args = parser.parse_args();

    global debug, primaryonly, numprocesses, distancethreshold, maxh, multivariant;
    debug = args.debug;
    primaryonly = args.p;
    numprocesses = args.t;
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
