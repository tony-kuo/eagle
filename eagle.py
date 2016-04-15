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
from multiprocessing import Process, Queue, Pool;
from datetime import datetime;
from signal import signal, SIGPIPE, SIG_DFL;
signal(SIGPIPE,SIG_DFL);

# Constants
omega = 1E-6; # Prior probability read is from some "elsewhere" that is paralogous
lg = np.log10(omega);
e3 = np.log10(3);
l50 = np.log(0.5);
l10 = np.log(0.1);
l90 = np.log(0.9);
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
            if len(entry[varid[i]]) <= 1:
                if multivariant: del entry[varid[i]]; # Remove singletons if we only consider the multi-variant hypothesis
                continue;
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
    misc = {};
    varid = sorted(entry.keys());
    for i in range(0, len(varid)):
        if len(entry[varid[i]]) <= 1: continue; # Ignore singletons
        x = max([entry[varid[i]][j+1][0] - entry[varid[i]][j][0] for j in range(0, len(entry[varid[i]])-1)]); # If all variants in set are adjacent to each other
        y = max([max(len(j[1]), len(j[2])) for j in entry[varid[i]]]); # If all variants are SNPs
        if x == 1 and y == 1:
            ref = ''.join([j[1] for j in entry[varid[i]]]);
            alt = ''.join([j[2] for j in entry[varid[i]]]);
            misc[(varid[i][0], varid[i][1], ref, alt)] = [(varid[i][1], ref, alt)]; # Merge and delete extraneous entries
            del entry[varid[i]];
    if misc: entry.update(misc);
    # Organize the sets so that every variant of interest is the first in set while the rest contain nearby variants
    #if not multivariant:
        #misc = {};
        #varid = sorted(entry.keys());
        #for i in range(0, len(varid)):
            #if len(entry[varid[i]]) <= 1: continue; # Ignore singletons
            #for j in entry[varid[i]]:
                #newid = (varid[i][0], j[0], j[1], j[2]);
                #misc[newid] = [j];
                #for k in entry[varid[i]]:
                    #if j != k and abs(j[0] - k[0]) <= distancethreshold: misc[newid].append(k);
            #del entry[varid[i]];
        #if misc: entry.update(misc);

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
                if len(seqid) < 3: seqid = line[1:]; # >1 or >chr1
                seq[seqid] = [];
            else:
                seq[seqid].append(line);
    fh.close;
    for seqid in seq: 
        seq[seqid] = ''.join(seq[seqid]).encode('utf-8'); # Explicit unicode for python3, needed for cffi char*
        seqlength[seqid] = len(seq[seqid]);
    return(seq, seqlength);

def readPYSAM(files, var_list, outfile):
    entry = [];
    for fn in files: 
        print("Start:\t{0}\t{1}".format(fn, datetime.now()), file=sys.stderr);
        varid = list(var_list.keys());
        numvariants = len(varid);
        #entry = [evaluateVariant(fn, varid[n], var_list[varid[n]]) for n in range(0, numvariants)];

        try:
            pool = Pool(processes=numprocesses);
            results = [pool.apply_async(evaluateVariant, args=(fn, varid[n], var_list[varid[n]])) for n in range(0, numvariants)];
        finally:
            pool.close();
            pool.join();
        for p in results:
            for j in p.get():
                if j: entry.append(j); # Only keep non-empty results

        #n = 0;
        #outqueue = Queue();
        #while n < numvariants:
            #procs = [];
            #for i in range(0, numprocesses):
                #if n >= numvariants: break;
                #p = Process(target=evaluateVariant, args=(outqueue, fn, varid[n], var_list[varid[n]]));
                #procs.append(p);
                #p.start();
                #print("Read {0}:\t{1} of {2} variants\t{3}\t{4}\t{5}".format(fn, n, numvariants, varid[n], i, datetime.now()), file=sys.stderr);
                #n += 1;
            #for i in procs: # Get results from each process from queue
                #for j in outqueue.get(): 
                    #if j: entry.append(j); # Only keep non-empty results
            #for i in procs: i.join();
    if len(outfile) > 0: fh = open(outfile, 'w');
    else: fh = sys.stdout;
    print('#SEQ\tPOS\tREF\tALT\tReads\tAltReads\tProb(log10)\tOdds(log10)\tVarSet', file=fh);
    for i in naturalSort(entry): print(i, file=fh);
    fh.close();
    print("Done:\t{0}\t{1}".format(fn, datetime.now()), file=sys.stderr);

def evaluateVariant(fn, varid, var_set):
    refentry = {};
    altentry = {};
    readentry = {};
    readentry[varid] = {};
    refentry[varid] = {};
    altentry[varid] = {};
    varstart = min([a[0] for a in var_set]);
    varend = max([a[0] for a in var_set]);

    if multivariant: hypotheses = list(combinations(var_set, len(var_set)));
    elif len(var_set) > maxk: 
        hypotheses = list(combinations(var_set, 1)); # Solo variant hypotheses
        hypotheses.extend(list(combinations(var_set, len(var_set)))); # All variant co-occurs hypothesis
        for i in range(2, int(np.sqrt(len(var_set)))): hypotheses.extend(list(combinations(var_set, i))); # n choose k variant combination hypotheses, up to n choose sqrt(n)
    else: hypotheses = chain(*map(lambda x: combinations(var_set, x), range(1, len(var_set)+1))); # powerset of variants in set excluding empty set

    setid = 0;
    for currentset in hypotheses:
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
                ref = chr(refseq[varid[0]][ i[0]-1 ]);
                alt = ref+i[2];
            elif i[2] == '-': # Account for '-' representation of variants, deletion
                alt = chr(refseq[varid[0]][ i[0]-2 ]);
                ref = alt + i[1];
            else:
                ref = i[1];
                alt = i[2];
            offset += len(alt) - len(ref);
            altseq = altseq[:pos] + alt.encode('utf-8') + altseq[(pos+len(ref)):]; # Explicit unicode for python3, needed for cffi char*
        altseqlength = len(altseq);
        # Fetch reads that cross the variant location
        samfile = pysam.AlignmentFile(fn, "rb");
        readset = samfile.fetch(reference=varid[0], start=varstart-1, end=varend+1);
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
            notbase = callerror - e3; #log10(e/3)
            # Convert read sequence into a probability matrix based on base-call error
            readprobmatrix = ffi.new("double[]", readlength*5);
            C.setReadProbMatrix(read.query_sequence.encode('utf-8'), readlength, list(isbase), list(notbase), readprobmatrix);

            # Calculate the probability for "elsewhere", assuming the read is correct but is from somewhere paralogous
            elsewhereprobability = sum([logsumexp(np.array([isbase[a]*2, (callerror[a]*2)-e3]) / np.log10(np.e)) for a in range(0,len(callerror))]); # (1-e)^2 + e^2/3
            readentry[varid][setid][readid] = elsewhereprobability;

            # Calculate the probability given reference genome
            readprobability = calcReadProbability(refseq[samfile.getrname(read.reference_id)], reflength[samfile.getrname(read.reference_id)], read.reference_start, readlength, readprobmatrix, elsewhereprobability);
            if readid not in refentry[varid][setid]: refentry[varid][setid][readid] = readprobability;
            else: refentry[varid][setid][readid] = logsumexp([refentry[varid][setid][readid], readprobability]);

            # Calculate the probability given alternate genome
            readprobability = calcReadProbability(altseq, altseqlength, read.reference_start, readlength, readprobmatrix, elsewhereprobability);
            if readid not in altentry[varid][setid]: altentry[varid][setid][readid] = readprobability;
            else: altentry[varid][setid][readid] = logsumexp([altentry[varid][setid][readid], readprobability]);
            if debug: print('{0}\t{1}\t{2}\t{3}\t{4}'.format(refentry[varid][setid][readid], altentry[varid][setid][readid], readentry[varid][setid][readid], read, currentset));

            # Multi-mapped alignments
            if not primaryonly and read.has_tag('XA'):
                for j in read.get_tag('XA').split(';'):
                    if len(j) <= 0: break; # Because last element is empty
                    t = j.split(',');
                    xa_pos = int(t[1]);
                    if (read.is_reverse == False and xa_pos < 0) or (read.is_reverse == True and xa_pos > 0): # If strand is opposite of that from primary alignment
                        newreadseq = ''.join(complement.get(base,base) for base in reversed(read.query_sequence)).encode('utf-8');
                        newreadprobmatrix = ffi.new("double[]", readlength*5);
                        C.setReadProbMatrix(newreadseq, readlength, list(isbase[::-1]), list(notbase[::-1]), newreadprobmatrix);
                    else: newreadprobmatrix = readprobmatrix;
                    xa_pos = abs(xa_pos);

                    # The more multi-mapped, the more likely it the read is from the outside
                    readentry[varid][setid][readid] += elsewhereprobability;
                    # Probability given reference genome
                    readprobability = calcReadProbability(refseq[t[0]], reflength[t[0]], xa_pos, readlength, newreadprobmatrix, elsewhereprobability);
                    refentry[varid][setid][readid] = logsumexp([refentry[varid][setid][readid], readprobability]);
                    if debug: print(readprobability, end='\t');
                    # Probability given alternate genome
                    if t[0] == samfile.getrname(read.reference_id): # Secondary alignments are in same chromosome (ie. contains variant thus has modified coordinates), just in case multi-mapped position also crosses the variant position
                        if xa_pos > varid[1]-1: xa_pos += offset;
                        readprobability = calcReadProbability(altseq, altseqlength, xa_pos, readlength, newreadprobmatrix, elsewhereprobability);
                    altentry[varid][setid][readid] = logsumexp([altentry[varid][setid][readid], readprobability]);
                    if debug: print('{0}\t{1}\t{2}'.format(readprobability, refentry[varid][setid][readid], altentry[varid][setid][readid]));
        setid += 1;

    outlist = [];
    if readentry[varid]: # Compile probabilities from read data if exists
        ref = {};
        alt = {};
        het = {};
        het10 = {};
        het90 = {};
        elsewhere = {};
        refcount = {};
        altcount = {};
        currentset = ();

        refprior = 0.5;
        elsewhereprior = omega;
        if multivariant: altprior = np.log(1-refprior-elsewhereprior); # multivariant is one hypothesis for all variants as a haplotype
        else: altprior = np.log((1-refprior-elsewhereprior) / float(len(var_set))); # remainder divided evenly among the variant hypotheses
        refprior = np.log(refprior);
        elsewhereprior = np.log(elsewhereprior);
        for setid in readentry[varid]:
            currentset = readentry[varid][setid]['_VARSET_'];
            if currentset not in ref: 
                ref[currentset] = float(0);
                alt[currentset] = float(0);
                het[currentset] = float(0);
                het10[currentset] = float(0);
                het90[currentset] = float(0);
                elsewhere[currentset] = float(0);
                refcount[currentset] = 0;
                altcount[currentset] = 0;

            for readid in refentry[varid][setid]:
                prgx = refentry[varid][setid][readid] + refprior; # ln of P(r|Gx), where Gx is original genome variant region
                prgv = altentry[varid][setid][readid] + altprior; # ln of P(r|Gv), where Gv is mutated genome variant region
                pelsewhere = readentry[varid][setid][readid] + elsewhereprior;
                # Mixture model: ln of (mu)(P(r|Gv)) + (1-mu)(P(r|Gx))
                phet = logsumexp([l50 + prgv, l50 + prgx]) + altprior;
                phet10 = logsumexp([l10 + prgv, l90 + prgx]) + altprior;
                phet90 = logsumexp([l90 + prgv, l10 + prgx]) + altprior;

                # Read count is only incremented when the difference in probability is not ambiguous
                max_prgv = max([prgv, phet, phet10, phet90]);
                if max_prgv > prgx and max_prgv > pelsewhere and max_prgv-prgx > 0.69: altcount[currentset] += 1; # about ln(2) difference
                elif prgx > max_prgv and prgx > pelsewhere and prgx-max_prgv > 0.69: refcount[currentset] += 1;
                
                ref[currentset] += prgx;
                alt[currentset] += prgv;
                het[currentset] += phet;
                het10[currentset] += phet10;
                het90[currentset] += phet90;
                elsewhere[currentset] += pelsewhere;
                if debug: print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}'.format(prgx, phet, prgv, pelsewhere, varid[0], currentset, readid, altcount[currentset])); # ln likelihoods
            if debug: print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}'.format(ref[currentset], het10[currentset], het[currentset], het90[currentset], alt[currentset], elsewhere[currentset], varid[0], currentset, altcount[currentset])); # ln likelihoods

        total = logsumexp(list(ref.values()) + list(alt.values()) + list(het.values()) + list(het10.values()) + list(het90.values()) + list(elsewhere.values()));
        not_alt = list(ref.values()) + list(elsewhere.values());
        for i in var_set:
            #i = var_set[0]; # Variant of interest is the first in set
            marginal_alt = [];
            for v in alt:
                if i in v: marginal_alt.extend([alt[v], het[v], het10[v], het90[v]]);
                else: not_alt.extend([alt[v], het[v], het10[v], het90[v]]);
            if multivariant: outstr = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t'.format(varid[0], i[0], i[1], i[2], altcount[currentset]+refcount[currentset], altcount[currentset]);
            else: outstr = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t'.format(varid[0], i[0], i[1], i[2], max(altcount.values())+max(refcount.values()), max(altcount.values()));
            # Probability and odds in log10
            outstr += '{0}\t{1}\t'.format((logsumexp(marginal_alt) - total) / np.log(10), (logsumexp(marginal_alt) - logsumexp(not_alt)) / np.log(10)); 
            # Print the variant set entries if exists    
            if len(var_set) > 1: outstr += '{0}'.format(var_set);
            else: outstr += '[]';
            outlist.append(outstr.strip());
    return(outlist);

ffi = FFI();
ffi.cdef ("""
void initAlphaMap(void);
void setReadProbMatrix(char*, int, double*, double*, double*);
double getReadProb(char*, int, int, int, double*);
void getReadProbList(char*, int, int, int, double*, double*);
""")
C = ffi.verify ("""
static int alphabetval[26];
static void initAlphaMap(void) {
    memset(alphabetval, 0, sizeof(alphabetval)); // Zero out array
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
double getReadProb(char *seq, int seqlength, int refpos, int readlength, double *matrix) {
    int b; //array[width * row + col] = value
    double probability = 0;
    //printf("%d\\t", refpos);
    for ( b = refpos;  b < refpos+readlength; b++ ) {
        if ( b < 0) continue; // Skip if position is before the start of reference seq
        if ( b >= seqlength ) break; // Stop if it reaches the end of reference seq
        //printf("%c", seq[b]);
        probability += matrix[5*(b-refpos)+alphabetval[toupper(seq[b])-'A']]; 
    }
    //printf("\\t%f\\n", probability);
    return (probability);
}
void getReadProbList(char *seq, int seqlength, int refpos, int readlength, double *matrix, double *probabilityarray) {
    memset(probabilityarray, 0, readlength*2*sizeof(double)); // Zero out array
    int i;
    int n = 0;
    int slidelength = readlength;
    for ( i = refpos-slidelength; i <= refpos+slidelength; i++ ) {
        if ( i + readlength < 0 ) continue;
        if ( i >= seqlength ) break;
        probabilityarray[n] = getReadProb(seq, seqlength, i, readlength, matrix);
        n++;
    }
}
""")
C.initAlphaMap(); # Initialize alphabet to int mapping table
def calcReadProbability(refseq, reflength, refpos, readlength, readprobmatrix, elsewhereprobability):
    probabilityarray = ffi.new("double[]", readlength*2);
    C.getReadProbList(refseq, reflength, refpos, readlength, readprobmatrix, probabilityarray);
    overall_probability = np.frombuffer(ffi.buffer(probabilityarray));
    overall_probability = logsumexp(overall_probability[np.nonzero(overall_probability)] / np.log10(np.e)); # sum of probabilities as ln exponential (converted from log10)
    #if overall_probability < elsewhereprobability + lg: overall_probability = elsewhereprobability + lg; # floor probability as a fraction of probability from "elsewhere"
    return(overall_probability);

# Global Variables
debug = False;
primaryonly = False;
numprocesses = 1;
distancethreshold = 30;
maxk = 12;
multivariant = False;
refseq = {};
reflength = {};
def main():
    parser = argparse.ArgumentParser(description='Evaluate the significance of alternative genome using an explicit probabilistic model');
    parser.add_argument('-v', help='variants list from getVariantGenome.pl');
    parser.add_argument('-a', nargs='+', help='alignment data bam files');
    parser.add_argument('-r', help='reference sequence fasta file');
    parser.add_argument('-o', type=str, default='', help='output file (default: stdout)');
    parser.add_argument('-n', type=int, default=10, help='consider nearby variants within n bases apart as a set (off: 0, default: 10)');
    parser.add_argument('-k', type=int, default=10, help='maximum number of variants above which test 2^sqrt(n) instead of 2^n combinations (default: 10)');
    parser.add_argument('-mvh', action='store_true', help='consider nearby variants as one multi-variant hypothesis only');
    parser.add_argument('-p', action='store_true', help='consider only primary alignments');
    parser.add_argument('-t', type=int, default=1, help='number of processes to use (default: 1)');
    parser.add_argument('-debug', action='store_true', help='debug mode');
    if len(sys.argv) == 1:
        parser.print_help();
        sys.exit(1);
    args = parser.parse_args();

    global debug, primaryonly, numprocesses, distancethreshold, maxk, multivariant;
    debug = args.debug;
    primaryonly = args.p;
    numprocesses = args.t;
    distancethreshold = args.n;
    maxk = args.k;
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
