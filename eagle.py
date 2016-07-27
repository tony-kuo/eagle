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
OMEGA = 1E-4; # Prior probability of read originating from an outside paralogous source
ALPHA = 1.3; # Factor to account for longer read lengths lowering the probability a sequence matching an outside paralogous source

# Precalculated log values 
LG3 = np.log10(3);
LG50 = np.log10(0.5);
LG10 = np.log10(0.1);
LG90 = np.log10(0.9);
LGOMEGA = np.log10(OMEGA);
LG1_OMEGA = np.log10(1-OMEGA);
LGALPHA = np.log10(ALPHA);

complement = { "A": "T", "C": "G", "G": "C", "T": "A", "a": "t", "c": "g", "g": "c", "t": "a" };

ffi = FFI();
ffi.cdef ("""
void initAlphaMap(void);
double log10addexp(double a, double b);
double log10sumexp(const double *n, size_t count);
void setProbMatrix(double *matrix, const char *seq, int readlength, const double *ismatch, const double *nomatch);
double calcProb(const char *seq, int seqlength, int pos, const double *matrix, int readlength, double baseline);
double calcProbDist(const char *seq, int seqlength, int pos, const double *matrix, int readlength);
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
    if (a > b) max_exp = a;
    else max_exp = b;
    return log10(pow(10, a - max_exp) + pow(10, b - max_exp)) + max_exp;
}
double log10sumexp(const double *n, size_t size) {
    size_t i;
    double max_exp = n[0]; 
    for (i = 1; i < size; i++) { 
        if (n[i] > max_exp) max_exp = n[i]; 
    }
    double sum = 0;
    for (i = 0; i < size; i++) sum += pow(10, n[i] - max_exp);
    return log10(sum) + max_exp;
}
void setProbMatrix(double *matrix, const char *seq, int readlength, const double *ismatch, const double *nomatch) {
    int i, b; // array[width * row + col] = value
    memset(matrix, 0, readlength * 5 * sizeof matrix);
    for (b = 0; b < readlength; ++b) {
        for (i = 0; i < 5; ++i) matrix[5 * b + i] = nomatch[b];
        matrix[5 * b +alphabetval[toupper(seq[b]) - 'A']] = ismatch[b];
    }
}
double calcProb(const char *seq, int seqlength, int pos, const double *matrix, int readlength, double baseline) {
    int b; // array[width * row + col] = value
    double probability = 0;
    for (b = pos;  b < pos+readlength; b++) {
        if (b < 0) continue;
        if (b >= seqlength) break;
        probability += matrix[5 * (b - pos) + alphabetval[toupper(seq[b]) - 'A']]; 
        if (probability < baseline - 10) break; // stop if less than 1% contribution to baseline (best/highest) probability mass
    }
    return probability;
}
double calcProbDist(const char *seq, int seqlength, int pos, const double *matrix, int readlength) {
    int i;
    double probability = 0;
    double baseline = calcProb(seq, seqlength, pos, matrix, readlength, -1000); // first probability at given pos, likely the highest, for initial baseline
    for (i = pos-readlength; i <= pos+readlength; ++i) {
        if (i + readlength < 0) continue;
        if (i >= seqlength) break;
        probability = probability == 0 ? calcProb(seq, seqlength, i, matrix, readlength, baseline) : log10addexp(probability, calcProb(seq, seqlength, i, matrix, readlength, baseline));
        if (probability > baseline) baseline = probability;
    }
    return probability;
}
""")
C.initAlphaMap(); # Initialize alphabet to int mapping table

def naturalSort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower();
    alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", key)];
    return (sorted(l, key=alphanum_key))

def readVCF(filename):
    entry = {};
    with open(filename, "r") as fh:
        for line in fh:
            if line[0] == "#" or line.strip() == "": continue;
            var = line.strip().split("\t");
            pos = int(var[1]);
            ref = var[3].upper().split(",");
            alt = var[4].upper().split(",");
            # Account for double heterozygous non-reference or entries with the same position
            for i in ref:
                for j in alt:
                    entry[(var[0], pos, i, j)] = [(pos, i, j)]; # key = (chr, pos, ref, alt)
    fh.close;
    print("Read VCF:\t{0}\t{1} entries\t{2}".format(filename, len(entry.keys()), datetime.now()), file=sys.stderr);
    return (entry);

def groupVariants(entry, distlim):
    # Check if current variant is within distance limit of the previous on the same chromosome, merge if so
    varid = sorted(entry.keys());
    skip2ind = -1;
    for i in range(0, len(varid)-1):
        if i <= skip2ind: continue;
        for j in range(i+1, len(varid)):
            if varid[j][0] != varid[j-1][0] or varid[j][1] - varid[j-1][1] > distlim: break;
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
                    misc[(varid[i][0], varid[i][1], varid[i][2], varid[i][3], j)] = list(entry[varid[i]]); # Copy the list of entries
                    del entry[varid[i]][j]; # Delete one of the same position entry from list
                    del misc[(varid[i][0], varid[i][1], varid[i][2], varid[i][3], j)][j+1]; # Delete the other same position entry from list
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
        #z = [j[1] for j in entry[varid[i]]] + [j[2] for j in entry[varid[i]]]; # Account for "-" representation
        #if x == 1 and y == 1 and "-" not in z:
            #ref = "".join([j[1] for j in entry[varid[i]]]);
            #alt = "".join([j[2] for j in entry[varid[i]]]);
            #misc[(varid[i][0], varid[i][1], ref, alt)] = [(varid[i][1], ref, alt)]; # Merge and delete extraneous entries
            #del entry[varid[i]];
    #if misc: entry.update(misc);
    #for i in entry:
        #if len(entry[i]) > maxk: print(i, len(entry[i]))
    #sys.exit(0);
    print("Group within:\t{0} bp\t{1} entries \t{2}".format(distlim, len(entry.keys()), datetime.now()), file=sys.stderr);
    return (entry);

def readFasta(filename):
    print("Read genome:\t{0}\t{1}".format(filename, datetime.now()), file=sys.stderr);
    seqid = "";
    seq = {};
    seqlength = {};
    with open(filename, "r") as fh:
        for line in fh:
            if line[0] == ">": #re.match("^>.+", line):
                seqid = re.split(">| ", line.strip())[1]; # >chr1 1 -> ["", "chr1", "1"]
                seq[seqid] = [];
            else:
                seq[seqid].append(line.strip().upper());
    fh.close;
    for seqid in seq: 
        seq[seqid] = "".join(seq[seqid]).encode("utf-8"); # Explicit byte string for python3, needed for cffi char*
        seqlength[seqid] = len(seq[seqid]);
    return (seq, seqlength);

def processWork(inqueue, results):
    while not inqueue.empty(): 
        args = inqueue.get();
        results.extend(evaluateVariant(args));

def processVariants(fn, var_list, outfile, numproc):
    print("Start:\t{0} procs \t{1}\t{2}".format(numproc, fn, datetime.now()), file=sys.stderr);

    results = [];
    varid = sorted(list(var_list.keys()));

    # Testing speed with Process, seems slower than Pool
    #manager = Manager();
    #work = manager.Queue(numproc);
    #results = manager.list();
    #pool = [];
    #for i in range(0,numproc):
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
            #results.extend(evaluateVariant((fn, varid[n], var_list[varid[n]]))); # single process, for testing and nicer error messages
    try:
        pool = Pool(processes=numproc);
        poolresults = pool.map_async(evaluateVariant, args);
        for p in poolresults.get(): results.extend(p);
    finally:
        pool.close();
        pool.join();

    fh = open(outfile, "w") if len(outfile) > 0 else sys.stdout;
    print("#SEQ\tPOS\tREF\tALT\tReads\tAltReads\tlog10_Prob\tlog10_Odds\tVarSet", file=fh);
    for i in naturalSort(results): print(i, file=fh);
    fh.close();
    print("Done:\t{0}\t{1}".format(fn, datetime.now()), file=sys.stderr);

#from readset import evaluateVariantCython; # Testing Cython: python setup.py build_ext --inplace
def evaluateVariant(args):
    #return (evaluateVariantCython(args, refseq, reflength, hetbias, maxh, multivariant, primaryonly, debug));
    (fn, varid, varset) = args;
    varstart = min([a[0] for a in varset]);
    varend = max([a[0] for a in varset]);

    # Fetch reads that cross the variant location
    samfile = pysam.AlignmentFile(fn, "rb");
    readset = list(samfile.fetch(reference=varid[0], start=max(0,varstart-1), end=min(varend+1,reflength[varid[0]])));

    v = [1] if len(varset) == 1 else [1, len(varset)];
    if not multivariant: v.extend(range(2, len(varset)));

    setid = 0;
    setentry = {};
    for currentset in chain(*map(lambda x: combinations(varset, x), v)): # For each set combination in powerset of variants in set excluding empty set
        if setid > 0 and setid >= maxh+len(varset)+1 and len(setentry[setid-1]) != len(currentset): break; # Stop if combinations limit reached and finished current k in n choose k
        setentry[setid] = currentset;
        setid += 1;

    refprior = np.log10(0.5);
    if len(varset) == 1 or multivariant: # Eeither one variant or multiple variants as a haplotype, homozygous & non-homozygous
        altprior = np.log10(0.5 * (1-hetbias));
        hetprior = np.log10(0.5 * hetbias);
    else:  # Variant set: divided evenly among the variant hypotheses, homozygous & non-homozygous
        altprior = np.log10(0.5 * (1-hetbias) / len(setentry));
        hetprior = np.log10(0.5 * hetbias / len(setentry));

    ref = 0.0;
    alt = {};
    het = {};
    refcount = {};
    altcount = {};
    prgu = {};
    pout = {};
    for setid in setentry:
        currentset = setentry[setid];
        alt[setid] = 0.0;
        het[setid] = 0.0;
        refcount[setid] = 0;
        altcount[setid] = 0;

        # Construct the variant sequence
        offset = 0;
        altseq = refseq[varid[0]];
        for v in sorted(currentset):
            pos = v[0] - 1 + offset;
            if v[1] == "-": # Account for "-" representation of variants, insertion
                pos += 1;
                varref = "";
                varalt = v[2];
            elif v[2] == "-": # Account for "-" representation of variants, deletion
                varalt = "";
                varref = v[1];
            else:
                varref = v[1];
                varalt = v[2];
            offset += len(varalt) - len(varref);
            altseq = altseq[:pos] + varalt.encode("utf-8") + altseq[(pos+len(varref)):]; # Explicit byte string for python3, needed for cffi char*
        altseqlength = len(altseq);

        for read in readset: 
            if read.is_unmapped: continue;
            if primaryonly and read.is_secondary: continue;

            readid = "{0}/1".format(read.query_name) if read.is_read1 else "{0}/2".format(read.query_name);

            # Convert read sequence into a probability matrix based on base-call error
            callerror = np.array(read.query_qualities) / -10.0; # Already converted to number by pysam
            callerror[callerror == 0] = -0.01; # Ensure there are no zeros, defaulting to a very high base-call error probability
            ismatch = np.log10(1 - np.power(10, callerror)); # log10(1-e);
            nomatch = callerror - LG3; #log10(e/3)

            readlength = len(read.query_sequence);
            readprobmatrix = np.zeros(readlength*5);
            p_readprobmatrix = ffi.cast("double *", readprobmatrix.ctypes.data); # Pointer to read probability matrix
            C.setProbMatrix(p_readprobmatrix, read.query_sequence.encode("utf-8"), readlength, list(ismatch), list(nomatch));

            # Reference genome probability and "elsewhere" probability only needs to be calculated once per readid
            if setid == 0:
                # Calculate the probability for "elsewhere", assuming the read is from somewhere paralogous. Approximate probability distribution by accounting for the bulk with:
                #   perfect match = prod[ (1-e) ]
                #   hamming/edit distance 1 = prod[ (1-e) ] * sum[ (e/3) / (1-e) ]
                # We also account for if reads have different lengths (hard clipped), where longer reads should have a relatively lower probability of originating from some paralogous elsewhere 
                #   lengthfactor = alpha ^ (readlength - expected readlength)
                # P(elsewhere) = (perfect + hamming) / lengthfactor
                elsewhere = C.log10addexp(sum(ismatch), sum(ismatch) + C.log10sumexp(list(nomatch - ismatch), len(ismatch))) - (LGALPHA * (readlength - read.infer_query_length())); 
                pout[readid] = elsewhere;
                # Calculate the probability given reference genome, once per varid for setid 0
                readprobability = C.calcProbDist(refseq[samfile.getrname(read.reference_id)], reflength[samfile.getrname(read.reference_id)], read.reference_start, p_readprobmatrix, readlength);
                prgu[readid] = readprobability;

            # Calculate the probability given alternate genome
            prgv = C.calcProbDist(altseq, altseqlength, read.reference_start, p_readprobmatrix, readlength);
            if debug: print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(setid, prgu[readid], prgv, pout[readid], read, currentset));

            # Multi-mapped alignments
            if not primaryonly and read.has_tag("XA"):
                for j in read.get_tag("XA").split(";"):
                    if len(j) <= 0: break; # Because last element is empty
                    t = j.split(",");
                    if t[0] not in refseq: continue;
                    xa_pos = int(t[1]);
                    if (read.is_reverse == False and xa_pos < 0) or (read.is_reverse == True and xa_pos > 0): # If aligned strand is opposite of that from primary alignment
                        newreadseq = "".join(complement.get(base,base) for base in reversed(read.query_sequence)).encode("utf-8"); # Reverse complement read sequence
                        newreadprobmatrix = np.zeros(readlength*5);
                        p_readprobmatrix = ffi.cast("double *", newreadprobmatrix.ctypes.data); # Pointer to read probability matrix
                        C.setProbMatrix(p_readprobmatrix, newreadseq, readlength, list(ismatch[::-1]), list(nomatch[::-1]));
                    else: 
                        p_readprobmatrix = ffi.cast("double *", readprobmatrix.ctypes.data); # Pointer to read probability matrix
                    xa_pos = abs(xa_pos);

                    readprobability = C.calcProbDist(refseq[t[0]], reflength[t[0]], xa_pos, p_readprobmatrix, readlength);
                    if setid == 0: 
                        # Probability given reference genome
                        prgu[readid] = C.log10addexp(prgu[readid], readprobability);
                        # The more multi-mapped, the more likely it is the read is from elsewhere (paralogous), hence it scales (multiplied) with the number of multi-mapped locations
                        pout[readid] = C.log10addexp(pout[readid], elsewhere);
                    if debug: print(readprobability, end="\t");
                    # Probability given alternate genome
                    if t[0] == samfile.getrname(read.reference_id): # If secondary alignments are in same chromosome (ie. contains variant thus has modified coordinates), in case it also crosses the variant position, otherwise is the same as probability given reference
                        if abs(xa_pos - varid[1]-1) < 50: readprobability = C.calcProbDist(altseq, altseqlength, xa_pos, p_readprobmatrix, readlength);
                    prgv = C.log10addexp(prgv, readprobability);
                    if debug: print("{0}\t{1}\t{2}".format(readprobability, prgu[readid], prgv));

            # Mixture model: probability that the read is from elsewhere, outside paralogous source
            if setid == 0: prgu[readid] = C.log10addexp(LGOMEGA - LG1_OMEGA + pout[readid], prgu[readid]);
            prgv = C.log10addexp(LGOMEGA - LG1_OMEGA + pout[readid], prgv);

            # Mixture model: heterozygosity or heterogeneity as explicit allele frequency mu such that P(r|GuGv) = (mu)(P(r|Gv)) + (1-mu)(P(r|Gu))
            phet = C.log10addexp(LG50 + prgv, LG50 + prgu[readid]);
            phet10 = C.log10addexp(LG10 + prgv, LG90 + prgu[readid]);
            phet90 = C.log10addexp(LG90 + prgv, LG10 + prgu[readid]);
            phet = max([phet, phet10, phet90]); # Use the best allele frequency probability

            # Read count is only incremented when the difference in probability is not ambiguous, ~log(2) difference
            if prgv > prgu[readid] and prgv - prgu[readid] > 0.3: altcount[setid] += 1;
            elif prgu[readid] > prgv and prgu[readid] - prgv > 0.3: refcount[setid] += 1;
            
            if setid == 0: ref += prgu[readid] + refprior; # Only one reference hypothesis
            alt[setid] += prgv + altprior;
            het[setid] += phet + hetprior;
            if debug: print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(prgu[readid], phet, prgv, pout[readid], altcount[setid], varid[0], readid, setentry[setid])); # log likelihoods
        if debug: print("{0}=\t{1}\t{2}\t{3}\t{4}\t{5}".format(setid, ref, het[setid], alt[setid], altcount[setid], setentry[setid])); # log likelihoods
    if not prgu: return ([]); # Return empty list if no read data

    outlist = [];
    total = [ref] + list(alt.values()) + list(het.values());
    total = C.log10sumexp(total, len(total));
    for v in varset:
        marginal_alt = 0.0;
        not_alt = ref;
        marginal_count = 0;
        for setid in alt:
            if v in setentry[setid]: 
                marginal_alt = C.log10addexp(alt[setid], het[setid]) if marginal_alt == 0.0 else C.log10addexp(marginal_alt, C.log10addexp(alt[setid], het[setid]));
                if altcount[setid] > marginal_count: marginal_count = altcount[setid];
            else: 
                not_alt = C.log10addexp(not_alt, C.log10addexp(alt[setid], het[setid]));
        outstr = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t".format(varid[0], v[0], v[1], v[2], max(refcount.values())+max(altcount.values()), marginal_count);
        # Probability and odds in log10
        outstr += "{0}\t{1}\t".format(marginal_alt - total, marginal_alt - not_alt); 
        # Print the variant set entries if exists    
        outstr += "{0}".format(varset) if len(varset) > 1 else "[]";
        outlist.append(outstr.strip());
    return (outlist);

# Global Variables
debug = False;
primaryonly = False;
hetbias = 0.5;
maxh = 1024;
multivariant = False;
refseq = {};
reflength = {};
def main():
    parser = argparse.ArgumentParser(description="Explicit alternative genome likelihood evaluator.");
    parser.add_argument("-v", help="variants VCF file");
    parser.add_argument("-a", help="alignment data bam files (index and ref coord sorted)");
    parser.add_argument("-r", help="reference sequence fasta file");
    parser.add_argument("-o", type=str, default="", help="output file (default: stdout)");
    parser.add_argument("-t", type=int, default=1, help="number of processes to use (default: 1)");
    parser.add_argument("-n", type=int, default=10, help="consider nearby variants within n bases in the set of hypotheses (off: 0, default: 10)");
    parser.add_argument("--maxh", type=int, default=1024, help="the maximum number of combinations in the set of hypotheses, instead of all 2^n (default: 2^10 = 1024)");
    parser.add_argument("--mvh", action="store_true", help="consider nearby variants as *one* multi-variant hypothesis");
    parser.add_argument("--hetbias", type=float, default=0.5, help="prior probability bias towards non-homozygous mutations (value between [0,1], default: 0.5 unbiased)");
    parser.add_argument("--pao", action="store_true", help="consider primary alignments only");
    parser.add_argument("--debug", action="store_true", help="debug mode, printing information on every read for every variant");
    if len(sys.argv) == 1:
        parser.print_help();
        sys.exit(1);
    args = parser.parse_args();

    global maxh, multivariant, hetbias, primaryonly, debug;
    numproc = max(1, args.t); # Main process counts as 1
    distlim = max(0, args.n);
    maxh = args.maxh;
    multivariant = args.mvh;
    hetbias = args.hetbias;
    primaryonly = args.pao;
    debug = args.debug;
    if hetbias < 0 or hetbias > 1: hetbias = 0.5;
    if maxh < 0: maxh = 1024;

    var_list = readVCF(args.v);
    if distlim > 0: var_list = groupVariants(var_list, distlim);

    global refseq, reflength;
    (refseq, reflength) = readFasta(args.r);
    processVariants(args.a, var_list, args.o, numproc);

if __name__ == "__main__":
    try:
        main();
    except KeyboardInterrupt: 
        try:
            sys.exit(0);
        except SystemExit:
            os._exit(0);
