#!/usr/bin/python

# Copyright 2016 Tony Kuo
# This program is distributed under the terms of the GNU General Public License

# Utility Script
# Classifies reads as support reference, alternate, or reference-like-allele
# Can optionally provide a SAM file for splitting into the different classes

# Used in conjunction with EAGLE + debug mode -1 + mvh switches
# Recommended that you compile EAGLE with OMEGA=1E-40 to be less strict on the floor probability

# Example:
# eagle40 -t 8 -a reads.bam -r ref.fa -v variants.vcf --mvh --pao --isc -d -1 1> out.txt 2> out.readinfo.txt
# samtools view -h -F 4 reads.bam | python readClassify.py -o out -v out.txt -s - out.readinfo.txt > out.list

from __future__ import print_function;
import argparse;
import sys;
import numpy as np;
from re import split, match;
from datetime import datetime;
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 

def naturalSort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower();
    alphanum_key = lambda key: [convert(c) for c in split('([0-9]+)', key)];
    return(sorted(l, key=alphanum_key));

def readVar(fn):
    entry = [];
    with open(fn, 'r') as fh:
        for line in fh:
            if match('^#', line): continue;

            t = line.strip().split('\t');
            if t[9] == "[]":
                entry.extend(["{},{},{},{}".format(t[0], t[1], t[2], t[3])]);
            else:
                entry.extend(["{},{}".format(t[0], a) for a in split(';', t[9][1:-1])[0:-1]]);

    print("Read:\t{}\t{} variants\t{}".format(fn, len(entry), datetime.now()), file=sys.stderr);
    return(entry);

def readEAGLE(fn, varlist):
    #B_genome_base-36324     Chr1    999796  -70.812579      -60.964306      -92.313537      150M    0       PAIRED,PROPER_PAIR,REVERSE,READ2        Chr1,999943,C,T;
    fh = open(fn);
    tabletype = np.dtype([('rid','O'), ('chr','O'), ('pos','O'), ('prgu',float), ('prgv',float), ('pout',float), ('cigar','O'), ('multimap','O'), ('flags','O'), ('var','O')])
    table = np.loadtxt(fh, dtype=tabletype, delimiter='\t');

    entry = {};
    for i in range(0, len(table)):
        rid = table['rid'][i];
        if rid not in entry:
            entry[rid] = {};
            entry[rid]['var'] = {};
            entry[rid]['prgu'] = table['prgu'][i];
            entry[rid]['prgv'] = table['prgv'][i];
            entry[rid]['pout'] = table['pout'][i];
            entry[rid]['mult'] = False;
        else:
            entry[rid]['prgu'] = np.logaddexp(table['prgu'][i], entry[rid]['prgu']);
            entry[rid]['prgu'] = np.logaddexp(table['prgv'][i], entry[rid]['prgv']);
        entry[rid]['var'][table['var'][i]] = 0;
        
        if varlist and entry[rid]['mult'] == False: # check for the "ref-like allele" in multi-allelic variants
            for v in table['var'][i].split(';')[0:-1]:
                entry[rid]['mult'] = True;
                break;

    print("Read:\t{}\t{} entries\t{}".format(fn, len(entry), datetime.now()), file=sys.stderr);
    return(entry);

def classifyReads(entry):
    ref = {};
    alt = {};
    unk = {};
    mul = {};
    for rid in entry:
        var = entry[rid]['var'].keys();
        if entry[rid]['mult'] == True: # "ref-like allele", assumed to be closer to the reference due to smaller likelihood ratio over all reads
                ref[rid] = "{}\t{}\t{}\t{}".format(entry[rid]['prgu'], entry[rid]['prgv'], entry[rid]['pout'], " ".join(var));
                print("RLA=\t{}\t{}".format(rid, ref[rid]));
                continue;

        if len(var) == 2:
            pos = {};
            pos_alt = {};
            t = split('[,;]', "".join(var));
            for i in range(1, len(t), 4): 
                pos[t[i]] = 0;
                pos_alt["{},{}".format(t[i], t[i + 2])] = 0;
            if len(pos) != len(pos_alt): # duplicate position found
                mul[rid] = "{}\t{}\t{}\t{}".format(entry[rid]['prgu'], entry[rid]['prgv'], entry[rid]['pout'], " ".join(var));
                print("MUL=\t{}\t{}".format(rid, mul[rid]));
                continue;

        if (entry[rid]['prgu'] > entry[rid]['prgv'] and entry[rid]['prgu'] - entry[rid]['prgv'] >= 0.69):
            ref[rid] = "{}\t{}\t{}\t{}".format(entry[rid]['prgu'], entry[rid]['prgv'], entry[rid]['pout'], " ".join(var));
            print("REF=\t{}\t{}".format(rid, ref[rid]));
        elif (entry[rid]['prgv'] > entry[rid]['prgu'] and entry[rid]['prgv'] - entry[rid]['prgu'] >= 0.69):
            alt[rid] = "{}\t{}\t{}\t{}".format(entry[rid]['prgu'], entry[rid]['prgv'], entry[rid]['pout'], " ".join(var));
            print("ALT=\t{}\t{}".format(rid, alt[rid]));
        else:
            unk[rid] = "{}\t{}\t{}\t{}".format(entry[rid]['prgu'], entry[rid]['prgv'], entry[rid]['pout'], " ".join(var));
            print("UNK=\t{}\t{}".format(rid, unk[rid]));

    return(ref, alt, unk, mul);

def readLIST(fn):
    ref = {};
    alt = {};
    unk = {};
    mul = {};
    n = 0;
    with open(fn, 'r') as fh:
        for line in fh:
            if match('^#', line): continue;

            t = line.strip().split('\t');
            if t[0] == "REF=" or t[0] == "RLA":
                ref[t[1]] = line;
            elif t[0] == "ALT=":
                alt[t[1]] = line;
            elif t[0] == "UNK=":
                unk[t[1]] = line;
            elif t[0] == "MUL=":
                mul[t[1]] = line;
            n += 1;
    print("Read:\t{0}\t{1} entries".format(files, n), file=sys.stderr);
    return(ref, alt, unk, mul);

def classifySAM(sam, prefix, ref, alt, unk, mul):
    if sam == '-':
        fh = sys.stdin;
    else:
        fh = open(filename, 'r');

    fh_ref = open(prefix+".ref.sam", 'w');
    fh_alt = open(prefix+".alt.sam", 'w');
    fh_mul = open(prefix+".multiallele.sam", 'w');
    fh_common = open(prefix+".common.sam", 'w');
    for line in fh:
        if match('^@', line): 
            print(line, end='', file=fh_ref);
            print(line, end='', file=fh_alt);
            print(line, end='', file=fh_mul);
            print(line, end='', file=fh_common);
            continue;

        ID = line.strip().split('\t')[0];
        if ID in ref: 
            print(line, end='', file=fh_ref);
        elif ID in alt: 
            print(line, end='', file=fh_alt);
        elif ID in mul: 
            print(line, end='', file=fh_mul);
        elif ID not in unk:
            print(line, end='', file=fh_common);

    fh_ref.close();
    fh_alt.close();
    fh_mul.close();
    fh_common.close();
    if fh is not sys.stdin:
        fh.close();

def main():
    parser = argparse.ArgumentParser(description="Classify reads based on EAGLE read likelihoods from -d -1 mode");
    parser.add_argument('files', help="EAGLE read info file");
    parser.add_argument('-o', type=str, help="output files prefix for SAM files");
    parser.add_argument('-s', type=str, default="", help="corresponding SAM file of reads to be classified (- for piped input)");
    parser.add_argument('-v', type=str, default="", help="corresponding EAGLE results (with --mvh) for variant phase information");
    parser.add_argument('-l', action='store_true', help="input is previously processed read list file rather than raw EAGLE read info");
    if len(sys.argv) == 1:
        parser.print_help();
        sys.exit(1);
    args = parser.parse_args();

    if args.l:
        (ref, alt, unk, mul) = readLIST(args.files);
    else:
        varlist = [];
        if args.v:
            varlist = readVar(args.v);
        entry = readEAGLE(args.files, varlist);
        (ref, alt, unk, mul) = classifyReads(entry);
        
    if args.s:
        classifySAM(args.s, args.o, ref, alt, unk, mul);

if __name__ == '__main__':
    main()

