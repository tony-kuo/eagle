#!/usr/bin/python

# Copyright 2017 Tony Kuo
# This program is distributed under the terms of the GNU General Public License

# Utility Script
# Compiles a consensus classification based on the results of readclassify, for hexaploid.
# Classification is determined by the reference with the highest probability.
# ex)
#   eagle -t 8 -a chrA/refsort.bam -r A.fa -v AB.vcf --rc --splice 1> A.vs.B.txt 2> A.vs.B.readinfo.txt
#   eagle-rc --listonly -a chrA/refsort.bam -o A.vs.B -v A.vs.B.txt A.vs.B.readinfo.txt > A.vs.B.list
#
#   eagle -t 8 -a chrA/refsort.bam -r A.fa -v AD.vcf --rc --splice 1> A.vs.D.txt 2> A.vs.D.readinfo.txt
#   eagle-rc --listonly -a chrA/refsort.bam -o A.vs.D -v A.vs.D.txt A.vs.D.readinfo.txt > A.vs.D.list
#
#   eagle -t 8 -a chrB/refsort.bam -r B.fa -v BA.vcf --rc --splice 1> B.vs.A.txt 2> B.vs.A.readinfo.txt
#   eagle-rc --listonly -a chrB/refsort.bam -o B.vs.A -v B.vs.A.txt B.vs.A.readinfo.txt > B.vs.A.list
#
#   eagle -t 8 -a chrB/refsort.bam -r B.fa -v BD.vcf --rc --splice 1> B.vs.D.txt 2> B.vs.D.readinfo.txt
#   eagle-rc --listonly -a chrB/refsort.bam -o B.vs.D -v B.vs.D.txt B.vs.D.readinfo.txt > B.vs.D.list
#
#   eagle -t 8 -a chrD/refsort.bam -r D.fa -v DA.vcf --rc --splice 1> D.vs.A.txt 2> D.vs.A.readinfo.txt
#   eagle-rc --listonly -a chrD/refsort.bam -o D.vs.A -v D.vs.A.txt D.vs.A.readinfo.txt > D.vs.A.list
#
#   eagle -t 8 -a chrD/refsort.bam -r D.fa -v DB.vcf --rc --splice 1> D.vs.B.txt 2> D.vs.B.readinfo.txt
#   eagle-rc --listonly -a chrD/refsort.bam -o D.vs.B -v D.vs.B.txt D.vs.B.readinfo.txt > D.vs.B.list
#
#   python ref_consensus.py --pe -u -d -o sample -A A.vs.B.list A.vs.D.list -B B.vs.A.list B.vs.D.list -D D.vs.A.list D.vs.B.list
#   eagle-rc --refonly --readlist -a chrA/refsort.bam -o sample.chrA sample.chrA.list
#   eagle-rc --refonly --readlist -a chrB/refsort.bam -o sample.chrB sample.chrB.list
#   eagle-rc --refonly --readlist -a chrD/refsort.bam -o sample.chrD sample.chrD.list

from __future__ import print_function
import argparse
import re, os, sys
import numpy as np
from scipy.special import logsumexp
from datetime import datetime
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 

def readFile(fn, entry):
    with open(fn, 'r') as fh:
        for line in fh:
            if re.match('^#', line): continue

            t = line.strip().split('\t')
            if len(t) < 8:
                t.append(' ')

            key = '{}\t{}'.format(t[0], t[7])
            pos = '{}\t{}'.format(t[2], t[3])

            if key not in entry: # pos, prgu, total, pout, n
                entry[key] = (pos, float(t[4]), float(t[5]), float(t[6]), 0)
            elif key in entry:
                entry[key] = (pos, entry[key][1], np.logaddexp(entry[key][2], float(t[5])), entry[key][3], entry[key][4] + 1)
    fh.close
    print('Read:\t{}\t{}'.format(fn, datetime.now()), file=sys.stderr)
    return(entry)

def combinePE(data):
    entry = {}
    for key in data:
        t = key.strip().split('\t')
        if t[0] not in entry:
            entry[t[0]] = data[key]
        elif t[0] in entry:
            entry[t[0]] = (entry[t[0]][0], entry[t[0]][1] + data[key][1], entry[t[0]][2] + data[key][2], entry[t[0]][3] + data[key][3], max(entry[t[0]][4], data[key][4]))
    return(entry)

def classifySingle(key, chrA, fh, p_threshold):
    if chrA[key][4] > 0 and chrA[key][1] - chrA[key][2] >= p_threshold: c = 'REF'
    else: c = 'UNK'
    t = key.strip().split('\t')
    if len(t) > 1: f = t[1]
    else: f = '-'
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(t[0], c, chrA[key][0], chrA[key][1], chrA[key][2], chrA[key][3], f), file=fh)

def classifyDouble(key, chrA, chrB, idx, fh, p_threshold, m_threshold):
    pos = [chrA[key][0], chrB[key][0]]
    x = [chrA[key][1], chrB[key][1]]
    y = [chrA[key][2], chrB[key][2]]
    z = [chrA[key][3], chrB[key][3]] 
    n = [chrA[key][4], chrB[key][4]] # number of pair-wise analysis that cross variants

    p = [x[i] - y[i] for i in range(len(x))]
    total = logsumexp(p)
    m = [p[i] - total for i in range(len(p))]
    i = max(range(len(p)), key=p.__getitem__)

    if n[i] > 0 and p[i] >= p_threshold and m[i] >= m_threshold: c = 'REF'
    else: c = 'UNK'
    t = key.strip().split('\t')
    if len(t) > 1: f = t[1]
    else: f = '-'
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(t[0], c, pos[i], x[i], y[i], z[i], f), file=fh[idx[i]])

def writeTable(chrA, chrB, chrD, doubles, unique_reads, out_prefix):
    fhA = open(out_prefix + '.chrA.list', 'w')
    fhB = open(out_prefix + '.chrB.list', 'w')
    fhD = open(out_prefix + '.chrD.list', 'w')
    fh = [fhA, fhB, fhD]

    p_threshold = np.log(2)
    m_threshold = np.log(0.51)
    for key in chrA:
        if key not in chrB or key not in chrD: continue

        pos = [chrA[key][0], chrB[key][0], chrD[key][0]]
        x = [chrA[key][1], chrB[key][1], chrD[key][1]]
        y = [chrA[key][2], chrB[key][2], chrD[key][2]]
        z = [chrA[key][3], chrB[key][3], chrD[key][3]]
        n = [chrA[key][4], chrB[key][4], chrD[key][4]] # number of pair-wise analysis that cross variants

        p = [x[i] - y[i] for i in range(len(x))]
        total = logsumexp(p)
        m = [p[i] - total for i in range(len(p))]
        i = max(range(len(p)), key=p.__getitem__)

        if n[i] > 0 and p[i] >= p_threshold and m[i] >= m_threshold: c = 'REF'
        else: c = 'UNK'
        t = key.strip().split('\t')
        if len(t) > 1: f = t[1]
        else: f = '-'
        print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(t[0], c, pos[i], x[i], y[i], z[i], f), file=fh[i])

    if doubles:
        for key in chrA:
            if key in chrB and key not in chrD: classifyDouble(key, chrA, chrB, (0, 1), fh, p_threshold, m_threshold)
            elif key not in chrB and key in chrD: classifyDouble(key, chrA, chrD, (0, 2), fh, p_threshold, m_threshold)
        for key in chrB:
            if key not in chrA and key in chrD: classifyDouble(key, chrB, chrD, (1, 2), fh, p_threshold, m_threshold)

    if unique_reads:
        for key in chrA:
            if key not in chrB and key not in chrD: classifySingle(key, chrA, fhA, p_threshold)
        for key in chrB:
            if key not in chrA and key not in chrD: classifySingle(key, chrB, fhB, p_threshold)
        for key in chrD:
            if key not in chrA and key not in chrB: classifySingle(key, chrD, fhD, p_threshold)

    fhA.close()
    fhB.close()
    fhD.close()
    print('Done:\t{}'.format(datetime.now()), file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description='Determine read classification REF = A, B, D.')
    parser.add_argument('-A', nargs='+', required=True, help='2 list files: from readclassify with A as reference')
    parser.add_argument('-B', nargs='+', required=True, help='2 list files: from readclassify with B as reference')
    parser.add_argument('-D', nargs='+', required=True, help='2 list files: from readclassify with D as reference')
    parser.add_argument('-o', type=str, required=True, help='output file prefix')
    parser.add_argument('-u', action='store_true', help='include reads that map uniquely to one subgenome')
    parser.add_argument('-d', action='store_true', help='classify reads that map to double copy homeologs')
    parser.add_argument('--pe', action='store_true', help='consider paired-end reads together')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    print('Start:\t{0}'.format(datetime.now()), file=sys.stderr)
    chrA = {}
    chrA = readFile(args.A[0], chrA) # file 1
    chrA = readFile(args.A[1], chrA) # file 2
    chrB = {}
    chrB = readFile(args.B[0], chrB) # file 1
    chrB = readFile(args.B[1], chrB) # file 2
    chrD = {}
    chrD = readFile(args.D[0], chrD) # file 1
    chrD = readFile(args.D[1], chrD) # file 2
    if args.pe:
        chrA = combinePE(chrA)
        chrB = combinePE(chrB)
        chrD = combinePE(chrD)
    writeTable(chrA, chrB, chrD, args.d, args.u, args.o)

if __name__ == '__main__':
    main()
    #os._exit(1)
