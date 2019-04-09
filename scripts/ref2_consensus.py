#!/usr/bin/python

# Copyright 2017 Tony Kuo
# This program is distributed under the terms of the GNU General Public License

# Utility Script
# Compiles a consensus classification based on the results of readclassify, for tetraploid.
# Classification is determined by the reference with the highest probability.
# ex)
#   eagle -t 8 -a chrA/refsort.bam -r A.fa -v AB.vcf --rc --splice 1> A.vs.B.txt 2> A.vs.B.readinfo.txt
#   eagle-rc --listonly -a chrA/refsort.bam -o A.vs.B -v A.vs.B.txt A.vs.B.readinfo.txt > A.vs.B.list
#
#   eagle -t 8 -a chrB/refsort.bam -r B.fa -v BA.vcf --rc --splice 1> B.vs.A.txt 2> B.vs.A.readinfo.txt
#   eagle-rc --listonly -a chrB/refsort.bam -o B.vs.A -v B.vs.A.txt B.vs.A.readinfo.txt > B.vs.A.list
#
#   python ref2_consensus.py --pe -u -o sample -A A.vs.B.list -B B.vs.A.list
#   eagle-rc --refonly --readlist -a chrA/refsort.bam -o sample.chrA sample.chrA.list
#   eagle-rc --refonly --readlist -a chrB/refsort.bam -o sample.chrB sample.chrB.list

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
            entry[key] = (pos, float(t[4]), float(t[5]), float(t[6]))
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
            entry[t[0]] = (entry[t[0]][0], entry[t[0]][1] + data[key][1], entry[t[0]][2] + data[key][2], entry[t[0]][3] + data[key][3])
    return(entry)

def classifySingle(key, chrA, fh, p_threshold):
    if chrA[key][1] - chrA[key][2] >= p_threshold: c = 'REF'
    else: c = 'UNK'
    t = key.strip().split('\t')
    if len(t) > 1: f = t[1]
    else: f = '-'
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(t[0], c, chrA[key][0], chrA[key][1], chrA[key][2], chrA[key][3], f), file=fh)

def writeTable(chrA, chrB, unique_reads, out_prefix):
    fhA = open(out_prefix + '.chrA.list', 'w')
    fhB = open(out_prefix + '.chrB.list', 'w')
    fh = [fhA, fhB]

    p_threshold = np.log(2)
    m_threshold = np.log(0.51)
    for key in chrA:
        if key not in chrB: continue

        pos = [chrA[key][0], chrB[key][0]]
        x = [chrA[key][1], chrB[key][1]]
        y = [chrA[key][2], chrB[key][2]]
        z = [chrA[key][3], chrB[key][3]] 

        p = [x[i] - y[i] for i in range(len(x))]
        total = logsumexp(p)
        m = [p[i] - total for i in range(len(p))]
        i = max(range(len(p)), key=p.__getitem__)

        if p[i] >= p_threshold and m[i] > m_threshold: c = 'REF'
        else: c = 'UNK'
        t = key.strip().split('\t')
        if len(t) > 1: f = t[1]
        else: f = '-'
        print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(t[0], c, pos[i], x[i], y[i], z[i], f), file=fh[i])

    if unique_reads:
        for key in chrA:
            if key not in chrB: classifySingle(key, chrA, fhA, p_threshold)
        for key in chrB:
            if key not in chrA: classifySingle(key, chrB, fhB, p_threshold)

    fhA.close()
    fhB.close()
    print('Done:\t{}'.format(datetime.now()), file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description='Determine read classification REF = A, B.')
    parser.add_argument('-A', nargs='+', required=True, help='2 list files: from readclassify with A as reference')
    parser.add_argument('-B', nargs='+', required=True, help='2 list files: from readclassify with B as reference')
    parser.add_argument('-o', type=str, required=True, help='output file prefix')
    parser.add_argument('-u', action='store_true', help='include reads that map uniquely to one reference genome')
    parser.add_argument('--pe', action='store_true', help='consider paired-end reads together')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    print('Start:\t{0}'.format(datetime.now()), file=sys.stderr)
    chrA = {}
    chrA = readFile(args.A[0], chrA)
    chrB = {}
    chrB = readFile(args.B[0], chrB)
    if args.pe:
        chrA = combinePE(chrA)
        chrB = combinePE(chrB)
    writeTable(chrA, chrB, args.u, args.o)

if __name__ == '__main__':
    main()
    #os._exit(1)#!/usr/bin/python

