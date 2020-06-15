#!/usr/bin/python

# Copyright 2017 Tony Kuo
# This program is distributed under the terms of the GNU General Public License

# Utility Script
# Compiles a consensus classification based on the results of eagle-rc with --ngi mode, for hexaploid.
# Classification is determined by the reference with the highest probability.
# ex)
#   ~/eagle/eagle-rc --ngi --listonly --isc --ref1=$REF.fa.A.fa --ref2=$REF.fa.B.fa --bam1=$F.chrA.refsort.bam --bam2=$F.chrB.refsort.bam > $F.AvsB.list
#   ~/eagle/eagle-rc --ngi --listonly --isc --ref1=$REF.fa.A.fa --ref2=$REF.fa.D.fa --bam1=$F.chrA.refsort.bam --bam2=$F.chrD.refsort.bam > $F.AvsD.list
#   ~/eagle/eagle-rc --ngi --listonly --isc --ref1=$REF.fa.B.fa --ref2=$REF.fa.D.fa --bam1=$F.chrB.refsort.bam --bam2=$F.chrD.refsort.bam > $F.BvsD.list
#   python ref_consensus.py --pe -u -d -o sample -AB A.vs.B.list -AD A.vs.D.list -BD B.vs.D.list
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

def readFile(fn, entry1, entry2):
    with open(fn, 'r') as fh:
        for line in fh:
            if re.match('^#', line): continue

            t = line.strip().split('\t')
            if len(t) >= 8 and 'READ1' in t[7] and 'READ2' in t[7]:
                flag = 3
            elif len(t) >= 8 and 'READ2' in t[7]:
                flag = 2
            else:
                flag = 1

            key = '{}\t{}'.format(t[0], flag)
            pos = '{}\t{}'.format(t[2], t[3])

            # pos, prgu, total, pout, n
            if key not in entry1:
                entry1[key] = (pos, float(t[4]), float(t[5]), float(t[6]), 0)
            elif key in entry1:
                entry1[key] = (pos, np.logaddexp(entry1[key][1], float(t[4])), np.logaddexp(entry1[key][2], float(t[5])), np.logaddexp(entry1[key][3], float(t[6])), entry1[key][4] + 1)
            if key not in entry2:
                entry2[key] = (pos, float(t[5]), float(t[4]), float(t[6]), 0)
            elif key in entry2:
                entry2[key] = (pos, np.logaddexp(entry2[key][1], float(t[5])), np.logaddexp(entry2[key][2], float(t[4])), np.logaddexp(entry2[key][3], float(t[6])), entry2[key][4] + 1)
    fh.close
    print('Read:\t{}\t{}'.format(fn, datetime.now()), file=sys.stderr)
    return(entry1, entry2)

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
    parser.add_argument('-AB', required=True, help='AB classified list')
    parser.add_argument('-AD', required=True, help='AD classified list')
    parser.add_argument('-BD', required=True, help='BD classified list')
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
    chrB = {}
    chrD = {}
    (chrA, chrB) = readFile(args.AB, chrA, chrB)
    (chrA, chrD) = readFile(args.AD, chrA, chrD)
    (chrB, chrD) = readFile(args.BD, chrB, chrD)
    if args.pe:
        chrA = combinePE(chrA)
        chrB = combinePE(chrB)
        chrD = combinePE(chrD)
    writeTable(chrA, chrB, chrD, args.d, args.u, args.o)

if __name__ == '__main__':
    main()
    #os._exit(1)
