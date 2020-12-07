#!/usr/bin/python

# Copyright 2017 Tony Kuo
# This program is distributed under the terms of the GNU General Public License

# Utility Script
# Compiles a consensus classification based on the results of readclassify, for octoploid.
# Classification is determined by the reference with the highest probability.

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

    if n[i] > 1 and p[i] >= p_threshold and m[i] >= m_threshold: c = 'REF'
    else: c = 'UNK'
    t = key.strip().split('\t')
    if len(t) > 1: f = t[1]
    else: f = '-'
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(t[0], c, pos[i], x[i], y[i], z[i], f), file=fh[idx[i]])

def classifyTriple(key, chrA, chrB, chrC, idx, fh, p_threshold, m_threshold):
    pos = [chrA[key][0], chrB[key][0], chrC[key][0]]
    x = [chrA[key][1], chrB[key][1], chrC[key][1]]
    y = [chrA[key][2], chrB[key][2], chrC[key][2]]
    z = [chrA[key][3], chrB[key][3], chrC[key][3]]
    n = [chrA[key][4], chrB[key][4], chrC[key][4]] # number of pair-wise analysis that cross variants

    p = [x[i] - y[i] for i in range(len(x))]
    total = logsumexp(p)
    m = [p[i] - total for i in range(len(p))]
    i = max(range(len(p)), key=p.__getitem__)

    if n[i] > 1 and p[i] >= p_threshold and m[i] >= m_threshold: c = 'REF'
    else: c = 'UNK'
    t = key.strip().split('\t')
    if len(t) > 1: f = t[1]
    else: f = '-'
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(t[0], c, pos[i], x[i], y[i], z[i], f), file=fh[idx[i]])

def writeTable(chrA, chrB, chrC, chrD, triples, doubles, unique_reads, out_prefix):
    fhA = open(out_prefix + '.chrA.list', 'w')
    fhB = open(out_prefix + '.chrB.list', 'w')
    fhC = open(out_prefix + '.chrC.list', 'w')
    fhD = open(out_prefix + '.chrD.list', 'w')
    fh = [fhA, fhB, fhC, fhD]

    p_threshold = np.log(2)
    m_threshold = np.log(0.51)
    for key in chrA:
        if key not in chrB or key not in chrC or key not in chrD: continue

        pos = [chrA[key][0], chrB[key][0], chrC[key][0], chrD[key][0]]
        x = [chrA[key][1], chrB[key][1], chrC[key][1], chrD[key][1]]
        y = [chrA[key][2], chrB[key][2], chrC[key][2], chrD[key][2]]
        z = [chrA[key][3], chrB[key][3], chrC[key][3], chrD[key][3]]
        n = [chrA[key][4], chrB[key][4], chrC[key][4], chrD[key][4]] # number of pair-wise analysis that cross variants

        p = [x[i] - y[i] for i in range(len(x))]
        total = logsumexp(p)
        m = [p[i] - total for i in range(len(p))]
        i = max(range(len(p)), key=p.__getitem__)

        if n[i] > 1 and p[i] >= p_threshold and m[i] >= m_threshold: c = 'REF'
        else: c = 'UNK'
        t = key.strip().split('\t')
        if len(t) > 1: f = t[1]
        else: f = '-'
        print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(t[0], c, pos[i], x[i], y[i], z[i], f), file=fh[i])

    if triples:
        for key in chrA:
            if key in chrB and key in chrC and key not in chrD: classifyTriple(key, chrA, chrB, chrC, (0, 1, 2), fh, p_threshold, m_threshold)
            if key in chrB and key not in chrC and key in chrD: classifyTriple(key, chrA, chrB, chrD, (0, 1, 3), fh, p_threshold, m_threshold)
            if key not in chrB and key in chrC and key in chrD: classifyTriple(key, chrA, chrC, chrD, (0, 2, 3), fh, p_threshold, m_threshold)
        for key in chrB:
            if key not in chrA and key in chrC and key in chrD: classifyTriple(key, chrB, chrC, chrD, (1, 2, 3), fh, p_threshold, m_threshold)

    if doubles:
        for key in chrA:
            if key in chrB and key not in chrC and key not in chrD: classifyDouble(key, chrA, chrB, (0, 1), fh, p_threshold, m_threshold)
            elif key not in chrB and key in chrC and key not in chrD: classifyDouble(key, chrA, chrC, (0, 2), fh, p_threshold, m_threshold)
            elif key not in chrB and key not in chrC and key in chrD: classifyDouble(key, chrA, chrD, (0, 3), fh, p_threshold, m_threshold)
        for key in chrB:
            if key not in chrA and key in chrC and key not in chrD: classifyDouble(key, chrB, chrC, (1, 2), fh, p_threshold, m_threshold)
            elif key not in chrA and key not in chrC and key in chrD: classifyDouble(key, chrB, chrD, (1, 3), fh, p_threshold, m_threshold)
        for key in chrC:
            if key not in chrA and key not in chrB and key in chrD: classifyDouble(key, chrC, chrD, (2, 3), fh, p_threshold, m_threshold)

    if unique_reads:
        for key in chrA:
            if key not in chrB and key not in chrC and key not in chrD: classifySingle(key, chrA, fhA, p_threshold)
        for key in chrB:
            if key not in chrA and key not in chrC and key not in chrD: classifySingle(key, chrB, fhB, p_threshold)
        for key in chrC:
            if key not in chrA and key not in chrB and key not in chrD: classifySingle(key, chrC, fhC, p_threshold)
        for key in chrD:
            if key not in chrA and key not in chrB and key not in chrC: classifySingle(key, chrD, fhD, p_threshold)

    fhA.close()
    fhB.close()
    fhC.close()
    fhD.close()
    print('Done:\t{}'.format(datetime.now()), file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description='Determine read classification REF = A, B, D.')
    parser.add_argument('-A', nargs='+', required=True, help='3 list files: from readclassify with A as reference')
    parser.add_argument('-B', nargs='+', required=True, help='3 list files: from readclassify with B as reference')
    parser.add_argument('-C', nargs='+', required=True, help='3 list files: from readclassify with C as reference')
    parser.add_argument('-D', nargs='+', required=True, help='3 list files: from readclassify with D as reference')
    parser.add_argument('-o', type=str, required=True, help='output file prefix')
    parser.add_argument('-u', action='store_true', help='include reads that map uniquely to one subgenome')
    parser.add_argument('-d', action='store_true', help='classify reads that map to double copy homeologs')
    parser.add_argument('-t', action='store_true', help='classify reads that map to triple copy homeologs')
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
    chrA = readFile(args.A[2], chrA) # file 3
    chrB = {}
    chrB = readFile(args.B[0], chrB) # file 1
    chrB = readFile(args.B[1], chrB) # file 2
    chrB = readFile(args.B[2], chrB) # file 3
    chrC = {}
    chrC = readFile(args.C[0], chrC) # file 1
    chrC = readFile(args.C[1], chrC) # file 2
    chrC = readFile(args.C[2], chrC) # file 3
    chrD = {}
    chrD = readFile(args.D[0], chrD) # file 1
    chrD = readFile(args.D[1], chrD) # file 2
    chrD = readFile(args.D[2], chrD) # file 3
    if args.pe:
        chrA = combinePE(chrA)
        chrB = combinePE(chrB)
        chrC = combinePE(chrC)
        chrD = combinePE(chrD)
    writeTable(chrA, chrB, chrC, chrD, args.t, args.d, args.u, args.o)

if __name__ == '__main__':
    main()
    #os._exit(1)
