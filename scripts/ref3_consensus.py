#!/usr/bin/python

# Copyright 2017 Tony Kuo
# This program is distributed under the terms of the GNU General Public License

# Utility Script
# Compiles a consensus classification based on the results of readclassify, for hexaploid.
# Classification is determined by the reference with the highest probability.
# ex)
#   eagle -t 8 -a chrA/refsort.bam -r A.fa -v AB.vcf --omega=1e-40 --mvh --isc --splice --verbose 1> A.vs.B.txt 2> A.vs.B.readinfo.txt
#   eagle-rc -a chrA/refsort.bam --listonly -o A.vs.B A.vs.B.txt A.vs.B.readinfo.txt > A.vs.B.list
#
#   eagle -t 8 -a chrA/refsort.bam -r A.fa -v AD.vcf --omega=1e-40 --mvh --isc --splice --verbose 1> A.vs.D.txt 2> A.vs.D.readinfo.txt
#   eagle-rc -a chrA/refsort.bam --listonly -o A.vs.D A.vs.D.txt A.vs.D.readinfo.txt > A.vs.D.list
#
#   eagle -t 8 -a chrB/refsort.bam -r B.fa -v BA.vcf --omega=1e-40 --mvh --isc --splice --verbose 1> B.vs.A.txt 2> B.vs.A.readinfo.txt
#   eagle-rc -a chrB/refsort.bam --listonly -o B.vs.A B.vs.A.txt B.vs.A.readinfo.txt > B.vs.A.list
#
#   eagle -t 8 -a chrB/refsort.bam -r B.fa -v BD.vcf --omega=1e-40 --mvh --isc --splice --verbose 1> B.vs.D.txt 2> B.vs.D.readinfo.txt
#   eagle-rc -a chrB/refsort.bam --listonly -o B.vs.D B.vs.D.txt B.vs.D.readinfo.txt > B.vs.D.list
#
#   eagle -t 8 -a chrD/refsort.bam -r D.fa -v DA.vcf --omega=1e-40 --mvh --isc --splice --verbose 1> D.vs.A.txt 2> D.vs.A.readinfo.txt
#   eagle-rc -a chrD/refsort.bam --listonly -o D.vs.A D.vs.A.txt D.vs.A.readinfo.txt > D.vs.A.list
#
#   eagle -t 8 -a chrD/refsort.bam -r D.fa -v DB.vcf --omega=1e-40 --mvh --isc --splice --verbose 1> D.vs.B.txt 2> D.vs.B.readinfo.txt
#   eagle-rc -a chrD/refsort.bam --listonly -o D.vs.B D.vs.B.txt D.vs.B.readinfo.txt > D.vs.B.list
#
#   python ref_consensus.py -o sample -A A.vs.B.list A.vs.D.list -B B.vs.A.list B.vs.D.list -D D.vs.A.list D.vs.B.list
#   eagle-rc -a chrA/refsort.bam --refonly --readlist -o sample.chrA sample.chrA.list
#   eagle-rc -a chrB/refsort.bam --refonly --readlist -o sample.chrB sample.chrB.list
#   eagle-rc -a chrD/refsort.bam --refonly --readlist -o sample.chrD sample.chrD.list

from __future__ import print_function
import argparse
import re, os, sys
import numpy as np
from datetime import datetime
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 

def readFile(fn, entry, n):
    with open(fn, 'r') as fh:
        for line in fh:
            if re.match('^#', line): continue

            t = line.strip().split('\t')
            key = t[0]

            if key not in entry: 
                entry[key] = (float(t[4]), np.logaddexp(np.logaddexp(float(t[4]), float(t[5])), float(t[6])), 0) # prgu + prgv + pout
            elif key in entry:
                prgu = entry[key][0] + float(t[4]) # ref[N0] * ref[N1]
                #prgv = entry[key][1] + float(t[5]) # alt[N0] * alt[N1]
                total = entry[key][1] + np.logaddexp(np.logaddexp(float(t[4]), float(t[5])), float(t[6])) # total[N0] * total[N1]
                entry[key] = (prgu, total, n)
    fh.close
    print("Read:\t{}\t{}".format(fn, datetime.now()), file=sys.stderr)
    return(entry)

def writeTable(chrA, chrB, chrD, unique_reads, out_prefix):
    fhA = open(out_prefix + '.chrA.list', 'w')
    fhB = open(out_prefix + '.chrB.list', 'w')
    fhD = open(out_prefix + '.chrD.list', 'w')
    fh = [fhA, fhB, fhD]

    l2 = np.log(2)
    threshold = np.log(0.95)
    for key in chrA:
        if key not in chrB or key not in chrD: continue

        x = [chrA[key][0], chrB[key][0], chrD[key][0]] # numerator
        y = [chrA[key][1], chrB[key][1], chrD[key][1]] # denominator
        z = [chrA[key][2], chrB[key][2], chrD[key][2]]

        for i in range(len(z)):
            if z[i] == 0: y[i] += l2 # divide by 2 if a pair-wise analysis did not cross variants, ie. common between a pairing = ambiguous

        i = max(range(len(x)), key=x.__getitem__)
        #p = [x[0] - y[0], x[1] - y[1], x[2] - y[2]]
        #d = [np.exp(p[i]) - np.exp(p[j]) for j in range(len(p)) if i != j]

        if x[i] - y[i] > threshold: c = "REF"
        else: c = "UNK"
        print("{}\t{}\t-\t-\t{}\t{}\t-".format(key, c, x[i], y[i]), file=fh[i])

    if unique_reads:
        for key in chrA:
            if key in chrB or key in chrD: continue
            x = chrA[key][0]
            y = chrA[key][1]
            if chrA[key][2] == 0:
                y += l2 # divide by 2 if a pair-wise analysis did not cross variants, ie. common between a pairing = ambiguous
            if x - y > threshold: c = "REF"
            else: c = "UNK"
            print("{}\t{}\t-\t-\t{}\t{}\t-".format(key, c, x, y), file=fhA)

        for key in chrB:
            if key in chrA or key in chrD: continue
            x = chrB[key][0]
            y = chrB[key][1]
            if chrB[key][2] == 0:
                y += l2 # divide by 2 if a pair-wise analysis did not cross variants, ie. common between a pairing = ambiguous
            if x - y > threshold: c = "REF"
            else: c = "UNK"
            print("{}\t{}\t-\t-\t{}\t{}\t-".format(key, c, x, y), file=fhB)

        for key in chrD:
            if key in chrA or key in chrB: continue
            x = chrD[key][0]
            y = chrD[key][1]
            if chrD[key][2] == 0:
                y += l2 # divide by 2 if a pair-wise analysis did not cross variants, ie. common between a pairing = ambiguous
            if x - y > threshold: c = "REF"
            else: c = "UNK"
            print("{}\t{}\t-\t-\t{}\t{}\t-".format(key, c, x, y), file=fhD)

    fhA.close()
    fhB.close()
    fhD.close()
    print("Done:\t{}".format(datetime.now()), file=sys.stderr)

def main():
    #python ref_consensus.py -o $F.ref -A $F.A.vs.B.con.list $F.A.vs.D.con.list -B $F.B.vs.A.con.list $F.B.vs.D.con.list -D $F.D.vs.A.con.list $F.D.vs.A.con.list
    parser = argparse.ArgumentParser(description='Determine read classification REF = A, B, D.  Classification is determined by log likelihood ratio')
    parser.add_argument('-A', nargs='+', required=True, help='2 list files: from readclassify with A as reference followed by mirror consensus')
    parser.add_argument('-B', nargs='+', required=True, help='2 list files: from readclassify with B as reference followed by mirror consensus')
    parser.add_argument('-D', nargs='+', required=True, help='2 list files: from readclassify with D as reference followed by mirror consensus')
    parser.add_argument('-o', type=str, required=True, help='output file prefix')
    parser.add_argument('-u', action='store_true', help='include reads that map uniquely to one reference genome')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    print("Start:\t{0}".format(datetime.now()), file=sys.stderr)
    chrA = {}
    chrA = readFile(args.A[0], chrA, 0) # file 1
    chrA = readFile(args.A[1], chrA, 1) # file 2
    chrB = {}
    chrB = readFile(args.B[0], chrB, 0) # file 1
    chrB = readFile(args.B[1], chrB, 1) # file 2
    chrD = {}
    chrD = readFile(args.D[0], chrD, 0) # file 1
    chrD = readFile(args.D[1], chrD, 1) # file 2
    writeTable(chrA, chrB, chrD, args.u, args.o)

if __name__ == '__main__':
    main()
    #os._exit(1)
