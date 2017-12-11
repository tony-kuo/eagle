#!/usr/bin/python

# Copyright 2017 Tony Kuo
# This program is distributed under the terms of the GNU General Public License

# Utility Script
# Compiles a consensus classification based on the results of readclassify using mirrored references.
# Classification is determined by the log likelihood sum over all variants that a read crosses.
# ex)
#   eagle -t 8 -a refsort.bam -r A.fa -v AB.vcf --omega=1e-40 --mvh --isc --splice --verbose 1> A.vs.B.txt 2> A.vs.B.readinfo.txt
#   readclassify -a refsort.bam --listonly -o A.vs.B A.vs.B.txt A.vs.B.readinfo.txt > A.vs.B.list
#
#   eagle -t 8 -a refsort.bam -r B.fa -v BA.vcf --omega=1e-40 --mvh --isc --splice --verbose 1> B.vs.A.txt 2> B.vs.A.readinfo.txt
#   readclassify -a refsort.bam --listonly -o B.vs.A B.vs.A.txt B.vs.A.readinfo.txt > B.vs.A.list
#
#   python mirror_consensus.py A.vs.B.list B.vs.A.list > consensus.A.vs.B.list
#   readclassify -a refsort.bam -o consensus_A --readlist consensus.A.list
#
#   python mirror_consensus.py B.vs.A.list A.vs.B.list > consensus.B.vs.A.list
#   readclassify -a refsort.bam -o consensus_A --readlist consensus.B.list

from __future__ import print_function
import argparse
import re
import os, sys
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
                entry[key] = (float(t[4]), float(t[5]), 0)
            elif key in entry and n > 0:
                prgu = float(t[5]) + entry[key][0] # alt[N1] + ref[N0]
                prgv = float(t[4]) + entry[key][1] # ref[N1] + alt[N0]
                entry[key] = (prgu, prgv, n)
    fh.close
    print("Read:\t{}\t{}".format(fn, datetime.now()), file=sys.stderr)
    return(entry)

def writeTable(entry):
    for key in entry:
        if entry[key][2] > 0:
            if entry[key][0] > entry[key][1] and entry[key][0] - entry[key][1] > 0.69: c = "REF"
            elif entry[key][1] > entry[key][0] and entry[key][1] - entry[key][0] > 0.69: c = "ALT"
            else: c = "UNK"
            print("{}\t{}\t-\t-\t{}\t{}\t-".format(key, c, entry[key][0], entry[key][1]))
    print("Done:\t{}".format(datetime.now()), file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description='Consensus of paired read classifications with mirrored references (REF1 = ALT2 && ATL1 = REF2).  Classification is determined by the log likelihood sum over all variants that a read crosses.  REF uses the first file as the template.')
    parser.add_argument('files', nargs='+', help='2 list files from readclassify')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    print("Start:\t{0}".format(datetime.now()), file=sys.stderr)
    entry = {}
    entry = readFile(args.files[0], entry, 0) # file 1
    entry = readFile(args.files[1], entry, 1) # file 2
    writeTable(entry)

if __name__ == '__main__':
    main()
    #os._exit(1)
