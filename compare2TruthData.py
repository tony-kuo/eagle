#!/usr/bin/python

# 2016 Tony Kuo

# Utility script to compare results of EAGLE to VCF truth data
# True positives are written to stdout
# False positives are written to file

from __future__ import print_function
import argparse, re, sys, os;
from signal import signal, SIGPIPE, SIG_DFL;
signal(SIGPIPE,SIG_DFL) 

def naturalSort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)] 
    return sorted(l, key=alphanum_key)

def readFiles(filename):
    entry = {};
    with open(filename, 'r') as fh:
        for line in fh:
            if re.match('^#', line) or len(line.strip())==0: continue;
            t = line.strip().split('\t');
            key = '\t'.join(t[0:4]);
            entry[key] = line.strip();
    fh.close;
    return(entry);

def readVCF(filename):
    entry = {};
    with open(filename, 'r') as fh:
        for line in fh:
            if re.match('^#', line): continue;
            var = line.strip().split('\t');
            pos = int(var[1]);
            ref = var[3].split(',');
            alt = var[4].split(',');
            if var[0] not in entry: entry[var[0]] = {};
            # Account for double heterozygous non-reference or entries with the same position
            for i in ref:
                for j in alt:
                    (s1, s2) = removeCommonPrefix(i, j);
                    entry[var[0]][(pos, s1, s2)] = (pos, i, j);
    fh.close;
    return(entry);

def removeCommonPrefix(s1, s2):
    if len(s1) == 1 and len(s2) == 1: return(s1, s2); # Return as is if SNP
    i = 0;
    while (i < min(len(s1),len(s2)) and s1[i] == s2[i]): i += 1;
    if i > 0: return(s1[i-1:], s2[i-1:]);
    return(s1, s2); # Return as is if SNP

def outputResults(lr, vcf, false_filename, within_dist):
    with open(false_filename, 'w') as fh:
        for key in naturalSort(lr):
            t = key.strip().split('\t');
            refid = t[0];
            pos = int(t[1]);
            (ref, alt) = removeCommonPrefix(t[2], t[3]);
            if refid in vcf and within_dist > 0:
                pos_list = [a[0] for a in vcf[refid].keys()];
                if pos in pos_list: print('{0}'.format(lr[key]));
                elif min(pos_list, key=lambda x:abs(x-pos)) <= within_dist: print('{0}'.format(lr[key]));
                else: print('{0}'.format(lr[key]), file=fh);
            elif refid in vcf and (pos, ref, alt) in vcf[refid]: print('{0}'.format(lr[key]));
            else: print('{0}'.format(lr[key]), file=fh);
    fh.close;

def main():
    parser = argparse.ArgumentParser(description='Compare to VCF Truth Data');
    parser.add_argument('-l', help='output from EAGLE');
    parser.add_argument('-v', help='vcf truth data');
    parser.add_argument('-n', type=int, default=0, help='positional correctness only, within N bases is considered true (default: 0 [off])');
    parser.add_argument('-f', help='false positive file name');
    args = parser.parse_args();

    lr = readFiles(args.l); 
    vcf = readVCF(args.v); 
    outputResults(lr, vcf, args.f, args.n);

if __name__ == '__main__':
    try:
        main();
    except KeyboardInterrupt: 
        try:
            sys.exit(0);
        except SystemExit:
            os._exit(0);
