#!/usr/bin/python

# 2016 Tony Kuo

# Compiles multiple results from EAGLE
# Positive results must be > minimum for all samples
# Negative results must be < maximum for all samples
# A somatic result requires existence in both negative and positive samples

from __future__ import print_function
import argparse, re, sys, os;
import numpy as np;
from scipy.misc import logsumexp;
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 

def naturalSort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)] 
    return sorted(l, key=alphanum_key)

def readFiles(files):
    entry = {};
    for fn in files: 
        with open(fn, 'r') as fh:
            for line in fh:
                if re.match('^#', line) or len(line.strip())==0: continue;
                t = line.strip().split('\t');
                key = '\t'.join(t[0:4]);
                depth = int(t[4]);
                if depth > 0:
                    af = float(t[5])/depth;
                    prob = float(t[6]);
                    lr = float(t[7]);
                    in_set = t[8];
                    if key not in entry: 
                        entry[key] = {};
                        entry[key][fn] = [];
                    entry[key][fn].append((depth, af, lr, in_set));
        fh.close;
    return(entry);

def compileEntries(entry, threshold, keepone, isnegative):
    new_entry = {};
    for key in entry:
        for fn in entry[key]: entry[key][fn] = sorted(entry[key][fn], key=lambda tup:tup[2], reverse=True)[0]; # Max LR per file for entries with same key
        current = sorted(entry[key].values(), key=lambda tup:tup[2], reverse=True)[0]; # Max LR across files
        valid = False;
        if isnegative and current[2] < threshold and current[0] > 0 and current[1] <= 0: valid = True; # Negative samples require: LR < threshold, depth > 0, af <= 0
        elif not isnegative and current[2] > threshold and current[0] > 0 and current[1] > 0: valid = True; # Positive samples require: LR > threshold, depth and af > than 0
        if valid:
            if keepone: new_entry[key] = current;
            else: new_entry[key] = entry[key];
    return(new_entry);

def outputResults(pos_entry, neg_entry, somatic, pos_files, neg_files, pos_keepone, neg_keepone):
    header = 0;
    for key in naturalSort(pos_entry):
        if somatic and (not neg_entry or key not in neg_entry): continue;
        if not neg_entry or key in neg_entry:
            if header == 0:
                header = 1;
                outstr = '#\t\t\t\t';
                if not pos_keepone: outstr += '\t\t\t\t'.join(pos_files);
                else: outstr += '\t\t\t\t';
                if neg_files and not neg_keepone: outstr += '\t\t\t\t' + '\t\t\t\t'.join(neg_files);
                print(outstr);

            outstr = '{0}\t'.format(key);
            if pos_keepone: outstr += '{0}\t{1:.4}\t{2:.4}\t{3}\t'.format(pos_entry[key][0], pos_entry[key][1], pos_entry[key][2], pos_entry[key][3]);
            else:
                for fn in pos_files: 
                    if fn in pos_entry[key]: outstr += '{0}\t{1:.4}\t{2:.4}\t{3}\t'.format(pos_entry[key][fn][0], pos_entry[key][fn][1], pos_entry[key][fn][2], pos_entry[key][fn][3]);
                    else: outstr += '\t\t\t\t';
            if neg_files:
                if neg_keepone: outstr += '{0}\t{1:.4}\t{2:.4}\t{3}'.format(neg_entry[key][0], neg_entry[key][1], neg_entry[key][2], neg_entry[key][3]);
                else:
                    for fn in neg_files: 
                        fn = fn.split(',')[0];
                        if fn in neg_entry[key]: outstr += '{0}\t{1:.4}\t{2:.4}\t{3}\t'.format(neg_entry[key][fn][0], neg_entry[key][fn][1], neg_entry[key][fn][2], neg_entry[key][fn][3]);
                        else: outstr += '\t\t\t\t';
            print(outstr.strip());

def main():
    parser = argparse.ArgumentParser(description='Compile results from output of evalVariant [multiple, positive/negative]');
    parser.add_argument('-p', nargs='+', help='positive samples [f1 f2...]');
    parser.add_argument('-n', nargs='+', help='negative samples [f1 f2...]');
    parser.add_argument('-min', type=float, default=0, help='threshold for minimum log likelihood ratio for positive samples (default: 0)');
    parser.add_argument('-max', type=float, default=0, help='threshold for maximum log likelihood ratio for negative samples (default: 0)');
    parser.add_argument('-s', action='store_true', help='somatic, in that it must exist in the negative sample but log likelihood ratio < -max');
    parser.add_argument('-p1', action='store_true', help='print only one positive entry among many (entry with max likelihood ratio)');
    parser.add_argument('-n1', action='store_true', help='print only one negative entry among many (entry with max likelihood ratio)');
    args = parser.parse_args();

    pos_entry = readFiles(args.p); 
    pos_entry = compileEntries(pos_entry, args.mn, args.p1, False);
    neg_entry = {};
    if args.n:
        neg_entry = readFiles(args.n);
        neg_entry = compileEntries(neg_entry, args.mx, args.n1, True);

    outputResults(pos_entry, neg_entry, args.s, args.p, args.n, args.p1, args.n1);

if __name__ == '__main__':
    try:
        main();
    except KeyboardInterrupt: 
        try:
            sys.exit(0);
        except SystemExit:
            os._exit(0);
