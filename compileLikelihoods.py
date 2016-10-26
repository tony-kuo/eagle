#!/usr/bin/python

# Copyright 2016 Tony Kuo
# This program is distributed under the terms of the GNU General Public License

# Utility Script
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
                ref = t[2].split(',');
                alt = t[3].split(',');
                # Account for double heterozygous non-reference or entries with the same position
                for i in ref:
                    for j in alt:
                        key = t[0]+'\t'+t[1]+'\t'+i+'\t'+j;
                        depth = int(t[5]) + int(t[6]);
                        if depth > 0:
                            af = float(t[6])/depth;
                            prob = float(t[7]);
                            lr = float(t[8]);
                            in_set = t[9];
                            if key not in entry: entry[key] = {};
                            if fn not in entry[key]: entry[key][fn] = [];
                            entry[key][fn].append((depth, af, lr, prob, in_set));
        fh.close;
    return(entry);

def compileEntries(entry, likelihood, af, keepone, isnegative):
    new_entry = {};
    for key in entry:
        # Check LR
        for fn in entry[key]: entry[key][fn] = sorted(entry[key][fn], key=lambda tup:tup[2], reverse=True)[0]; # Max LR per file for entries with same key
        current = sorted(entry[key].values(), key=lambda tup:tup[2], reverse=True)[0]; # Max LR across files
        valid = True;
        if isnegative and current[2] > likelihood: valid = False;
        elif not isnegative and current[2] < likelihood: valid = False;
        # Check AF
        current = sorted(entry[key].values(), key=lambda tup:tup[1], reverse=True)[0]; # Max AF across files
        if isnegative and current[1] > af: valid = False;
        elif not isnegative and current[1] < af: valid = False;

        # Positive samples require: LR >= threshold, AF <= threshold
        # Negative samples require: LR <= threshold, AF >= threshold
        if valid:
            if keepone: new_entry[key] = current;
            else: new_entry[key] = entry[key];
    return(new_entry);

def outputResults(pos_entry, neg_entry, somatic, pos_files, neg_files, pos_keepone, neg_keepone):
    header = 0;
    for key in naturalSort(pos_entry):
        if somatic:
            if not neg_entry or key not in neg_entry: continue;
            pos_prob = sorted(pos_entry[key].values(), key=lambda tup:tup[3], reverse=True)[0]; # Max probability across files
            neg_prob = sorted(neg_entry[key].values(), key=lambda tup:tup[3], reverse=True)[0]; # Max probability across files
            if np.power(10, pos_prob[3]) * (1-np.power(10, neg_prob[3])) < 0.99: 
                #print(pos_prob, neg_prob, np.power(10, pos_prob[3]) * (1-np.power(10, neg_prob[3])));
                continue;

        if header == 0:
            header = 1;
            outstr = '#\t\t\t\t';
            if not pos_keepone: outstr += '\t\t\t\t'.join(pos_files);
            else: outstr += '\t\t\t\t';
            if neg_files and not neg_keepone: outstr += '\t\t\t\t' + '\t\t\t\t'.join(neg_files);
            print(outstr);

        outstr = '{0}\t'.format(key);
        if pos_keepone: outstr += '{0}\t{1:.4}\t{2:.4}\t{3}\t'.format(pos_entry[key][0], pos_entry[key][1], pos_entry[key][2], pos_entry[key][4]);
        else:
            for fn in pos_files: 
                if fn in pos_entry[key]: outstr += '{0}\t{1:.4}\t{2:.4}\t{3}\t'.format(pos_entry[key][fn][0], pos_entry[key][fn][1], pos_entry[key][fn][2], pos_entry[key][fn][4]);
                else: outstr += '\t\t\t\t';
        if neg_entry and key in neg_entry:
            if neg_keepone: outstr += '{0}\t{1:.4}\t{2:.4}\t{3}'.format(neg_entry[key][0], neg_entry[key][1], neg_entry[key][2], neg_entry[key][4]);
            else:
                for fn in neg_files: 
                    fn = fn.split(',')[0];
                    if fn in neg_entry[key]: outstr += '{0}\t{1:.4}\t{2:.4}\t{3}\t'.format(neg_entry[key][fn][0], neg_entry[key][fn][1], neg_entry[key][fn][2], neg_entry[key][fn][4]);
                    else: outstr += '\t\t\t\t';
        print(outstr.strip());

def main():
    parser = argparse.ArgumentParser(description='Compile results from output of evalVariant [multiple, positive/negative]');
    parser.add_argument('-p', nargs='+', help='positive samples [f1 f2...]');
    parser.add_argument('-n', nargs='+', help='negative samples [f1 f2...]');
    parser.add_argument('-minlr', type=float, default=3, help='threshold for minimum log likelihood ratio for positive samples (default: 3)');
    parser.add_argument('-maxlr', type=float, default=-3, help='threshold for maximum log likelihood ratio for negative samples (default: -3)');
    parser.add_argument('-minaf', type=float, default=0.05, help='minimum allele frequency for positive samples (default: 0.05)');
    parser.add_argument('-maxaf', type=float, default=0.02, help='maximum allele frequency for negative samples (default: 0.02)');
    parser.add_argument('-s', action='store_true', help='somatic, it must exist in the negative sample and [Pr_positive * (1-Pr_negative) >= 0.99]');
    parser.add_argument('-p1', action='store_true', help='print only one positive entry among many (entry with max likelihood ratio)');
    parser.add_argument('-n1', action='store_true', help='print only one negative entry among many (entry with max likelihood ratio)');
    args = parser.parse_args();

    pos_entry = readFiles(args.p); 
    pos_entry = compileEntries(pos_entry, args.minlr, args.minaf, args.p1, False);
    neg_entry = {};
    if args.n:
        neg_entry = readFiles(args.n);
        neg_entry = compileEntries(neg_entry, args.maxlr, args.maxaf, args.n1, True);

    outputResults(pos_entry, neg_entry, args.s, args.p, args.n, args.p1, args.n1);

if __name__ == '__main__':
    try:
        main();
    except KeyboardInterrupt: 
        try:
            sys.exit(0);
        except SystemExit:
            os._exit(0);
