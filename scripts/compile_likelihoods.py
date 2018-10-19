#!/usr/bin/python

# Copyright 2016 Tony Kuo
# This program is distributed under the terms of the GNU General Public License

# Utility Script
# Compiles multiple results from EAGLE
# Positive results must be > minimum for all samples
# Negative results must be < maximum for all samples
# A somatic result requires existence in both negative and positive samples

from __future__ import print_function
import argparse, re, sys, os
import numpy as np
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 

def naturalSort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)] 
    return sorted(l, key=alphanum_key)

def readFiles(files, reads_seen, minrp):
    entry = {}
    for fn in files: 
        with open(fn, 'r') as fh:
            for line in fh:
                if re.match('^#', line) or len(line.strip())==0: continue
                t = line.strip().split('\t')
                ref = t[2].split(',')
                alt = t[3].split(',')
                if int(t[5]) == 0 and int(t[6]) == 0: continue
                if (float(t[5]) + float(t[6])) / float(t[4]) < minrp: continue
                # Account for double heterozygous non-reference or entries with the same position
                for i in ref:
                    for j in alt:
                        key = t[0]+'\t'+t[1]+'\t'+i+'\t'+j
                        if reads_seen:
                            depth = int(t[4])
                        else:
                            depth = int(t[5]) + int(t[6])

                        if depth > 0:
                            af = float(t[6])/depth
                            prob = float(t[7])
                            lr = float(t[8])
                            in_set = t[9]
                            if key not in entry: entry[key] = {}
                            if fn not in entry[key]: entry[key][fn] = [(depth, af, lr, prob, in_set)]
                            else: entry[key][fn].append((depth, af, lr, prob, in_set))
        fh.close
    return(entry)

def compileEntries(entry, likelihood, af, mindepth, maxdepth, isnegative):
    new_entry = {}
    for key in entry:
        # Check LR
        for fn in entry[key]: entry[key][fn] = sorted(entry[key][fn], key=lambda tup:tup[2], reverse=True)[0] # Max LR per file for entries with same key
        current = sorted(entry[key].values(), key=lambda tup:tup[2], reverse=True)[0] # Max LR across files
        
        if isnegative and current[2] > likelihood: continue
        elif not isnegative and current[2] < likelihood: continue
        # Check AF
        current = sorted(entry[key].values(), key=lambda tup:tup[1], reverse=True)[0] # Max AF across files
        if isnegative and current[1] > af: continue
        elif not isnegative and current[1] < af: continue
        # Check read min depth
        current = sorted(entry[key].values(), key=lambda tup:tup[1])[0] # Min Read Depth across files
        if mindepth > 0 and current[0] < mindepth: continue
        # Check read max depth
        current = sorted(entry[key].values(), key=lambda tup:tup[1], reverse=True)[0] # Max Read Depth across files
        if maxdepth > 0 and current[0] > maxdepth: continue

        # Positive samples require: LR >= threshold, AF <= threshold
        # Negative samples require: LR <= threshold, AF >= threshold
        new_entry[key] = entry[key]
    return(new_entry)

def compileLOH(pos_entry, neg_entry, minlr, maxlr, depth):
    new_pos_entry = {}
    new_neg_entry = {}
    for key in pos_entry:
        if key not in neg_entry: continue

        current = sorted(pos_entry[key].values(), key=lambda tup:tup[2], reverse=True)[0] # Max LR across files
        if current[2] < minlr: continue
        if current[1] < 0.3 or current[1] > 0.7: continue # pos should be heterozygous variant
        current = sorted(pos_entry[key].values(), key=lambda tup:tup[1])[0] # Min Read Depth across files
        if depth > 0 and current[0] < depth: continue

        current = sorted(neg_entry[key].values(), key=lambda tup:tup[2], reverse=True)[0] # Max LR across files
        if maxlr < current[2] < minlr: continue # neg should be homozygous variant or reference, thus should be < maxlr or > minlr
        if current[1] >= 0.3 or current[1] <= 0.7: continue # neg should be homozygous variant or reference
        current = sorted(neg_entry[key].values(), key=lambda tup:tup[1])[0] # Min Read Depth across files
        if depth > 0 and current[0] < depth: continue

        new_pos_entry[key] = pos_entry[key]
        new_neg_entry[key] = neg_entry[key]
    return(new_pos_entry, new_neg_entry)

def outputResults(pos_entry, neg_entry, pos_files, neg_files, mincp):
    header = True
    for key in naturalSort(pos_entry):
        if neg_files:
            if not neg_entry or key not in neg_entry: continue
            pos_prob = sorted(pos_entry[key].values(), key=lambda tup:tup[3], reverse=True)[0] # Max probability across files
            neg_prob = sorted(neg_entry[key].values(), key=lambda tup:tup[3], reverse=True)[0] # Max probability across files
            if np.power(10, pos_prob[3]) * (1 - np.power(10, neg_prob[3])) < mincp: 
                #print(key, pos_prob, neg_prob, np.power(10, pos_prob[3]) * (1-np.power(10, neg_prob[3])))
                continue

        if header:
            header = False
            outstr = '#\t\t\t\t' + '\t\t\t\t'.join(pos_files)
            if neg_files: outstr += '\t\t\t\t' + '\t\t\t\t'.join(neg_files)
            print(outstr)

        outstr = '{0}\t'.format(key)
        for fn in pos_files: 
            if fn in pos_entry[key]: outstr += '{0}\t{1:.4}\t{2:.4}\t{3}\t'.format(pos_entry[key][fn][0], pos_entry[key][fn][1], pos_entry[key][fn][2], pos_entry[key][fn][4])
            else: outstr += '\t\t\t\t'
        if neg_entry and key in neg_entry:
            for fn in neg_files: 
                if fn in neg_entry[key]: outstr += '{0}\t{1:.4}\t{2:.4}\t{3}\t'.format(neg_entry[key][fn][0], neg_entry[key][fn][1], neg_entry[key][fn][2], neg_entry[key][fn][4])
                else: outstr += '\t\t\t\t'
            outstr += '\tSOM'
        print(outstr.strip())

def outputLOH(pos_entry, neg_entry, pos_files, neg_files):
    for key in naturalSort(pos_entry):
        if key not in neg_entry: continue

        outstr = '{0}\t'.format(key)
        for fn in pos_files: 
            if fn in pos_entry[key]: outstr += '{0}\t{1:.4}\t{2:.4}\t{3}\t'.format(pos_entry[key][fn][0], pos_entry[key][fn][1], pos_entry[key][fn][2], pos_entry[key][fn][4])
            else: outstr += '\t\t\t\t'
        for fn in neg_files: 
            if fn in neg_entry[key]: outstr += '{0}\t{1:.4}\t{2:.4}\t{3}\t'.format(neg_entry[key][fn][0], neg_entry[key][fn][1], neg_entry[key][fn][2], neg_entry[key][fn][4])
            else: outstr += '\t\t\t\t'
        outstr += '\tLOH'
        print(outstr.strip())

def main():
    parser = argparse.ArgumentParser(description='Compile results from output of EAGLE [multiple, positive/negative]. If negative samples provided, somatic mutations are compiled where data for variant must exist in the negative sample and not support the variant [so Pr_positive * (1-Pr_negative) >= 0.99]')
    parser.add_argument('-p', nargs='+', help='positive samples [f1 f2...]')
    parser.add_argument('-n', nargs='+', help='negative samples [f1 f2...]')
    parser.add_argument('-minlr', type=float, default=5, help='threshold for minimum log likelihood ratio for positive samples [default: 5]')
    parser.add_argument('-maxlr', type=float, default=-2, help='threshold for maximum log likelihood ratio for negative samples [default: -2]')
    parser.add_argument('-minaf', type=float, default=0.1, help='minimum allele frequency for positive samples [default: 0.1]')
    parser.add_argument('-maxaf', type=float, default=0.02, help='maximum allele frequency for negative samples [default: 0.02]')
    parser.add_argument('-mindepth', type=int, default=0, help='minimum read depth [default: 0 off]')
    parser.add_argument('-maxdepth', type=int, default=0, help='maximum read depth [default: 0 off]')
    parser.add_argument('-mincp', type=float, default=0.95, help='minimum combined probability: positive * (1-negative) [default: 0.95]')
    parser.add_argument('-minrp', type=float, default=0.75, help='minimum read proportion: reads seen that are unambiguously ref or alt [default: 0.75]')
    parser.add_argument('-seen', action='store_true', help='use the total number of reads seen at this position as the depth [instead of: ref + alt]')
    parser.add_argument('-loh', action='store_true', help='include mutations with loss of heterozygosity, labeled with LOH')
    args = parser.parse_args()

    pos = readFiles(args.p, args.seen, args.minrp) 
    pos_entry = compileEntries(pos, args.minlr, args.minaf, args.mindepth, args.maxdepth, False)
    neg_entry = {}
    pos_loh_entry = {}
    neg_loh_entry = {}
    if args.n:
        neg = readFiles(args.n, args.seen, args.minrp)
        neg_entry = compileEntries(neg, args.maxlr, args.maxaf, args.mindepth, args.maxdepth, True)
        if args.loh:
            (pos_loh_entry, neg_loh_entry) = compileLOH(pos, neg, args.minlr, args.maxlr, args.mindepth)

    outputResults(pos_entry, neg_entry, args.p, args.n, args.mincp)
    if args.loh:
        outputLOH(pos_loh_entry, neg_loh_entry, args.p, args.n)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt: 
        try:
            sys.exit(0)
        except SystemExit:
            os._exit(0)
