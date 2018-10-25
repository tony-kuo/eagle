#!/usr/bin/python

# Copyright 2016 Tony Kuo
# This program is distributed under the terms of the GNU General Public License

# Utility Script
# Compare results of EAGLE to VCF truth data
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

def readEAGLE(filename):
    entry = {};
    with open(filename, 'r') as fh:
        for line in fh:
            if re.match('^#', line): continue;
            t = line.strip().split('\t');
            entry[ '\t'.join(t[0:4]) ] = t;
    fh.close;
    return(entry);

def readVCF(filename, eagle, raw):
    entry = {};
    with open(filename, 'r') as fh:
        for line in fh:
            if re.match('^#', line): 
                print(line.strip());
                continue;

            t = line.strip().split('\t');
            key = "{}\t{}\t{}\t{}".format(t[0], t[1], t[3], t[4]);
            if key in eagle:
                if raw:
                    af = float(eagle[key][6]) / float(eagle[key][4]);
                    lr = float(eagle[key][8]);
                else:
                    af = float(eagle[key][5]);
                    lr = float(eagle[key][6]);
                t[7] += ";EAGLEAF={:.2f};EAGLELR={:.4f}".format(af, lr);
                print('\t'.join(t))

    fh.close;
    return(entry);

def main():
    parser = argparse.ArgumentParser(description='Combine VCF and EAGLE output');
    parser.add_argument('-v', help='VCF file');
    parser.add_argument('-e', help='EAGLE file');
    parser.add_argument('-raw', action='store_true', help='raw EAGLE output, otherwise output after compile_likelihoods');
    args = parser.parse_args();

    eagle = readEAGLE(args.e); 
    readVCF(args.v, eagle, args.raw); 

if __name__ == '__main__':
    try:
        main();
    except KeyboardInterrupt: 
        try:
            sys.exit(0);
        except SystemExit:
            os._exit(0);
