#!/usr/bin/python

# Copyright 2016 Tony Kuo
# This program is distributed under the terms of the GNU General Public License

# Utility Script
# Works with tab-delimited tables, where samples are in rows and variables are in columns
# Useful for finding the union or intersection of multiple table files, or extracting column values

from __future__ import print_function
import argparse
import re
import os, sys

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 

def naturalSort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)] 
    if sortdesc: return sorted(l, key=alphanum_key, reverse=True)
    return sorted(l, key=alphanum_key)

def sortbyMean(l, count):
    if sortdesc: meanvalues = [( k, sum([float(v[0]) for v in l[k].values()]) / count ) for k in l]
    else: meanvalues = [( k, sum([1/float(v[0]) for v in l[k].values()]) / count ) for k in l]
    return([a[0] for a in sorted(meanvalues, key=lambda tup:tup[1], reverse=True)])

def readFiles(files, idcol, valcol, delim, addtogether):
    entry = {}
    numval = {}
    if len(idcol) > 0:
        c = idcol.strip().split(',')
        idcol = []
        for i in c:
            i = i.split('-')
            if len(i) > 1: idcol.extend(range(int(i[0]), int(i[1])+1))
            else: idcol.append(int(i[0]))
        numid = len(idcol)-1
    else: numid = 1

    if len(valcol) > 0:
        c = valcol.strip().split(',')
        valcol = []
        for i in c:
            i = i.split('-')
            if len(i) > 1: valcol.extend(range(int(i[0]), int(i[1])+1))
            else: valcol.append(int(i[0]))

    for fn in files: 
        with open(fn, 'r') as fh:
            header_count = 0
            numlines = 0
            for line in fh:
                if re.match('^#', line): continue
                if header_count < skipheaders: # skip header lines
                    header_count += 1
                    continue

                numlines += 1
                if delim == 's+': t = line.strip().split()
                else: t = line.strip().split(delim)
                if len(idcol) > 1: key = '\t'.join([t[i] for i in idcol])
                else: key = t[idcol[0]]

                if len(valcol) != 0: 
                    values = [t[i] for i in valcol] # user provided value columns
                    numval[fn] = len(valcol)
                else: 
                    values = t[(idcol[-1]+1):] # start from last idcol + 1 for value columns
                    numval[fn] = len(t[(idcol[-1]+1):])-1

                if key not in entry:
                    entry[key] = {}
                    entry[key][fn] = values
                elif key in entry and fn in entry[key] and addtogether:
                    entry[key][fn] = [float(values[i]) + float(entry[key][fn][i]) for i in range(0, len(values))]
                else:
                    entry[key][fn] = values
        fh.close
        print("Read:\t{0}\t{1} entries".format(fn, numlines), file=sys.stderr)
    return(entry, numid, numval)

def writeTable(entry, numid, numval, files):
    if header:
        line = 'ID' + '\t'*numid
        for fn in files: 
            if fn in numval and numval[fn]: line += fn + '\t'*numval[fn]
        print(line.strip())

    if sortbymean: sortedkeys = sortbyMean(entry, len(numval.values()))
    else: sortedkeys = naturalSort(entry)

    f = len(files)
    for k1 in sortedkeys: #entry.items():
        #k1 = i[0] # entry id
        skipentry = False
        n = len(list(entry[k1].keys()))
        if notexistfirst: # check if entry does not exist in first files
            skipentry = True
            if files[0] not in entry[k1]: skipentry = False
        elif existone:
            if n > 1: skipentry = True
        elif notexistall: # check if entry does not exist in all files
            if n == f: skipentry = True
        elif existall: # check if entry exists in all files
            if n != f: skipentry = True
        elif existfirst: # check if entry exists in first file
            skipentry = True
            if files[0] in entry[k1]: skipentry = False

        if skipentry: continue

        line = k1 + '\t'
        for fn in files:
            if fn in entry[k1]: line += '\t'.join(map(str,entry[k1][fn]))+'\t'
            else: 
                if fn not in numval or (notexistfirst and fn == files[0]): continue
                line += (missing+'\t')*(numval[fn]+1)
        print(line.strip())

existall = False
existfirst = False
notexistall = False
notexistfirst = False
existone = False
skipheaders = 0
header = False
missing = ''
sortbymean = False
sortdesc = False
def main():
    parser = argparse.ArgumentParser(description='Data files should be tab-delimited column files: [id] [value1] [value2] ...')
    parser.add_argument('files', nargs='+', help='data file')
    parser.add_argument('-i', type=str, default='0', help='comma separated values to specify the id columns, 0-based indexing (default: 0)')
    parser.add_argument('-c', type=str, default='', help='comma separated values to specify the value columns to extract from each file, 0-based indexing (default: 1:end)')
    parser.add_argument('-a', action='store_true', help='only print entry if *all* data files contain said entry')
    parser.add_argument('-a0', action='store_true', help='print entries in *1st data file* and any others that contain said entry')
    parser.add_argument('-v', action='store_true', help='only print entry if *not all* data files contain said entry, over-rides [-a]')
    parser.add_argument('-v0', action='store_true', help='only print entry if *1st data file does not* contain said entry, over-rides [-a & -v]')
    parser.add_argument('-v1', action='store_true', help='only print entry if *one file* contains said entry, over-rides [-a & -v]')
    parser.add_argument('-skip', type=int, default=0, help='skip a number of header rows for all files (default: 0)')
    parser.add_argument('-header', action='store_true', help='print header')
    parser.add_argument('-delim', type=str, default='\t', help='delimiter for columns (default: \t, s+ for whitespace)')
    parser.add_argument('-miss', type=str, default='', help='symbol for missing values (default: '')')
    parser.add_argument('-mean', action='store_true', help='sort by mean of values (numerical data only), ascending order, rather than id')
    parser.add_argument('-desc', action='store_true', help='sort by descending order')
    parser.add_argument('-add', action='store_true', help='if duplicate id in the same file then add the values to aggregate (numerical data only)')
    args = parser.parse_args()

    global existall, existfirst, notexistall, existone, notexistfirst, header, missing, skipheaders, sortbymean, sortdesc
    existall = args.a
    existfirst = args.a0
    notexistall = args.v
    existone = args.v1
    notexistfirst = args.v0
    header = args.header
    missing = args.miss
    skipheaders = args.skip
    sortbymean = args.mean
    sortdesc = args.desc

    (entry, numid, numval) = readFiles(args.files, args.i, args.c, args.delim, args.add)
    #if addtogetherfiles:
    #    for key in entry:
    #        a = sum([float(x) for x in entry[key] for x in entry[key][x]])
    #        entry[key][files[0]] = [str(a)]
    #    args.files = [args.files[0]]
    writeTable(entry, numid, numval, args.files)

if __name__ == '__main__':
    main()
    #os._exit(1)
