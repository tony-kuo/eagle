#!/usr/bin/python

# Copyright 2016 Tony Kuo
# This program is distributed under the terms of the GNU General Public License

# Utility Script
# Finds reciprocal best hit mRNA transcripts from LAST MAF output
# Determines sequences differences as putative candidate variants / genotypes
# Writes genotypes into VCFs in the RNA and DNA coordinates of the "database" reference genome

# Example:
## LAST on transcripts, database is the REF's transcripts
# lastdb -uNEAR -R01 database database.fa
# lastal -D10000000000 database -P8 query.fa | last-map-probs -m 0.49 > d1.maf
# lastdb -uNEAR -R01 query query.fa
# lastal -D10000000000 query -P8 database.fa | last-map-probs -m 0.49 > d2.maf

#python ~/homeolog_genotypes.py -o output -g $GTF d1.maf d2.maf # coordinates based on database.fa
#cut -f 1,6 output.rna.vcf | sort -k1 | uniq > database.query.txt

from __future__ import print_function;
import argparse;
import sys;
from re import split, match;
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 

def naturalSort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower();
    alphanum_key = lambda key: [convert(c) for c in split('([0-9]+)', key)];
    return(sorted(l, key=alphanum_key));

complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'Y':'R', 'R':'Y', 'W':'W', 'S':'S', 'K':'M', 'M':'K', 'D':'H', 'H':'D', 'V':'B', 'B':'V', 'N':'N', '-':'-'};

def readGTF(fn):
    entry = {};
    with open(fn, 'r') as fh:
        for line in fh:
            if match('^#', line): continue;

            t = line.strip().split('\t');
            if t[2] != "exon": continue;

            m = match('.*transcript_id "(\S+)";', t[8]);
            if not m: continue;

            if m.group(1) not in entry:
                entry[m.group(1)] = {};
                entry[m.group(1)]['id'] = t[0];
                entry[m.group(1)]['sense'] = t[6];
                entry[m.group(1)]['offset'] = [0];
                entry[m.group(1)]['exon'] = [(int(t[3]), int(t[4]))];
            else:
                entry[m.group(1)]['offset'].append(int(t[3]) - entry[m.group(1)]['exon'][-1][-1] - 1);
                entry[m.group(1)]['exon'].append((int(t[3]), int(t[4])));
    return(entry);

def readMAF(fn, reciprocal):
    entry = {};
    with open(fn, 'r') as fh:
        rSeq = qSeq = "";
        for line in fh:
            if match('^#', line): continue;

            t = line.strip().split();
            if line[0] == "a":
                a = line;
            elif line[0] == "s":
                if not rSeq:
                    rId = t[1];
                    rStart = int(t[2]) + 1; # switch to 1-index
                    rEnd = rStart + int(t[3]) - 1;
                    rLength = int(t[5]);
                    rSeq = t[6].upper();
                elif not qSeq:
                    qId = t[1];
                    qStart = int(t[2]) + 1; # switch to 1-index
                    qEnd = qStart + int(t[3]) - 1;
                    qLength = int(t[5]);
                    qSeq = t[6].upper();
                    sense = t[4];
            elif not t:
                if reciprocal: # flip reference and query as the ID
                    ID = rId + " " + qId;
                else:
                    ID = qId + " " + rId;

                entry[ID] = (a, rId, rStart, rEnd, rLength, rSeq, qId, qStart, qEnd, qLength, qSeq, sense);
                rSeq = qSeq = "";
    return(entry);

def reciprocalBestHit(e1, e2, prefix, gtf):
    entry = {};
    dna_entry = {};

    gen = (x for x in e1 if x in e2);
    for ID in gen:
        # Find variants and coordinates, always from the coordinate system of e1 reference (which is always + sense)
        a = e1[ID][5];
        b = e1[ID][10];
        aa = []
        bb = [];
        offset_a = 0;
        for i in range(0, len(a)):
            if a[i] == "-": offset_a += 1;
            if a[i] != b[i]:
                aa.append(a[i]);
                bb.append(b[i]);
            else: 
                if aa and bb:
                    rna_pos = e1[ID][2] + i - len(aa) - offset_a;

                    if aa[0] == '-': aa = ['-'];
                    elif bb[0] == '-': bb = ['-'];

                    entry["{}\t{}\t.\t{}\t{}\t{}".format(e1[ID][1], rna_pos, "".join(aa), "".join(bb), e2[ID][1])] = 0;

                    rID = e1[ID][1];
                    if rID in gtf:
                        if gtf[rID]['sense'] == "-": # reverse complement base if negative sense
                            aa = list(reversed([complement[x] for x in aa]));
                            bb = list(reversed([complement[x] for x in bb]));
                            rna_pos = e1[ID][4] - rna_pos + 1; # find position from back instead, since back is the front

                        dna_pos = gtf[rID]['exon'][0][0] + rna_pos - 1; # reference from start of first exon
                        for i in range(0, len(gtf[rID]['exon'])):
                            dna_pos += gtf[rID]['offset'][i];
                            if dna_pos <= gtf[rID]['exon'][i][-1]: break;
                        dna_entry["{}\t{}\t.\t{}\t{}".format(gtf[rID]['id'], dna_pos, "".join(aa), "".join(bb))] = "{}\t{}".format(e1[ID][1], e2[ID][1]);


                aa = [];
                bb = [];

    with open(prefix+".rna.vcf", 'w') as fh:
        for i in naturalSort(entry):
            print(i, file=fh);
    if gtf:
        with open(prefix+".dna.vcf", 'w') as fh:
            for i in naturalSort(dna_entry):
                print("{}\t{}".format(i, dna_entry[i]), file=fh);

def main():
    parser = argparse.ArgumentParser(description='Find homeologs from reciprocal best hit');
    parser.add_argument('files', nargs='+', help='2 MAF alignment files');
    parser.add_argument('-o', type=str, help='output file prefix');
    parser.add_argument('-g', type=str, default="", help='gtf file to translate into reference genome coordinates');
    if len(sys.argv) == 1:
        parser.print_help();
        sys.exit(1);
    args = parser.parse_args();

    gtf = {};
    if args.g:
        gtf = readGTF(args.g);
        
    e1 = readMAF(args.files[0], True);
    e2 = readMAF(args.files[1], False);
    reciprocalBestHit(e1, e2, args.o, gtf);

if __name__ == '__main__':
    main()

