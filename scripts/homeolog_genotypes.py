#!/usr/bin/python

# Copyright 2016 Tony Kuo
# This program is distributed under the terms of the GNU General Public License

# Utility Script
# Finds reciprocal best hit transcripts from LAST (http://last.cbrc.jp/) MAF output
# Determines sequences differences as putative candidate variants / genotypes
# Writes genotypes into VCFs in the RNA and DNA coordinates of the "database" reference genome

## If you need transcripts and have the annotation GFF, using gffread (Cufflinks: https://github.com/cole-trapnell-lab/cufflinks.git):
#    gffread -MKQ -w transcripts.fa -g REF.fa GTF.gff

## LAST reciprocal best hit to find homeologs from transcripts
# S1 as reference:
#   lastdb -uNEAR -R01 S1 S1.fa
#   lastal -D10000000000 S1 -P8 S2.fa | last-map-probs -m 0.49 > d1.maf
# S2 as reference:
#   lastdb -uNEAR -R01 S2 S2.fa
#   lastal -D10000000000 S2 -P8 S1.fa | last-map-probs -m 0.49 > d2.maf

## Given GFF annotation:
#   gffread S1.gff -T -o S1.gtf
#   python ~/homeolog_genotypes.py -o output -f exon -g S1.gtf d1.maf d2.maf # coordinates based on S1.fa

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

def readGTF(fn, field):
    entry = {};
    with open(fn, 'r') as fh:
        for line in fh:
            if match('^#', line): continue;

            t = line.strip().split('\t');
            if t[2] != field: continue;

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

def readMAF(fn, refFirst):
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
                    rStart = int(t[2]);
                    rAlnLen = int(t[3]);
                    rLength = int(t[5]);
                    rSeq = t[6].upper();
                elif not qSeq:
                    qId = t[1];
                    qStart = int(t[2]);
                    qAlnLen = int(t[3]);
                    qLength = int(t[5]);
                    qSeq = t[6].upper();
                    sense = t[4];
            elif not t:
                if refFirst: # flip reference and query as the ID
                    ID = rId + "\t" + qId;
                else:
                    ID = qId + "\t" + rId;

                if rAlnLen > 200 and qAlnLen > 200:
                    entry[ID] = (a, rId, rStart, rAlnLen, rLength, rSeq, qId, qStart, qAlnLen, qLength, qSeq, sense);
                rSeq = qSeq = "";
    return(entry);

def reciprocalBestHit(e1, e2, prefix, gtf):
    entry = {};
    dna_entry = {};

    gen = (x for x in e1 if x in e2);
    with open(prefix+".reciprocal_best", 'w') as fh:
        for ID in gen:
            # Find variants and coordinates, always from the coordinate system of e1 transcriptome (which is always + sense)
            print(ID, file=fh);

            a = e1[ID][5];
            b = e1[ID][10];
            aa = []
            bb = [];
            offset_a = 0;
            for i in range(0, len(a)):
                if a[i] != b[i]:
                    if a[i] == "-": offset_a += 1;
                    aa.append(a[i]);
                    bb.append(b[i]);
                else: 
                    if aa and bb:
                        rna_pos = e1[ID][2] + i - offset_a - (len(aa) - aa.count('-')) + 1;

                        if aa[0] == '-' and aa.count('-') == len(aa): 
                            aa = ['-'];
                        elif bb[0] == '-' and bb.count('-') == len(bb): 
                            bb = ['-'];

                        if len(aa) > 1:
                            aa = [x for x in aa if x != '-'];
                        if len(bb) > 1:
                            bb = [x for x in bb if x != '-'];

                        rProportion = float(e1[ID][3]) / e1[ID][4];
                        qProportion = float(e1[ID][8]) / e1[ID][9];
                        entry["{}\t{}\t.\t{}\t{}".format(e1[ID][1], rna_pos, "".join(aa), "".join(bb))] = "{};rProp={:.4f};qProp={:.4f}".format(e2[ID][1], rProportion, qProportion);

                        rID = e1[ID][1];
                        if rID in gtf:
                            dna_pos = gtf[rID]['exon'][0][0] + rna_pos - 1; # reference from start of first exon
                            if gtf[rID]['sense'] == "-": # reverse complement base if negative sense transcript
                                aa = list(reversed([complement[x] for x in aa]));
                                bb = list(reversed([complement[x] for x in bb]));
                                rna_pos = e1[ID][4] - rna_pos + 1; # find position from back instead, since back is the front
                                dna_pos = gtf[rID]['exon'][0][0] + rna_pos - 1; # reference from start of first exon
                                dna_pos = dna_pos - len(aa) + 1; # compensate for 3'-5' to 5'-3' switch

                            for i in range(0, len(gtf[rID]['exon'])):
                                dna_pos += gtf[rID]['offset'][i];
                                if dna_pos <= gtf[rID]['exon'][i][-1]: break;

                            dna_entry["{}\t{}\t.\t{}\t{}".format(gtf[rID]['id'], dna_pos, "".join(aa), "".join(bb))] = "{};{};rProp={:.4f};qProp={:.4f}".format(e1[ID][1], e2[ID][1], rProportion, qProportion);


                    aa = [];
                    bb = [];

    with open(prefix+".raw.vcf", 'w') as fh:
        for i in naturalSort(entry):
            print("{}\t{}".format(i, entry[i]), file=fh);
    if gtf:
        with open(prefix+".gtf.vcf", 'w') as fh:
            for i in naturalSort(dna_entry):
                print("{}\t{}".format(i, dna_entry[i]), file=fh);

def main():
    parser = argparse.ArgumentParser(description='Find homeologs from reciprocal best hit. Reference coordinates based on the reference in first MAF file.');
    parser.add_argument('files', nargs='+', help='MAF alignment files');
    parser.add_argument('-o', type=str, help='output file prefix');
    parser.add_argument('-g', type=str, default="", help='gtf file to translate into reference genome coordinates, using (default CDS) entries.  Check the fasta for consistency.');
    parser.add_argument('-f', type=str, default="CDS", help='use another field for coding entries, i.e. exon');
    if len(sys.argv) == 1:
        parser.print_help();
        sys.exit(1);
    args = parser.parse_args();

    gtf = {};
    if args.g:
        gtf = readGTF(args.g, args.f);
        
    e1 = readMAF(args.files[0], True);
    e2 = readMAF(args.files[1], False);
    reciprocalBestHit(e1, e2, args.o, gtf);

if __name__ == '__main__':
    main()

