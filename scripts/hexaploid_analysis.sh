#! /bin/bash
# Example workflow using EAGLE-RC for hexaploid wheat
# gffread at: https://github.com/gpertea/gffread.  
# LAST at: http://last.cbrc.jp/ 
# STAR at: https://github.com/alexdobin/STAR
# featureCounts at: http://bioinf.wehi.edu.au/featureCounts/

## Homeolog identification

# The lyrata gene and scaffold ids should be modified to prepend Alyr_ or some other way make ids unique
REF=Taes_genome_2017_05
GTF=iwgsc_refseqv1.0_HighConf_UTR_2017May05
CPU=8

# Extract transcript sequences
gffread -T -o refseq.gtf $GTF.gff3 # gff to gtf
gffread -g $REF.fa -w refseq.fa $GTF.gff3

grep '^chr.A' refseq.gtf > refseq.chrA.gtf
grep '^chr.B' refseq.gtf > refseq.chrB.gtf
grep '^chr.D' refseq.gtf > refseq.chrD.gtf

gffread -g $REF.fa -w chrA.exon.fa refseq.chrA.gtf
gffread -g $REF.fa -w chrB.exon.fa refseq.chrB.gtf
gffread -g $REF.fa -w chrD.exon.fa refseq.chrD.gtf

# Reciprocal best hit
lastdb -uNEAR -R01 chrA_db chrA.exon.fa
lastdb -uNEAR -R01 chrB_db chrB.exon.fa
lastdb -uNEAR -R01 chrD_db chrD.exon.fa

lastal chrA_db -P$CPU -D10000000000 chrB.exon.fa | last-map-probs -m 0.49 > A.B.maf
lastal chrB_db -P$CPU -D10000000000 chrA.exon.fa | last-map-probs -m 0.49 > B.A.maf

lastal chrB_db -P$CPU -D10000000000 chrD.exon.fa | last-map-probs -m 0.49 > B.D.maf
lastal chrD_db -P$CPU -D10000000000 chrB.exon.fa | last-map-probs -m 0.49 > D.B.maf

lastal chrA_db -P$CPU -D10000000000 chrD.exon.fa | last-map-probs -m 0.49 > A.D.maf
lastal chrD_db -P$CPU -D10000000000 chrA.exon.fa | last-map-probs -m 0.49 > D.A.maf

# Create VCFs based on genotype differences between homeologs
python scripts/homeolog_genotypes.py -o A.vs.B -f exon -g refseq.gtf A.B.maf B.A.maf # coordinates based on A
python scripts/homeolog_genotypes.py -o B.vs.A -f exon -g refseq.gtf B.A.maf A.B.maf # coordinates based on B

python scripts/homeolog_genotypes.py -o B.vs.D -f exon -g refseq.gtf B.D.maf D.B.maf # coordinates based on B
python scripts/homeolog_genotypes.py -o D.vs.B -f exon -g refseq.gtf D.B.maf B.D.maf # coordinates based on D

python scripts/homeolog_genotypes.py -o A.vs.D -f exon -g refseq.gtf A.D.maf D.A.maf # coordinates based on A
python scripts/homeolog_genotypes.py -o D.vs.A -f exon -g refseq.gtf D.A.maf A.D.maf # coordinates based on D

# Triple copy homeologs
perl triple_homeolog.pl A.vs.B.reciprocal_best B.vs.D.reciprocal_best A.vs.D.reciprocal_best > homeolog.ABD.list

# Subgenome unique transcripts
cat A.vs.B.reciprocal_best A.vs.D.reciprocal_best | cut -f1 | sort | uniq > A.vs.all.list
cat B.vs.A.reciprocal_best B.vs.D.reciprocal_best | cut -f1 | sort | uniq > B.vs.all.list
cat D.vs.A.reciprocal_best D.vs.B.reciprocal_best | cut -f1 | sort | uniq > D.vs.all.list
grep $'mRNA\t' $GTF.gff3 | grep $'chr.A\t' | perl -ne 'chomp; m/ID=(.*?);/; print "$1\n";' > chrA.exon.list
grep $'mRNA\t' $GTF.gff3 | grep $'chr.B\t' | perl -ne 'chomp; m/ID=(.*?);/; print "$1\n";' > chrB.exon.list
grep $'mRNA\t' $GTF.gff3 | grep $'chr.D\t' | perl -ne 'chomp; m/ID=(.*?);/; print "$1\n";' > chrD.exon.list
python scripts/tablize.py -v0 A.vs.all.list chrA.exon.list > chrA.only.list
python scripts/tablize.py -v0 B.vs.all.list chrB.exon.list > chrB.only.list
python scripts/tablize.py -v0 D.vs.all.list chrD.exon.list > chrD.only.list

## Origin specific alignment with STAR
GENDIR=/project/wheat/stargenome
CHR="chrA chrB chrD"
for n in $CHR; do
    mkdir -p $GENDIR\_$n
    STAR --runMode genomeGenerate --genomeDir $GENDIR\_$n --genomeFastaFiles $REF.$n.fa --sjdbGTFfile $GTF.gtf --runThreadN $CPU

    for i in `ls *_R1.fastq.gz`; do 
        F=`basename $i _R1.fastq.gz`
        mkdir -p ./$n/star_$F
        STAR --genomeDir $GENDIR\_$n --readFilesCommand zcat --readFilesIn $F\_R1.fastq.gz $F\_R2.fastq.gz \
            --outFileNamePrefix star_$F- --runThreadN $CPU --genomeLoad NoSharedMemory \
            --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
            --outSJfilterCountUniqueMin 3 2 2 2 --outMultimapperOrder Random --outFilterType BySJout \
            --outStd SAM | samtools view -Shb - > ./$n/$F.bam
        samtools sort -o ./$n/$F.refsort.bam ./$n/$F.bam 
        samtools index -c ./$n/$F.refsort.bam
        mv star_$F-* ./$n/star_$F
    done
done

## EAGLE-RC: read classification and the quantification with featureCounts
# Put the appropriate vcfs to the corresponding dir, i.e. A.vs.*.gtf.vcf in chrA
for n in $CHR; do
    cd $n
    for i in `ls *.refsort.bam`; do 
        F=`basename $i .refsort.bam`
        for j in `ls *.gtf.vcf`; do
            V=`basename $j .gtf.vcf`
            eagle -t 8 -a $F.refsort.bam -r ../$REF.$n.fa -v $V.gtf.vcf --splice --rc 1> $F.$V.txt 2> $F.$V.readinfo.txt
            eagle-rc --listonly -a $F.refsort.bam -o $F.$V -v $F.$V.txt $F.$V.readinfo.txt > $F.$V.list
        done
    done
    cd ..
done

mkdir -p eagle
for i in `ls *_R1.fastq.gz`; do
    F=`basename $i _R1.fastq.gz`
    python scripts/ref3_consensus.py --pe -u -d -o eagle/$F.ref \
        -A chrA/$F.A.vs.B.list chrA/$F.A.vs.D.list \
        -B chrB/$F.B.vs.A.list chrB/$F.B.vs.D.list \
        -D chrD/$F.D.vs.A.list chrD/$F.D.vs.B.list
    eagle-rc --refonly --readlist -a chrA/$F.refsort.bam -o eagle/$F.chrA eagle/$F.ref.chrA.list
    eagle-rc --refonly --readlist -a chrB/$F.refsort.bam -o eagle/$F.chrB eagle/$F.ref.chrB.list
    eagle-rc --refonly --readlist -a chrD/$F.refsort.bam -o eagle/$F.chrD eagle/$F.ref.chrD.list
    featureCounts -T 8 -t exon -g transcript_id -a refseq.chrA.gtf -o eagle/$F.chrA.counts.txt eagle/$F.chrA.ref.bam
    featureCounts -T 8 -t exon -g transcript_id -a refseq.chrB.gtf -o eagle/$F.chrB.counts.txt eagle/$F.chrB.ref.bam
    featureCounts -T 8 -t exon -g transcript_id -a refseq.chrD.gtf -o eagle/$F.chrD.counts.txt eagle/$F.chrD.ref.bam
done

# Double and triple homeolog counts
python scripts/tablize.py -skip 1 -a -i 0 -c 6 eagle/*.chrA.counts.txt > eagle.chrA.tsv
python scripts/tablize.py -skip 1 -a -i 0 -c 6 eagle/*.chrB.counts.txt > eagle.chrB.tsv
python scripts/tablize.py -skip 1 -a -i 0 -c 6 eagle/*.chrD.counts.txt > eagle.chrD.tsv

# Triple homeologs counts in terms of chrA gene id, transcript level
cut -f 1 homeolog.ABD.list > homeolog.A.list
python scripts/tablize.py -a homeolog.A.list eagle.chrA.tsv | sort -k1 > eagle.chrA.homeolog.tsv
awk '{print $2"\t"$1;}' homeolog.ABD.list > homeolog.B.list
python scripts/tablize.py -a homeolog.B.list eagle.chrB.tsv | cut -f 2,3- | sort -k1 > eagle.chrB.homeolog.tsv
awk '{print $3"\t"$1;}' homeolog.ABD.list > homeolog.D.list
python scripts/tablize.py -a homeolog.D.list eagle.chrD.tsv | cut -f 2,3- | sort -k1 > eagle.chrD.homeolog.tsv

# Homeolog counts in terms of halleri gene id, gene level by summation of transcript counts
perl -ne 'chomp; @t=split(/\s+/); @i=split(/\./, $t[0]); $a=join("\t", @t[1..$#t]); print "$i[0]\t$a\n";' eagle.chrA.homeolog.tsv > temp.txt
python scripts/tablize.py -add temp.txt > eagle.chrA.homeolog.genelevel.tsv
perl -ne 'chomp; @t=split(/\s+/); @i=split(/\./, $t[0]); $a=join("\t", @t[1..$#t]); print "$i[0]\t$a\n";' eagle.chrB.homeolog.tsv > temp.txt
python scripts/tablize.py -add temp.txt > eagle.chrB.homeolog.genelevel.tsv
perl -ne 'chomp; @t=split(/\s+/); @i=split(/\./, $t[0]); $a=join("\t", @t[1..$#t]); print "$i[0]\t$a\n";' eagle.chrD.homeolog.tsv > temp.txt
python scripts/tablize.py -add temp.txt > eagle.chrD.homeolog.genelevel.tsv

# Subgenome unique mapped reads
for i in `ls *_R1.fastq.gz`; do
    F=`basename $i _R1.fastq.gz`
    echo "" > dummy.txt
    eagle-rc --refonly --readlist -a chrA/$F.refsort.bam -u chrB/$F.refsort.bam,chrD/$F.refsort.bam -o eagle/$F.chrA.only dummy.txt
    eagle-rc --refonly --readlist -a chrB/$F.refsort.bam -u chrA/$F.refsort.bam,chrD/$F.refsort.bam -o eagle/$F.chrB.only dummy.txt
    eagle-rc --refonly --readlist -a chrD/$F.refsort.bam -u chrA/$F.refsort.bam,chrB/$F.refsort.bam -o eagle/$F.chrD.only dummy.txt
    featureCounts -T 8 -t exon -g transcript_id -a refseq.chrA.gtf -o eagle/$F.chrA.only.counts.txt eagle/$F.chrA.only.ref.bam
    featureCounts -T 8 -t exon -g transcript_id -a refseq.chrB.gtf -o eagle/$F.chrB.only.counts.txt eagle/$F.chrB.only.ref.bam
    featureCounts -T 8 -t exon -g transcript_id -a refseq.chrD.gtf -o eagle/$F.chrD.only.counts.txt eagle/$F.chrD.only.ref.bam
done

python scripts/tablize.py -skip 1 -a -i 0 -c 6 eagle/*.chrA.only.counts.txt > eagle.chrA.only.tsv
python scripts/tablize.py -skip 1 -a -i 0 -c 6 eagle/*.chrB.only.counts.txt > eagle.chrB.only.tsv
python scripts/tablize.py -skip 1 -a -i 0 -c 6 eagle/*.chrD.only.counts.txt > eagle.chrD.only.tsv

# Subgenome unique mapped reads in subgenome unique genes (i.e. non-homeologs)
python scripts/tablize.py -a chrA.only.list eagle.chrA.tsv > eagle.chrA.only.unique.tsv
python scripts/tablize.py -a chrB.only.list eagle.chrB.tsv > eagle.chrB.only.unique.tsv
python scripts/tablize.py -a chrD.only.list eagle.chrD.tsv > eagle.chrD.only.unique.tsv
