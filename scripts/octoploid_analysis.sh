#! /bin/bash
# Example workflow using EAGLE-RC for octoploid
# LAST at: http://last.cbrc.jp/ 
# STAR at: https://github.com/alexdobin/STAR
# featureCounts at: http://bioinf.wehi.edu.au/featureCounts/

## Homeolog identification

REF=genome
GTF=gtf_file
CPU=8

gffread -T -o refseq.gtf $GTF.gff3 # gff to gtf
gffread -g $REF.fa -w refseq.fa $GTF.gff3

# Assuming we have the genome assemblies and annotation for the A, B, D subgenomes
grep '^chr.A' refseq.gtf > refseq.chrA.gtf
grep '^chr.B' refseq.gtf > refseq.chrB.gtf
grep '^chr.C' refseq.gtf > refseq.chrC.gtf
grep '^chr.D' refseq.gtf > refseq.chrD.gtf

gffread -g $REF.fa -w chrA.exon.fa refseq.chrA.gtf
gffread -g $REF.fa -w chrB.exon.fa refseq.chrB.gtf
gffread -g $REF.fa -w chrC.exon.fa refseq.chrC.gtf
gffread -g $REF.fa -w chrD.exon.fa refseq.chrD.gtf

# **IMPORTANT**
# Any other transcriptomes should be in similar format
# i.e. chrC.exon.fa


# Reciprocal best hit, assuming you have exon sequences
lastdb -uNEAR -R01 chrA_db chrA.exon.fa
lastdb -uNEAR -R01 chrB_db chrB.exon.fa
lastdb -uNEAR -R01 chrC_db chrC.exon.fa
lastdb -uNEAR -R01 chrD_db chrD.exon.fa

lastal chrA_db -P$CPU -D10000000000 chrB.exon.fa | last-map-probs -m 0.49 > A.B.maf
lastal chrA_db -P$CPU -D10000000000 chrC.exon.fa | last-map-probs -m 0.49 > A.C.maf
lastal chrA_db -P$CPU -D10000000000 chrD.exon.fa | last-map-probs -m 0.49 > A.D.maf

lastal chrB_db -P$CPU -D10000000000 chrA.exon.fa | last-map-probs -m 0.49 > B.A.maf
lastal chrB_db -P$CPU -D10000000000 chrC.exon.fa | last-map-probs -m 0.49 > B.C.maf
lastal chrB_db -P$CPU -D10000000000 chrD.exon.fa | last-map-probs -m 0.49 > B.D.maf

lastal chrC_db -P$CPU -D10000000000 chrA.exon.fa | last-map-probs -m 0.49 > C.A.maf
lastal chrC_db -P$CPU -D10000000000 chrB.exon.fa | last-map-probs -m 0.49 > C.B.maf
lastal chrC_db -P$CPU -D10000000000 chrD.exon.fa | last-map-probs -m 0.49 > C.D.maf

lastal chrD_db -P$CPU -D10000000000 chrA.exon.fa | last-map-probs -m 0.49 > D.A.maf
lastal chrD_db -P$CPU -D10000000000 chrB.exon.fa | last-map-probs -m 0.49 > D.B.maf
lastal chrD_db -P$CPU -D10000000000 chrC.exon.fa | last-map-probs -m 0.49 > D.C.maf

# Create VCFs based on genotype differences between homeologs
python scripts/homeolog_genotypes.py -o A.vs.B -f exon -g refseq.gtf A.B.maf B.A.maf # coordinates based on A
python scripts/homeolog_genotypes.py -o A.vs.C -f exon -g refseq.gtf A.C.maf C.A.maf # coordinates based on A
python scripts/homeolog_genotypes.py -o A.vs.D -f exon -g refseq.gtf A.D.maf D.A.maf # coordinates based on A

python scripts/homeolog_genotypes.py -o B.vs.A -f exon -g refseq.gtf B.A.maf A.B.maf # coordinates based on B
python scripts/homeolog_genotypes.py -o B.vs.C -f exon -g refseq.gtf B.C.maf C.B.maf # coordinates based on B
python scripts/homeolog_genotypes.py -o B.vs.D -f exon -g refseq.gtf B.D.maf D.B.maf # coordinates based on B

python scripts/homeolog_genotypes.py -o C.vs.A -f exon -g refseq.gtf C.A.maf A.C.maf # coordinates based on C, will not have valid gtf.vcf as C not in refseq.gtf
python scripts/homeolog_genotypes.py -o C.vs.B -f exon -g refseq.gtf C.B.maf B.C.maf # coordinates based on C, will not have valid gtf.vcf as C not in refseq.gtf
python scripts/homeolog_genotypes.py -o C.vs.D -f exon -g refseq.gtf C.D.maf D.C.maf # coordinates based on C, will not have valid gtf.vcf as C not in refseq.gtf

python scripts/homeolog_genotypes.py -o D.vs.A -f exon -g refseq.gtf D.A.maf A.D.maf # coordinates based on D
python scripts/homeolog_genotypes.py -o D.vs.B -f exon -g refseq.gtf D.B.maf B.D.maf # coordinates based on D
python scripts/homeolog_genotypes.py -o D.vs.C -f exon -g refseq.gtf D.C.maf C.D.maf # coordinates based on D

# Quad copy homeologs
perl quad_homeolog.pl A.vs.B.reciprocal_best A.vs.C.reciprocal_best A.vs.D.reciprocal_best B.vs.C.reciprocal_best B.vs.D.reciprocal_best C.vs.D.reciprocal_best > homeolog.ABCD.list

# Triple copy homeologs
perl triple_homeolog.pl A.vs.B.reciprocal_best A.vs.C.reciprocal_best B.vs.C.reciprocal_best > homeolog.ABC.list
perl triple_homeolog.pl A.vs.B.reciprocal_best A.vs.D.reciprocal_best B.vs.D.reciprocal_best > homeolog.ABD.list
perl triple_homeolog.pl A.vs.C.reciprocal_best A.vs.D.reciprocal_best C.vs.D.reciprocal_best > homeolog.ACD.list
perl triple_homeolog.pl B.vs.C.reciprocal_best B.vs.D.reciprocal_best C.vs.D.reciprocal_best > homeolog.BCD.list


## Origin specific alignment with STAR for RNA-seq to Genome alignment
#GENDIR=/project/wheat/stargenome
#CHR="chrA chrB chrD"
#for n in $CHR; do
#    mkdir -p $GENDIR\_$n
#    STAR --runMode genomeGenerate --genomeDir $GENDIR\_$n --genomeFastaFiles $REF.$n.fa --sjdbGTFfile $GTF.gtf --runThreadN $CPU
#
#    for i in `ls *_R1.fastq.gz`; do 
#        F=`basename $i _R1.fastq.gz`
#        mkdir -p ./$n/star_$F
#        STAR --genomeDir $GENDIR\_$n --readFilesCommand zcat --readFilesIn $F\_R1.fastq.gz $F\_R2.fastq.gz \
#            --outFileNamePrefix star_$F- --runThreadN $CPU --genomeLoad NoSharedMemory \
#            --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
#            --outSJfilterCountUniqueMin 3 2 2 2 --outMultimapperOrder Random --outFilterType BySJout \
#            --outStd SAM | samtools view -Shb - > ./$n/$F.bam
#        samtools sort -o ./$n/$F.refsort.bam ./$n/$F.bam 
#        samtools index -c ./$n/$F.refsort.bam
#        mv star_$F-* ./$n/star_$F
#    done
#done

## Origin specific alignment with Bowtie2 for RNA-seq to Transcriptome alignment
GENDIR=/project/wheat/bowtie2
CHR="chrA chrB chrD chrC"
for n in $CHR; do
    mkdir -p $GENDIR
    bowtie2-build -f $CHR.exon.fa $GENDIR/$n --threads $CPU

    for i in `ls *_R1.fastq.gz`; do 
        F=`basename $i _R1.fastq.gz`
        mkdir -p ./$n/bowtie2_$F
        bowtie2 -x $GENDIR/$n -1 $F\_R1.fastq.gz -2 $F\_R2.fastq.gz --threads $CPU | samtools view -SbF4 - > ./$n/$F.bam
        samtools sort -o ./$n/$F.refsort.bam ./$n/$F.bam 
        samtools index -c ./$n/$F.refsort.bam
    done
done

## EAGLE-RC: read classification and the quantification with featureCounts
# Put the appropriate vcfs to the corresponding dir, i.e. A.vs.*.gtf.vcf in chrA
#CHR="chrA chrB chrD"
#for n in $CHR; do
#    cd $n
#    for i in `ls *.refsort.bam`; do 
#        F=`basename $i .refsort.bam`
#        for j in `ls *.gtf.vcf`; do
#            V=`basename $j .gtf.vcf`
#            eagle -t 8 -a $F.refsort.bam -r ../$REF.$n.fa -v $V.gtf.vcf --splice --rc 1> $F.$V.txt 2> $F.$V.readinfo.txt
#            eagle-rc --listonly -a $F.refsort.bam -o $F.$V -v $F.$V.txt $F.$V.readinfo.txt > $F.$V.list
#        done
#    done
#    cd ..
#done

# For transcriptome alignments
CHR="chrA chrB chrD chrC"
for n in $CHR; do
    cd $n
    for i in `ls *.refsort.bam`; do 
        F=`basename $i .refsort.bam`
        for j in `ls *.gtf.vcf`; do
            V=`basename $j .gtf.vcf`
            eagle -t 8 -a $F.refsort.bam -r ../$CHR.exon.fa -v $V.raw.vcf --rc 1> $F.$V.txt 2> $F.$V.readinfo.txt
            eagle-rc --listonly -a $F.refsort.bam -o $F.$V -v $F.$V.txt $F.$V.readinfo.txt > $F.$V.list
        done
    done
    cd ..
done


mkdir -p eagle
for i in `ls *_R1.fastq.gz`; do
    F=`basename $i _R1.fastq.gz`
    python scripts/ref4_consensus.py --pe -u -d -o eagle/$F.ref \
        -A chrA/$F.A.vs.B.list chrA/$F.A.vs.C.list chrA/$F.A.vs.D.list \
        -B chrB/$F.B.vs.A.list chrB/$F.B.vs.C.list chrB/$F.B.vs.D.list \
        -C chrC/$F.C.vs.A.list chrC/$F.C.vs.B.list chrC/$F.C.vs.D.list \
        -D chrD/$F.D.vs.A.list chrD/$F.D.vs.B.list chrD/$F.D.vs.C.list
    eagle-rc --refonly --readlist -a chrA/$F.refsort.bam -o eagle/$F.chrA eagle/$F.ref.chrA.list
    eagle-rc --refonly --readlist -a chrB/$F.refsort.bam -o eagle/$F.chrB eagle/$F.ref.chrB.list
    eagle-rc --refonly --readlist -a chrC/$F.refsort.bam -o eagle/$F.chrC eagle/$F.ref.chrC.list
    eagle-rc --refonly --readlist -a chrD/$F.refsort.bam -o eagle/$F.chrD eagle/$F.ref.chrD.list
#    featureCounts -T 8 -t exon -g transcript_id -a refseq.chrA.gtf -o eagle/$F.chrA.counts.txt eagle/$F.chrA.ref.bam
#    featureCounts -T 8 -t exon -g transcript_id -a refseq.chrB.gtf -o eagle/$F.chrB.counts.txt eagle/$F.chrB.ref.bam
#    featureCounts -T 8 -t exon -g transcript_id -a refseq.chrC.gtf -o eagle/$F.chrC.counts.txt eagle/$F.chrC.ref.bam
#    featureCounts -T 8 -t exon -g transcript_id -a refseq.chrD.gtf -o eagle/$F.chrD.counts.txt eagle/$F.chrD.ref.bam
done
