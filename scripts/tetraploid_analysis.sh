#! /bin/bash
# Example workflow using EAGLE-RC for tetraploid Arabidopsis kamchatica (Arabidopsis halleri + Arabidopsis lyrata)
# gffread at: https://github.com/gpertea/gffread.  
# LAST at: http://last.cbrc.jp/ 
# STAR at: https://github.com/alexdobin/STAR
# featureCounts at: http://bioinf.wehi.edu.au/featureCounts/

## Homeolog identification

# The lyrata gene and scaffold ids should be modified to prepend Alyr_ or some other way make ids unique
HALREF=Ahal_v2_2
HALGTF=Ahal_v2_2
LYRREF=Alyr_v2_2
LYRGTF=Alyr_v2_2_1
CPU=8

# Extract transcript sequences
gffread -T -o $HALGTF.gtf $HALGTF.gff # gff to gtf
gffread -g $HALREF.fa -w halrefseq.fa $HALGTF.gff
gffread -T -o $LYRGTF.gtf $LYRGTF.gff # gff to gtf
gffread -g $LYRREF.fa -w lyrrefseq.fa $LYRGTF.gff

# Reciprocal best hit
lastdb -uNEAR -R01 hal_db halrefseq.fa
lastdb -uNEAR -R01 lyr_db lyrrefseq.fa

lastal hal_db -P$CPU -D10000000000 lyrrefseq.fa | last-map-probs -m 0.49 > H.L.maf
lastal lyr_db -P$CPU -D10000000000 halrefseq.fa | last-map-probs -m 0.49 > L.H.maf

# Create VCFs based on genotype differences between homeologs
python scripts/homeolog_genotypes.py -o H.vs.L -f exon -g $HALGTF.gtf H.L.maf L.H.maf # coordinates based on hal
python scripts/homeolog_genotypes.py -o L.vs.H -f exon -g $LYRGTF.gtf L.H.maf H.L.maf # coordinates based on lyr

# Subgenome unique transcripts
cat $HALGTF.gtf | perl -ne 'chomp; m/transcript_id "(.*?)";/; print "$1\n";' | sort | uniq > H.all.list
cat $LYRGTF.gtf | perl -ne 'chomp; m/transcript_id "(.*?)";/; print "$1\n";' | sort | uniq > L.all.list
python scripts/tablize.py -v0 H.vs.L.reciprocal_best H.all.list > H.only.list
python scripts/tablize.py -v0 L.vs.H.reciprocal_best L.all.list > L.only.list

## Origin specific alignment with STAR

HALGENDIR=/project/hal/stargenome
mkdir -p $HALGENDIR
STAR --runMode genomeGenerate --genomeDir $HALGENDIR --genomeFastaFiles $HALREF.fa --sjdbGTFfile $HALGTF.gtf --runThreadN $CPU 

LYRGENDIR=/project/lyr/stargenome
mkdir -p $LYRGENDIR
STAR --runMode genomeGenerate --genomeDir $LYRGENDIR --genomeFastaFiles $LYRREF.fa --sjdbGTFfile $LYRGTF.gtf --runThreadN $CPU

for i in `ls *_R1.fastq.gz`; do 
    F=`basename $i _R1.fastq.gz`
    mkdir -p ./hal/star_$F
    STAR --genomeDir $HALGENDIR --readFilesCommand zcat --readFilesIn $F\_R1.fastq.gz $F\_R2.fastq.gz \
        --outFileNamePrefix star_$F- --runThreadN $CPU --genomeLoad NoSharedMemory \
        --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
        --outSJfilterCountUniqueMin 3 2 2 2 --outMultimapperOrder Random \
        --outFilterType BySJout --outStd SAM | samtools view -Shb - > ./hal/$F.bam
    samtools sort -o ./hal/$F.refsort.bam ./hal/$F.bam 
    samtools index -c ./hal/$F.refsort.bam
    mv star_$F-* ./hal/star_$F

    mkdir -p ./lyr/star_$F
    STAR --genomeDir $LYRGENDIR --readFilesCommand zcat --readFilesIn $F\_R1.fastq.gz $F\_R2.fastq.gz \
        --outFileNamePrefix star_$F- --runThreadN $CPU --genomeLoad NoSharedMemory \
        --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
        --outSJfilterCountUniqueMin 3 2 2 2 --outMultimapperOrder Random \
        --outFilterType BySJout --outStd SAM | samtools view -Shb - > ./lyr/$F.bam
    samtools sort -o ./lyr/$F.refsort.bam ./lyr/$F.bam 
    samtools index -c ./lyr/$F.refsort.bam
    mv star_$F-* ./lyr/star_$F
done

## EAGLE-RC: read classification and the quantification with featureCounts

cd hal
for i in `ls *.refsort.bam`; do 
    F=`basename $i .refsort.bam`
    eagle -t 8 -a $F.refsort.bam -r ../$HALREF.fa -v ../H.vs.L.gtf.vcf --splice --rc 1> $F.H.vs.L.txt 2> $F.H.vs.L.readinfo.txt
    eagle-rc -a $F.refsort.bam --listonly -o $F.H.vs.L -v $F.H.vs.L.txt $F.H.vs.L.readinfo.txt > $F.H.vs.L.list
done
cd ..
cd lyr
for i in `ls *.refsort.bam`; do 
    F=`basename $i .refsort.bam`
    eagle -t 8 -a $F.refsort.bam -r $LYRREF.fa -v ../L.vs.H.gtf.vcf --splice --rc 1> $F.L.vs.H.txt 2> $F.L.vs.H.readinfo.txt
    eagle-rc -a $F.refsort.bam --listonly -o $F.L.vs.H -v $F.L.vs.H.txt $F.L.vs.H.readinfo.txt > $F.L.vs.H.list
done
cd ..

mkdir -p eagle
for i in `ls *_R1.fastq.gz`; do
    F=`basename $i _R1.fastq.gz`
    python scripts/ref2_consensus.py --pe -u -o eagle/$F.ref \
        -A hal/$F.H.vs.L.list \
        -B lyr/$F.L.vs.H.list
    eagle-rc --refonly --readlist -a hal/$F.refsort.bam -o eagle/$F.H eagle/$F.ref.chrA.list
    eagle-rc --refonly --readlist -a lyr/$F.refsort.bam -o eagle/$F.L eagle/$F.ref.chrB.list
    featureCounts -T 8 -t exon -g transcript_id -a $HALGTF.gtf -o eagle/$F.H.counts.txt eagle/$F.H.ref.bam 
    featureCounts -T 8 -t exon -g transcript_id -a $LYRGTF.gtf -o eagle/$F.L.counts.txt eagle/$F.L.ref.bam 
done

python scripts/tablize.py -skip 1 -a -i 0 -c 6 eagle/*.H.counts.txt > eagle.H.tsv
python scripts/tablize.py -skip 1 -a -i 0 -c 6 eagle/*.L.counts.txt > eagle.L.tsv

# Homeolog counts in terms of halleri gene id, transcript level
python scripts/tablize.py -a H.vs.L.reciprocal_best eagle.H.tsv | cut -f 1,3- | sort -k1 > eagle.H.homeolog.tsv
python scripts/tablize.py -a L.vs.H.reciprocal_best eagle.L.tsv | cut -f 2,3- | sort -k1 > eagle.L.homeolog.tsv

# Homeolog counts in terms of halleri gene id, gene level by summation of transcript counts
perl -ne 'chomp; @t=split(/\s+/); @i=split(/\./, $t[0]); $a=join("\t", @t[1..$#t]); print "$i[0]\t$a\n";' eagle.H.homeolog.tsv > temp.txt
python scripts/tablize.py -add temp.txt > eagle.H.homeolog.genelevel.tsv
perl -ne 'chomp; @t=split(/\s+/); @i=split(/\./, $t[0]); $a=join("\t", @t[1..$#t]); print "$i[0]\t$a\n";' eagle.L.homeolog.tsv > temp.txt
python scripts/tablize.py -add temp.txt > eagle.L.homeolog.genelevel.tsv

# Subgenome unique mapped reads
for i in `ls *_R1.fastq.gz`; do
    F=`basename $i _R1.fastq.gz`
    echo "" > dummy.txt
    eagle-rc --refonly --readlist -a hal/$F.refsort.bam -u lyr/$F.refsort.bam -o eagle/$F.H.only dummy.txt
    eagle-rc --refonly --readlist -a lyr/$F.refsort.bam -u hal/$F.refsort.bam -o eagle/$F.L.only dummy.txt
    featureCounts -T 8 -t exon -g transcript_id -a $HALGTF.gtf -o eagle/$F.H.only.counts.txt eagle/$F.H.only.ref.bam
    featureCounts -T 8 -t exon -g transcript_id -a $LYRGTF.gtf -o eagle/$F.L.only.counts.txt eagle/$F.L.only.ref.bam
done

python scripts/tablize.py -skip 1 -a -i 0 -c 6 eagle/*.H.only.counts.txt > eagle.H.only.tsv
python scripts/tablize.py -skip 1 -a -i 0 -c 6 eagle/*.L.only.counts.txt > eagle.L.only.tsv

# Subgenome unique mapped reads in subgenome unique genes
python scripts/tablize.py -a H.only.list eagle.H.only.tsv > eagle.H.only.unique.tsv
python scripts/tablize.py -a L.only.list eagle.L.only.tsv > eagle.L.only.unique.tsv
