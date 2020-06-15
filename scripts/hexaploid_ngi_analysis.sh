#! /bin/bash
# Example workflow using EAGLE-RC in ngi mode for hexaploid wheat
# LAST at: http://last.cbrc.jp/ 
# STAR at: https://github.com/alexdobin/STAR
# featureCounts at: http://bioinf.wehi.edu.au/featureCounts/

# The lyrata gene and scaffold ids should be modified to prepend Alyr_ or some other way make ids unique
REF=Taes_genome_2017_05
GTF=iwgsc_refseqv1.0_HighConf_UTR_2017May05
CPU=8

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
~/eagle/eagle-rc --ngi --listonly --splice --isc --ref1=$REF.chrA.fa --ref2=$REF.chrB.fa --bam1=$F.chrA.refsort.bam --bam2=$F.chrB.refsort.bam > $F.AvsB.list
~/eagle/eagle-rc --ngi --listonly --splice --isc --ref1=$REF.chrA.fa --ref2=$REF.chrD.fa --bam1=$F.chrA.refsort.bam --bam2=$F.chrD.refsort.bam > $F.AvsD.list
~/eagle/eagle-rc --ngi --listonly --splice --isc --ref1=$REF.chrB.fa --ref2=$REF.chrD.fa --bam1=$F.chrB.refsort.bam --bam2=$F.chrD.refsort.bam > $F.BvsD.list

mkdir -p eagle
for i in `ls *_R1.fastq.gz`; do
    F=`basename $i _R1.fastq.gz`
    python ref_ngi_consensus.py --pe -u -d -o eagle/$F.ref -AB chrA/$F.A.vs.B.list -AD chrA/$F.A.vs.D.list -BD chrB/$F.B.vs.D.list
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
