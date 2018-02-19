#! /bin/bash
## Example workflows using EAGLE

CPU=8

## -- Germline Variant Calling --
# Directly after variant calling, run eagle to calculate variant log likelihood ratios
eagle -t $CPU -v variants.vcf -a align.bam -r reference.fa > eagle.txt

# Filter variants based on eagle output
python scripts/compile_likelihoods.py -p eagle.txt -minlr 5 -minaf 0.3 -mindepth 10 -seen > eagle.filter.txt

# Filter the VCF and add eagle score as an annotation
python scripts/combine_vcf_eagle.py -v variants.vcf -e eagle.filter.txt > variants.eagle.filter.vcf

## -- Somatic Variant Calling --
# Run variant calling on tumor samples, 
# then run eagle to calculate log likelihood ratios for both normal and tumor samples
eagle -t $CPU -v tumor.vcf -a tumor.bam -r reference.fa > tumor.txt
eagle -t $CPU -v tumor.vcf -a normal.bam -r reference.fa > normal.txt

# Filter variants based on eagle outputs, where positive sample is tumor and negative sample is normal
python scripts/compile_likelihoods.py -p tumor.txt -n normal.txt -minlr 5 -minaf 0.05 -maxlr -2 -maxaf 0.04  -mindepth 10 -seen > somatic.filter.txt

# Filter the VCF and add eagle score as an annotation
python scripts/combine_vcf_eagle.py -v tumor.vcf -e somatic.filter.txt > variants.somatic.vcf

## -- Trio De Novo Variant Calling --
# Run variant calling on each of child and parents separately, then make a union vcf
python ~/scripts/tablize.py -i 0-4 child.vcf father.vcf mother.vcf | cut -f 1-5 > union.vcf
# Add headers to union.vcf [optional]
grep '^#' child.vcf | cat - union.vcf > tmp && mv tmp union.vcf

# then run eagle to calculate log likelihood ratios for both all samples using union.vcf
eagle -t $CPU -v union.vcf -a child.bam -r reference.fa > child.txt
eagle -t $CPU -v union.vcf -a father.bam -r reference.fa > father.txt
eagle -t $CPU -v union.vcf -a mother.bam -r reference.fa > mother.txt

# Filter variants based on eagle outputs, where positive sample is child and negative sample is parent
python scripts/compile_likelihoods.py -p child.txt -n father.txt -minlr 5 -minaf 0.25 -maxlr -2 -maxaf 0.04  -mindepth 10 -seen > child.vs.father.txt
python scripts/compile_likelihoods.py -p child.txt -n mother.txt -minlr 5 -minaf 0.25 -maxlr -2 -maxaf 0.04  -mindepth 10 -seen > child.vs.mother.txt

# Get the variants which are present only in the child
python scripts/tablize.py -a -i 0-7 child.vs.father.txt child.vs.mother.txt > child.vs.parents.txt

# Filter the VCF and add eagle score as an annotation
python scripts/combine_vcf_eagle.py -v child.vcf -e child.vs.parents.txt > child.denovo.vcf

## -- Read Classification for RNA-Seq expression quantification with homeologs --
## Requires LAST (http://last.cbrc.jp/)
## Align the sample data to each origin specific reference separately

## From origin specific transcripts in an allopolyploid, find homeologs and their genotype differences
lastdb -uNEAR -R01 A_origin A_transcripts.fa
lastdb -uNEAR -R01 B_origin B_transcripts.fa

lastal A_origin -P$CPU B_transcripts.fa | last-map-probs -m 0.49 > A_origin.maf
lastal B_origin -P$CPU A_transcripts.fa | last-map-probs -m 0.49 > B_origin.maf

# [Note] check to make sure your transcripts are based on exons or CDS in the annotation. If CDS then change below, -f exon to -f CDS
python scripts/homeolog_genotypes.py -f exon -o Ref_A -g annotation_A.gtf A_origin.maf B_origin.maf
python scripts/homeolog_genotypes.py -f exon -o Ref_B -g annotation_B.gtf B_origin.maf A_origin.maf

eagle -t $CPU -a data_align2_A.bam -r A.reference.fa -v Ref_A.gtf.vcf --omega=1e-40 --mvh --splice --isc --verbose 1> Ref_A.sample.txt 2> Ref_A.sample.readinfo.txt
eagle-rc -a data_align2_A.bam --listonly -o Ref_A.sample Ref_A.sample.txt Ref_A.sample.readinfo.txt > Ref_A.sample.list

eagle -t $CPU -a data_align2_B.bam -r B.reference.fa -v Ref_B.gtf.vcf --omega=1e-40 --mvh --splice --isc --verbose 1> Ref_B.sample.txt 2> Ref_B.sample.readinfo.txt
eagle-rc -a data_align2_B.bam --listonly -o Ref_B.sample Ref_B.sample.txt Ref_B.sample.readinfo.txt > Ref_B.sample.list

## Find the consensus classification based on likelihood
python scripts/ref2_consensus.py -A Ref_A.sample.list -B Ref_B.sample.list -o sample

## Write bam files based on consensus list, using A as the reference
eagle-rc -a data_align2_A.bam -o sample.chrA --readlist sample.chrA.list
eagle-rc -a data_align2_B.bam -o sample.chrB --readlist sample.chrB.list

## Perform read counting as you prefer, for example:
featureCounts -T $CPU -t exon -g transcript_id -a annotation_A.gtf -o sample.A.counts.txt sample.ref.bam
featureCounts -T $CPU -t exon -g transcript_id -a annotation_A.gtf -o sample.B.counts.txt sample.alt.bam
