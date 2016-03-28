# EAGLE: Explicit Alternative Genome Likelihood Evaluator

Implemented in python [2.7, 3]

**Requires**: pysam, cffi, numpy, scipy

**Inputs**

1. Alignment data in coordinated sorted indexed BAM format (.bam, .bam.bai)
2. Reference genome in FASTA format
3. Variant candidates in VCF format

**Output**

A tab-delimited text file with one row per variant and columns representing:

1. chromosom / sequence id
2. coordinate position
3. reference sequence
4. alternative sequence
5. total number of reads
6. number of reads supporting alternative sequence
7. log 10 probability
8. log 10 likelihood ratio (odds)
9. nearby variants in the set of hypotheses if any

