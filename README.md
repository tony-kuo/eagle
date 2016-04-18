# EAGLE: Explicit Alternative Genome Likelihood Evaluator

Implemented in python [Compatible with 2\.7\.\*, 3\.\*]

**Requires**: pysam, cffi, numpy, scipy

**Inputs**

1. Alignment data in coordinated sorted indexed BAM format (\*\.bam, \*\.bam.bai)
2. Reference genome in FASTA format
3. Variant candidates in VCF format

**Output**

A tab-delimited text file with one row per variant and columns representing:

1. chromosome / sequence id
2. coordinate position
3. reference sequence
4. alternative sequence
5. total number of reads
6. number of reads supporting alternative sequence
7. log 10 probability
8. log 10 likelihood ratio (odds)
9. nearby variants in the set of hypotheses if any

**Input/Output Parameters**

-v VCF file describing the variants, only the columns describing position and sequence are used [columns: 1,2,4,5].
-a BAM alignment data file, reference coordinated sorted with index [*filename*.bam.bai].
-r Reference genome [multi] fasta file.
-o Output file name, defaults to *stdout* if not specified.

**Program Parameters**

-t The number of threads to use.

-n Group nearby variants within *n* bp of each other to be considered in the set of hypotheses for marginal probability calculations.

-k Typically, 2^*n* combinations are checked for *n* variants in the hypotheses set. Limit the number of combinations to 2^sqrt(*n*) if the number of variants exceeds *k*.

-mvh Instead of considering the combinations of variants in the hypotheses set, consider that all variants in the set co-occur by testing **one** multi-variant hypothesis.

-p Use primary alignments only. This will ignoring multi-mapping considerations.
