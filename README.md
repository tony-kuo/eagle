# EAGLE: Explicit Alternative Genome Likelihood Evaluator

Implemented in python with ABI level C functions [Compatible with 2\.7\.\*, 3\.\*]

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
9. first variant in the set of variant hypotheses if any

**Input/Output Parameters**

-v VCF file describing the variants, only the columns describing position and sequence are used [columns: 1,2,4,5].

-a BAM alignment data file, reference coordinated sorted with index [*filename*.bam.bai].

-r Reference genome [multi] fasta file.

-o Output file name, defaults to *stdout* if not specified.

**Program Parameters**

-t The number of threads to use.

-n Group nearby variants within *n* bp of each other to be considered in the set of hypotheses for marginal probability calculations.

-maxh The maximum number of hypotheses to be tested.  Instead of looking at all 2^n hypotheses, if after the current *k* for *n choose k* combinations is finished and the number of hypotheses tested exceeds the maximum, then do not consider more combinations.  The solo variants and the "all variant" hypothesis are always tested first and do not count towards the maximum.  Default is 1025 (2^10) hypotheses.

-mvh Instead of considering the combinations of variants in the hypotheses set, consider that all variants in the set co-occur by testing **one** multi-variant hypothesis.

-p Use primary alignments only, as defined by the SAM flag. This will also ignore multi-mapping considerations.

**Usage Notes**

Heterozygous non-reference variants (VCF: comma separated multiple alternative sequences) are given as separate entries. Furthermore, if they are near other variants, the sets of variant combinations will separately consider each alternative sequence. This may result in entries with the first 4 columns being identical. The user needs to use their own judgement to interpret the relative effect sizes of duplicate entries in the output. The script *compileLikelihood.py* will retain the entry with the maximum probability.
