# EAGLE: Explicit Alternative Genome Likelihood Evaluator

**Requires**: htslib (http://www.htslib.org/). Set HTSDIR in the make file to the htslib folder and make.

Usage: 

`eagle -v variants.vcf -a alignment.bam -r reference.fasta > output.tab`

### For EAGLE read classification, see EAGLE-RC below

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
5. total number of reads seen
6. number of reads supporting reference sequence
7. number of reads supporting alternative sequence
8. log 10 probability
9. log 10 likelihood ratio (odds)
10. variants in the set of nearby variant if any, otherwise []

The read counts represent reads that are unambiguously for the reference or alternative sequence (2x the probability favoring one over the other), as opposed to aligned reads. Our model attempts to account for various uncertainties in the hypotheses.

**Input/Output Parameters**

-v  VCF file describing the variants, only the columns describing position and sequence are used [columns: 1,2,4,5].

-a  BAM alignment data file, reference coordinated sorted with index [*filename*.bam.bai].

-r  Reference genome [multi] fasta file.

-o  Output file name, defaults to *stdout* if not specified.

**Program Parameters**

-t INT  The number of processes to use. Default is 1.

-s INT  Group/chain nearby variants based on shared reads.  Default is 0 which is off and uses the distance based method instead.  Option 1 will group variants if any read that crosses the first variant also cross the variant under consideration.  Option 2 will group variants if any read that crosses any variant in the set also cross the variant under consideration.  Option 2 will chain together a greater number of variants, such that a full test of all combinations will exceed maxh. Option 1 is typically more computationally efficient and sufficient as a more aggressive chaining than the distance based method.  Alternatively, setting the distance limit to a bigger value (i.e. a third the read length, read length window limit) will yield similar results.

-n INT  Group/chain nearby variants within *n* bp of each other to be considered in the set of hypotheses for marginal probability calculations. Default is 10 bp (0 for off).

-w INT  Maximum number of bases between any two variants in the set of hypotheses. This sets a window size that allows for larger values of -n without chaining an excessive number of variants. Default is 0 bp (off).

--maxh INT  The maximum number of hypotheses to be tested.  Instead of looking at all 2^n combinations for a set of variants, if after the current *k* for *n choose k* combinations is finished and the number tested exceeds the maximum, then do not consider more combinations.  The solo variants and the "all variants" hypothesis are always tested first and do not count towards the maximum. Default is 1024 (2^10).

--mvh  Instead of the marginal probabilities over the hypotheses set, output only the maximum variant hypothesis (highest probability) among variant combinations in the set of hypotheses.  Keep in mind that *maxh* will limit the possible combinations tested.

--pao  Use primary alignments only, as defined by the SAM flag. This will also ignore multi-mapping considerations.

--isc  Ignore soft-clipped bases in reads when calculating the probabilities.

--verbose  Verbose mode. Output the likelihoods for every read seen for each hypothesis. Outputs to stderr.

--hetbias FLOAT  Bias the prior probability towards heterozygous or homozygous mutations. Value between [0,1] where 1 is towards heterozygosity. Default is 0.5 (unbiased).

**Usage Notes**

*compare2TruthData.py* is a simple script to separate false positives and true positives based on truth data given as a VCF. 

*compileLikelihoods.py* post-processes the probabilities calculated by EAGLE and can be used to find somatic mutations given positive (i.e. tumor) and negative (i.e. normal) results on the same set of variants. Likelihood ratio and allele frequency thresholds are then used to filter mutations.

If one expects that most mutations are not homozygous (i.e. in heterogenous tumor samples), then one can choose to skew the prior probability towards heterozygous/heterogeneous mutations. Otherwise, very low allele frequency mutations (~0.05) will have low probability. However, unless one has a good estimate of the cell mixture ratio in hetergenous samples and tune the bias appropriately, this will likely increase type I errors, depending on the threshold chosen.

Heterozygous non-reference variants (VCF: comma separated multiple alternative sequences) are output as separate entries. Furthermore, if they are near other variants, it will separately consider each alternative sequence in their own sets so that phasing is not an issue. This may result in entries with the first 4 columns being identical. The user will need to examine the variant set of these cases to determine which it belongs to. The script *compileLikelihood.py* will naively retain the one with the maximum probability, for entries with identical coordinates.

## EAGLE-RC

For read classification.  Use EAGLE to calculate likelihoods for each read and for each hypothesis (--verbose) as well as phased variants (--mvh) as output.  Then the program readclassify can take these inputs and classify the reads and optionally read in a bam file and split the reads into bam files for each class.  We also lower probability omega to be more tolerant to sequence differences outside the tested hypotheses.

Usage: 

`eagle -t 2 -v variants.vcf -a alignment.bam -r reference.fasta --omega=1e-40 --mvh --pao --isc --verbose 1> output.tab 2>readinfo.txt`

`readclassify -a alignment.bam -o out_prefix output.tab readinfo.txt > classified_reads.list`

**Program Parameters**

-o  prefix for output BAM files.

--listonly  print classified read list only (stdout) without processing BAM files

--readlist  read from classified read list instead of EAGLE outputs, perhaps from a previous run or a merged file, and process BAM files.


*old Python version*

I recommend using the C implementation as it is much faster, but the legacy Python implementation has been left in the repository [compatible with 2\.7\.\*, 3\.\*] with ABI level C functions via CFFI. I can't promise this will be kept updated.

**Requires**: pysam, cffi, numpy

Usage: `python eagle.py -v variants.vcf -a alignment.bam -r reference.fasta > output.tab`

