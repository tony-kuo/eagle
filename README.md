# EAGLE: Explicit Alternative Genome Likelihood Evaluator

**Requires**: htslib (http://www.htslib.org/). Set HTSDIR in the make file to the htslib folder and make.  Note that we merely call make on htslib, as such, its dependencies and system requirements need to be fulfilled.

Compile:

`git clone https://github.com/samtools/htslib.git`

`make`

Usage: 

`eagle -v variants.vcf -a alignment.bam -r reference.fasta > output.tab`

**scripts/example.sh describe some EAGLE workflows as a bash script**

**For EAGLE read classification options, see EAGLE-RC below**

### Inputs

1. Alignment data in coordinated sorted indexed BAM format (\*\.bam, \*\.bam.bai)
2. Reference genome in FASTA format
3. Variant candidates in VCF format

### Output

A tab-delimited text file with one row per variant and columns representing:

1. chromosome / sequence id
2. coordinate position
3. reference sequence
4. alternative sequence
5. total number of reads seen
6. number of reads supporting reference sequence
7. number of reads supporting alternative sequence
8. log10 probability
9. log10 likelihood ratio (odds)
10. variants in the set of nearby variant if any, otherwise []

The read counts represent reads that are unambiguously for the reference or alternative sequence (2x the probability favoring one over the other), as opposed to aligned reads. Our model attempts to account for various uncertainties in the hypotheses.

### Input/Output Parameters

**-v --vcf**  [FILE] VCF file describing the variants, only the columns describing position and sequence are used [columns: 1,2,4,5].

**-a --bam**  [FILE] BAM alignment data file, reference coordinated sorted with index [*filename*.bam.bai].

**-r --ref**  [FILE] Reference genome [multi] fasta file.

**-o --out**  [FILE] Output file name, defaults to *stdout* if not specified.

### Program Parameters

**-t --nthread** [INT]  The number of processes to use.  Default is 1.

**-s --sharedr** [INT]  Group/chain nearby variants based on shared reads.  Default is 0/off, which uses the distance based method.  Option 1 will group variants if any read that crosses the first variant also cross the variant under consideration.  Option 2 will group variants if any read that crosses any variant in the set also cross the variant under consideration.

**-n --distlim** [INT]  Group/chain nearby variants within *n* bp of each other to be considered in the set of hypotheses for marginal probability calculations.  Default is 10 bp (0 is off).

**-w --maxdist** [INT]  Maximum number of bases between any two variants in the set of hypotheses. This sets a window size that allows for larger values of -n without chaining an excessive number of variants.  Default is 0 bp (off).

**-m --maxh** [INT]  The maximum number of hypotheses to be tested.  Instead of looking at all 2^n combinations for a set of variants, if after the current *k* for *n choose k* combinations is considered and the number of tested hypthotheses exceeds the maximum, then do not consider more combinations.  The solo variants and the "all variants" hypothesis are always tested first and do not count towards the maximum, thus any value greater than 0 ensures all doubles are tested.  Default is 1024 (2^10).

**--mvh**  Instead of the marginal probabilities over the hypotheses set, output only the maximum variant hypothesis (highest probability of phased variants) in the set of hypotheses.  Note that **maxh** will limit the possible combinations tested.

**--pao**  Use primary alignments only, based on SAM flag.

**--isc**  Ignore soft-clipped bases in reads when calculating the probabilities, based on cigar string.

**--nodup**  Ignore marked duplicate reads, based on SAM flag.

**--splice**  Allow spliced reads, based on cigar string.

**--dp**  Instead of the short read model, which assumes no indel errors, use dynamic programming (short in long) to calculate the likelihood.  This allows handling of long reads which have higher rates of sequencing errors and indel errors.

**--gap\_op** [INT]  Gap open penalty for use with --dp.  Default is 6.  For long reads that contain indel errors, 2 may be a better.

**--gap\_ex** [INT]  Gap extend penalty for use with --dp.  Default is 1.

**--verbose**  Verbose mode.  Output the likelihoods for every read seen for every hypothesis to *stderr*.  Used in read classification with **eagle-rc**.

**--lowmem**  Low memory usage mode.  For SNPs, we use a method to quickly derive the alternative hypothesis probability from the reference hypothesis probability without constructing the alternative sequence in memory.  For indels, which can be treated as a series of SNPs, this method may not be faster depending on read depth due to the number of frameshifted bases to account for.  Though it will save memory which may allow for more threads without hitting some memory cap.

**--phred64**  Reads quality scores are in phred64.  Default is phred33.

**--hetbias** [FLOAT]  Prior probability bias towards heterozygous or homozygous mutations.  Value between [0,1] where 1 is towards heterozygosity.  Default is 0.5 (unbiased).

**--omega** [FLOAT]  Prior probability of originating from outside paralogous source (i.e. not from reference genome and also not from the candidate variant genome).  Value between [0,1].  Default is 1e-5.

**--rc**  Wrapper for read classification settings: --omega=1.0e-40 --isc --mvh --verbose --lowmem.

### Usage Notes

*compare2TruthData.py*: Separate false positives and true positives based on truth data given as a VCF. 

*compile\_likelihoods.py*: Post-process the probabilities calculated by EAGLE and compile the results into tab-delimited table format.  It can be used to find somatic mutations given positive (i.e. tumor) and negative (i.e. normal) results on the same set of variants. Likelihood ratio and allele frequency thresholds are then used to filter mutations.

*combine\_vcf\_eagle.py*: Annotate the vcf file with the results from *compile\_likelihoods.py*. Raw EAGLE output can also be used if specified in the options.

Heterozygous non-reference variants (VCF: comma separated multiple alternative sequences) are output as separate entries. Furthermore, if they are near other variants, it will separately consider each alternative sequence in their own sets so that phasing is not an issue. This may result in entries with the first 4 columns being identical. The user will need to examine the variant set of these cases to determine which it belongs to. The script *compileLikelihood.py* will naively retain the one with the maximum probability, for entries with identical coordinates.

## EAGLE-RC

EAGLE-RC is a method for classifying whether a read belongs to one genomic hypothesis or another, given a set of genotype differences between them.  This can be applicable for determining if reads originate from a specific allele or from a specific homeolog in allopolyploids.  

EAGLE-RC can also classify reads directly via calculating the likelihood of the read given two alignments to different (sub)genome references, without knowing the genotype differences explicitly, implemented as *no genotype info* (--ngi) mode.

First with explicit genotype difference information, use EAGLE to calculate likelihoods for each read and for each hypothesis (--verbose) as well as phased variants likelihoods (--mvh) as output.  Then the program *eagle-rc* can take these two inputs to classify reads and optionally split reads into bam files for each class.  We also use a lower omega to be more tolerant to sequence differences outside the tested hypotheses.  The wrapper option (--rc) takes care of these parameters.  Other options (such as --dp for long reads, --splice for RNA-seq, --isc to ignore soft clipped bases, etc.) may or may not be applicable depending on the use case.

Usage, with more details in *example.sh*: 

`eagle -t 2 -v variants.vcf -a alignment.bam -r reference.fasta --rc 1> output.tab 2>readinfo.txt`

`eagle-rc -a alignment.bam -o out_prefix -v output.tab readinfo.txt > classified_reads.list`

### Program Parameters

**-v --var**  [FILE] EAGLE output file, containing variant likelihoods

**-a --bam**  [FILE] BAM alignment data file to be processed, reference coordinated sorted with index

**-o --out**  [String] Output file name prefix

**-u --unique**  [FILE1,FILE2,...] Optionally, also output reads that are unique against other BAM files (comma separated list)

**--listonly**  Print classified read list only (stdout) without processing BAM files

**--readlist**  Read from classified read list instead of EAGLE outputs and process BAM files.  This is useful if you previously outputed with *listonly* and then performed post-processing on list files to obtain, for example, a consensus list.  This will allow the post-processed list to be used instead.

**--refonly**  Write REF classified reads only when processing BAM file.

**--paired**  Consider paired-end reads together.

**--pao**  Use primary alignments only, based on SAM flag.

For no genotype information classification, the options in the default mode listed above are also applicable. Usage, where the classification is from the point of view of ref1 as the reference hypothesis and ref2 as the alternative hypothesis:

`eagle-rc --ngi -a align1.bam -o out1 --ref1=ref1.fa --ref2=ref2.fa --bam1=align1.bam --bam2=align2.bam > classified_reads.1vs2.list`

### Program Parameters specific to no genotype information mode (--ngi)

**--ref1**  [FILE] Reference genome 1 fasta file

**--bam1**  [FILE] Alignments to reference genome 1, --ref1, bam file

**--ref2**  [FILE] Reference genome 2 fasta file, we recommend that the chromosome names are different between ref1 and ref2.

**--bam2**  [FILE] Alignments to reference genome 2, --ref2, bam file

**--isc**  Ignore soft-clipped bases in reads when calculating the probabilities, based on cigar string.

**--nodup**  Ignore marked duplicate reads, based on SAM flag.

**--splice**  Allow spliced reads, based on cigar string.

**--phred64**  Reads quality scores are in phred64.  Default is phred33.

Without completely rerunning, we can get the classified output bam files for ref2 and align2.bam using the --readlist option on the classified reads list.  Keep in mind here that the resulting -.alt.bam contains the ref2 classified alignments.  One may wish to switch the -.alt.bam and -.ref.bam names after they are written.  If we wish to get the classified list for ref2 for further processing, then it is possible to switch REF and ALT and columns 5 and 6 to get the probabilities, though the alignment positions will not be for ref2.  The simplest option is to rerun it:

`eagle-rc --ngi -a align2.bam -o out2 --ref2=ref1.fa --ref1=ref2.fa --bam2=align1.bam --bam1=align2.bam > classified_reads.2vs1.list`

## References
Tony Kuo and Martin C Frith and Jun Sese and Paul Horton. EAGLE: Explicit Alternative Genome Likelihood Evaluator. BMC Medical Genomics. 11(Suppl 2):28. https://doi.org/10.1186/s12920-018-0342-1
