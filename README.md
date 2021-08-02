# gendiv
These set of scripts take as input a BAM file and a list of SNPs (output from VCFtools --freq) and prints out a list of SNP-based haplotypes and haplotype counts. The script takes also as input a haplotype frequency treshold. The haplotype frequency threshold is defined in proportion to the maximum haplotype count for each chromosome in the BAM file.

Dependencies:
PySAM 0.16.0.1

Usage:

```
python haplotyper.2.py input.bam SNP_frq.txt threshold
```
