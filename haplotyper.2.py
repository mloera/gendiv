#!/usr/bin/python

import collections
import sys
from bam_parser import parse_bam
ingenotypes = sys.argv[1]
snps = sys.argv[2]

complement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}

snps_dict = {}

with open( snps, 'r') as infile:

        lineas = infile.readlines()

        for linea in lineas:

                locus = linea.split('\t')[1]
                snps_dict[locus] = []

with open( snps, 'r') as infile:

	lineas = infile.readlines()

	for linea in lineas:
		if 'POS' not in linea:
			locus = linea.split('\t')[1]
			pos = int(linea.split('\t')[2])
			major_allele = linea.split('\t')[5].split(":")[0]
			minor_allele = linea.split('\t')[6].split(":")[0]
			snps_dict[locus].append(tuple([pos, major_allele, minor_allele]))

# print(snps_dict)

haplotypes = {}
nucs = ['A', 'T', 'G', 'C']

bam = parse_bam( ingenotypes )

aligned_loci = set([ r.split(',')[0] for r in bam ])
# print(aligned_loci, dict(collections.Counter(aligned_loci)))

for locus in aligned_loci:
	haplotypes[ locus ] = []

for r in bam:

	locus = r.split(',')[0]
	read =  r.split(',')[1]
	aligned_seq = r.split(',')[2]
	offset = int(r.split(',')[3])

	haplotype = []
	
	if locus in snps_dict.keys():	
		
		for snp in snps_dict[ locus ] :
		
			try:
#			if snp[0] < len( aligned_seq )  and snp[0] - offset >= 0 and snp[0] - offset < len( aligned_seq ):
			
				#if aligned_seq[ snp[0] - 1] in [ snp[1], snp[2] ]:
				haplotype.append( aligned_seq[ snp[0] - offset -1] )

				#else:
				#	haplotype.append( '-' )
			except IndexError:
				haplotype.append('-')
	haplotypes[locus].append( ''.join( haplotype ) )



for hap in haplotypes:
	hap_counts = dict(collections.Counter( haplotypes[hap] ) )
	max_count = max( [hap_counts[c] for c in hap_counts ] )
	for r in hap_counts:
		if hap_counts[r] / max_count >= float(sys.argv[3]):
			print(hap, r, hap_counts[r] )


