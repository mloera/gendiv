#!/usr/bin/python

import sys

ingenotypes = sys.argv[1]
snps = sys.argv[2]
#indie = sys.argv[3]

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

haplotypes = []
nucs = ['A', 'T', 'G', 'C']

with open( ingenotypes, 'r') as infile:

	lineas = infile.readlines()

	for linea in lineas:
		indie = linea.split(' ')[0]
		haplotype = ''
		locus = linea.split(' ')[1]
		offset = int(linea.split(' ')[5])
		seq = linea.split(' ')[7].strip()
		read = linea.split(' ')[2]
		cigar = linea.split(' ')[6]
		flag = int(linea.split(' ')[4])

		cigar_n = cigar.replace('X', '=').replace('I', '=').replace('D', '=').replace('N', '=').replace('M', '=').replace('S', '=').replace('H', '=').split('=')
		cigar_l = [i for i in cigar if not i.isdigit()]
		# print(cigar_n, cigar_l)
		cigar_list = []
		for i in range(0, len(cigar_l)):

			cigar_element = tuple([cigar_n[i], cigar_l[i]])
			cigar_list.append(cigar_element)

		# print(cigar_list)
		
		seq_tr = ''
		# print(locus, seq)

	
		for c in cigar_list:
				# try:

			if c[1] == 'M' or c[1] == '=' or c[1] == 'X':
				seq_tr += seq[0:int(c[0])]
				seq = seq[int(c[0]) :len(seq)]
				#continue

			elif c[1] == 'D' or c[1] == 'N':
				seq_tr += int(c[0]) * '-'
				#continue

			elif c[1] == 'I' or c[1] == 'S' or c[1] == 'H' :
				seq = seq[int(c[0]) :len(seq)]						
				#continue
				# except ValueError:
					# print(cigar_list)
		# print(seq_tr)
				
		# break
		
		if flag == '16' or flag == '1040':
			try:
				seq_tr = "".join([complement[x] for x in seq_tr])
			except KeyError:
				pass
		else:
			seq_tr = seq_tr			
		for snp in snps_dict[locus]:
			# print('We are at the SNP pos loop')
			pos = snp[0]
			
			try:
				#if seq_tr[ pos - offset ] in (snp[1], snp[2]):
				#	haplotype += seq_tr[ pos - offset ]
				haplotype += seq_tr[ pos - offset ]
				#else:
				#	haplotype += '*'
			except IndexError:
				# print('IndexError')
				haplotype += '*'
		# if '-' not in haplotype:
		print(tuple([indie, locus, read, flag, offset, haplotype]))


#for hap in haplotypes:
#	print('We are at the last step')
#	if '-' not in hap[3]:
#	print(hap)

		
			
