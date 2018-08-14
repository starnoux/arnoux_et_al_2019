#!/usr/bin/python2.7

import os
import sys
import argparse
import re

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

### BEFORE EVERYTHING GET A GFF FILE CONTANING ONLY THE mRNA
# grep 'ID=mRNA' ../../../ref_genomes/tomato/ITAG3.2_gene_models.gff > ITAG3.2_gene_models.gff.gene

### TRY TO GET THE FORMAT TAB FOR THE 4P ## 
courriel = "stephanie.arnoux@inra.fr "
 
parser=argparse.ArgumentParser(                                                           #creation of a parser which gather the arguments of the program
    description=""" This program will generate sequences corresponding to the Ancestral Allele from the vcf file of an ancestral accession. """+courriel+"""\n\n""")

#parser.add_argument("-l",type=str,default="P",                                              help="prefix letter of the LD or LW or LO (default: %(default)s)")
parser.add_argument("-i",           type=str,       default="/data_nas1/wgs/solanaceae/Solution/starnoux/Tools/4p-master/",                         help="input directory (default: %(default)s)!")
parser.add_argument("-gff",           type=str,       default="ITAG3.2_gene_models.gff.gene",                         help="input gff name, here only the mRNA should be kept (default: %(default)s)!")
parser.add_argument("-mapin",           type=str,       default="LA_CSFILTER_MD_A.noMiss.pruned.map",                         help="output name (default: %(default)s)!")
parser.add_argument("-temp",           type=str,       default="LA_CSFILTER_MD_A.noMiss.pruned.map.temp",                         help="temporary file name name (default: %(default)s)!")
parser.add_argument("-fin",           type=str,       default="LA_CSFILTER_MD_A.noMiss.pruned.map.fin",                         help="output name (default: %(default)s)!")
parser.add_argument("-chrom",           type=str,       default="Solyc",                         help="output name (default: %(default)s)!")

args=parser.parse_args()

with open (args.i+args.gff, "r+b" ) as gff, open(args.i+args.mapin, 'r+b') as mapin,  open(args.i+args.temp, 'w') as temp:
	toto=list(word.split(' \n')[0].split("\t")[1].split("_")[0] for word in mapin)
	for line in gff:
		lstcol = line.split("\t")
		# Here we defined the column for each line. (list of columns in the line)
		POS = lstcol[3]
		ID = lstcol[8].split("=")[1].split(";")[0]
		line_out = ID + "\t" + str(int(POS)) + "\n"
		if any(t in ID for t in toto):
			#print(line_out)
			temp.write(line_out)

			
gff.close()
mapin.close()
temp.close()


with open(args.i+args.mapin, 'r+b') as mapin,  open(args.i+args.temp, 'r+b') as temp, open(args.i+args.fin, 'w') as fin:
	list_of_genes = list(word.split('\n')[0] for word in temp)
	#print(list_of_genes)
	location_ID = []
	location_OfTheGene = []
	for line in mapin:
		line = line.rstrip().split('\t')
		ID_line = line[1].split("_")[0]
		Loc_on_gene = line[1].split("_")[1]
		Nchr = re.search(args.chrom+'(.+?)g', ID_line)
		if Nchr :
			Nchrom = Nchr.group(1)
			for gene in list_of_genes :
				ID = gene.split('\t')[0]
				loc = gene.split('\t')[1]
				if ID == ID_line:
					sum = float(loc) + float(Loc_on_gene)
					line_out = Nchrom + '\t' + line[1] + '\t' + '0' + '\t' + str(int(sum)) + '\n'
					fin.write(str(line_out))
					#print(line_out)
		else : 
			print('You need to give the real Chromosome name in the reference e.g. Solyc or Capana')


mapin.close()
temp.close()
fin.close()

os.remove(args.temp)
