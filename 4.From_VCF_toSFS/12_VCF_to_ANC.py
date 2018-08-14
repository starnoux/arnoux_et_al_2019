#!/usr/bin/python2.7

import os
import sys
import argparse
import re

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

### TRY TO GET THE FORMAT TAB FOR THE 4P ## 
courriel = "stephanie.arnoux@inra.fr "
 
parser=argparse.ArgumentParser(                                                           #creation of a parser which gather the arguments of the program
    description=""" This program will generate sequences corresponding to the Ancestral Allele from the vcf file of an ancestral accession. """+courriel+"""\n\n""")

#parser.add_argument("-l",type=str,default="P",                                              help="prefix letter of the LD or LW or LO (default: %(default)s)")
parser.add_argument("-i",           type=str,       default="/data_nas1/wgs/solanaceae/Solution/starnoux/Tools/4p-master/OUTGROUPS/vcf_files/",                         help="input directory (default: %(default)s)!")
parser.add_argument("-tabout",           type=str,       default="LA_MD_LA2951.tab",                         help="input name (default: %(default)s)!")
parser.add_argument("-vcfin",           type=str,       default="LA_MD_LA2951.recode.vcf",                         help="output name (default: %(default)s)!")
parser.add_argument("-chrom",           type=str,       default="Solyc",                         help="output name (default: %(default)s)!")

args=parser.parse_args()

with open (args.i+args.vcfin, "r+b" ) as vcf :
	with open(args.i+args.tabout, 'w+') as target :
		for line in vcf:
			if line[0]=='#':
				pass
			else :
				#line_out = ""
				lstcol = line.split("\t")
				# Here we defined the column for each line. (list of columns in the line)
				chrom = lstcol[0]
				pos = lstcol[1]
				ref = lstcol[3]
				alt = lstcol[4]
				het = lstcol[9].split(":")[0]
				if args.chrom in chrom[:len(args.chrom)]:
					Nchrom = chrom[len(args.chrom):len(args.chrom)+2]
				#print(Nchrom)
					if het == "1/1" :
						line_out = Nchrom + '\t' + chrom + '_' + pos + '\t' + alt + '\n'
					elif het == "0/0" :
						line_out = Nchrom + '\t' + chrom + '_' + pos + '\t' + ref + '\n'
					else :
						line_out = ""
					#print(het)
					#print(line_out)
					target.write(line_out)


vcf.close()
target.close()
# time to close the files, it is done baby!
