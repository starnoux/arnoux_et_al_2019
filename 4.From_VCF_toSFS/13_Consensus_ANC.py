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
parser.add_argument("-i",           type=str,       default="/data_nas1/wgs/solanaceae/Solution/starnoux/Tools/4p-master/OUTGROUPS/tab_files/",                         help="input directory (default: %(default)s)!")
parser.add_argument("-ancA",           type=str,       default="LA_RMD_LA4126_OUT_B.tab",                         help="input name (default: %(default)s)!")
parser.add_argument("-ancB",           type=str,       default="LA_RMD_LA2951_OUT_B.tab",                         help="output name (default: %(default)s)!")
parser.add_argument("-ancC",           type=str,       default="LA_RMD_Consensus_OUT_B.tab",                         help="output name (default: %(default)s)!")


args=parser.parse_args()

with open (args.i+args.ancA, "r+b" ) as ancA, open(args.i+args.ancB, 'r+b') as ancB, open(args.i+args.ancC, 'w') as Cons:
	same = set(ancA).intersection(ancB)
	same.discard('\n')
	for line in same:
		Cons.write(line)


ancB.close()
ancA.close()
Cons.close()

with open(args.i+args.ancC) as fin:
    content = sorted(fin)

with open(args.i+args.ancC, "w") as fout:
    fout.writelines(content)

fin.close()
fout.close()

