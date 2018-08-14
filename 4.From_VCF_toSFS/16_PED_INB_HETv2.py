#!/usr/bin/python2.7

import os
import sys
import argparse
import re
from itertools import islice
import random 

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

### TRY TO GET HETEROZIGOUS SITES FOR INBREED POP ## 

courriel = "stephanie.arnoux@inra.fr "
 
parser=argparse.ArgumentParser(                           #creation of a parser which gather the arguments of the program
  description=""" This program will generate a ped file with non inbreed lines for .ped. """+courriel+"""\n\n""")

#parser.add_argument("-l",type=str,default="P",                       help="prefix letter of the LD or LW or LO (default: %(default)s)")
parser.add_argument("-i",      type=str,    default="/data_nas1/wgs/solanaceae/Solution/starnoux/Tools/4p-master/",             help="input directory (default: %(default)s)!")
parser.add_argument("-ped",      type=str,    default="LA_CSFILTER_MD_A.noMiss.pruned.ped",             help="input ped file with  (default: %(default)s)!")
parser.add_argument("-InPop",      type=str,    default="CROP",             help="Population inbreed to change(default: %(default)s)!")
parser.add_argument("-Nosex",      type=str,    default="LA_CSFILTER_MD_A.noMiss.pruned.nosex",             help="input nosex file with the family and IDs -> FID tab IDD ! (default: %(default)s)!")
parser.add_argument("-outputped",      type=str,    default="LA_CSFILTER_MD_A.noMiss.pruned.HET.ped",             help="output name (default: %(default)s)!")


args=parser.parse_args()

with open (args.i+args.ped, "r+b" ) as ped, open (args.i+args.Nosex, "r+b" ) as Nosex, open(args.i+args.outputped, 'w') as outputped:
	dictPED = {}
	for line in ped:
		IDs = line.rstrip().split('\t')[1]
		DNA = line.rstrip().split('\t')[0:len(line.split('\t'))]
		dictPED.setdefault(IDs, []).append(DNA)
	#keydictPED = dictPED.keys()
	#valdictPED = dictPED.values()
	## This dictionnary is made of the IDs only and the sequences!
	dictFAM = {}
	for line in Nosex:
		FAM = line.rstrip().split('\t')[0]
		IDD = line.rstrip().split('\t')[1]
		dictFAM.setdefault(FAM, []).append(IDD)
	Inb_Acc = dictFAM[args.InPop]
	for i in range(0,len(Inb_Acc),2):
		Acc1 = Inb_Acc[i]
		if i+1 < len(Inb_Acc):
			Acc2 = Inb_Acc[i+1]
			line_out = args.InPop + '\t' + Acc1 +'_' + Acc2 + '\t0\t0\t0\t-9' 
			for j in range(6, len(dictPED[Acc1][0][0:]), 2):
				line_out = line_out + '\t' + random.choice([dictPED[Acc1][0][j],dictPED[Acc1][0][j+1]]) + '\t' + random.choice([dictPED[Acc2][0][j],dictPED[Acc2][0][j+1]]) 
			#print(len(line_out.split('\t')))
		else:
			line_out = '\t'.join(map(str,dictPED[Acc1][0][0:]))
			#print(len(line_out.split('\t'))) # no i+! so print it all
		outputped.write(line_out + '\n')
	for key, value in dictFAM.items():
		if key in args.InPop:
			pass
		else:
			for x in value:
				line_out = '\t'.join(map(str,dictPED[x][0][0:]))
				#print(len(line_out.split('\t')))
				outputped.write(line_out + '\n')


ped.close()
Nosex.close()
outputped.close()




########################### DO THE LOOP NOW ! #####################
#def search(values, searchFor):
#	for k in values:
#		for v in values[k]:
#			if searchFor in v:
#				return k
#	return None
