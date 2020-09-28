# B. From a clean vcf to a 2D Allele Frequency Spectrum :  
## EXPLANATION UPDATED

### ref_genomes
 `*.fasta` or `*.fa`: 1 fasta == 1 gene.  
 `*.gff`: anotation file.  
 `*.gff.gene`: mRNA subset. Use dto find the START position of the gene on the chromosome, and produce the chromosome-coordinate system.  
  `*.tab`: within chromosomes coordinate system.  

#	step0a_____VCF_RAW_allSamples
	Coordinates: within genes.

	Raw VCFs.
	`satomato.snpSift.filtered_LOW.synonymous.vcf`
	`saeggplant.snpSift.filtered_LOW.synonymous.vcf`
	`sapepper.snpSift.filtered_LOW.synonymous.vcf`


#	step0b_____VCF_NoMiss_noChr0_allSamples
	Coordinates: within genes.

	VCFs excluding positions with at least one missing genotype (grep -v '\.\/\.') and non-aligned SNPs (grep -v 'SMEL_000g' or 'Capana00g' or 'Solyc00g').
	`grep -v '\.\/\.' ../VCF_RAW_allSamples/satomato.snpSift.filtered_LOW.synonymous.vcf > LA.synonymous.2019.noMiss.vcf`
	`grep -v 'Solyc00g' LA.synonymous.2019.noMiss.vcf > LA.synonymous.2019.noChr0_noMiss.vcf`
	#grep -v "#" *.noChr0_noMiss.vcf: 372521
	`grep -v '\.\/\.' ../VCF_RAW_allSamples/saeggplant.snpSift.filtered_LOW.synonymous.vcf > MM.synonymous.2019.noMiss.vcf`
	`grep -v 'SMEL_000g' MM.synonymous.2019.noMiss.vcf > MM.synonymous.2019.noChr0_noMiss.vcf`
	#grep -v "#" *.noChr0_noMiss.vcf: 46958
	`grep -v '\.\/\.' ../VCF_RAW_allSamples/sapepper.snpSift.filtered_LOW.synonymous.vcf > PM.synonymous.2019.noMiss.vcf`
	`grep -v 'Capana00g' PM.synonymous.2019.noMiss.vcf > PM.synonymous.2019.noChr0_noMiss.vcf`
	#grep -v "#" *.noChr0_noMiss.vcf: 76843


#	step1_____VCF_synonymous
	Coordinates: within genes ; vcf file: REF allele=#1 / ALT allele=#2.

	VCFs keeping only the list of Crop and Wild for each species (*.recode.vcf).
	`vcftools --vcf ../VCF_NoMiss_noChr0_allSamples/LA.synonymous.2019.noChr0_noMiss.vcf --keep sample_files/LA_dadi.list --out LA_synonymous_MD --recode`
	#grep -v "#" *.recode.vcf: 372521
	`vcftools --vcf ../VCF_NoMiss_noChr0_allSamples/MM.synonymous.2019.noChr0_noMiss.vcf --keep sample_files/MM_dadi.list --out MM_synonymous_MD --recode`
	#grep -v "# *.recode.vcf": 46958
	`vcftools --vcf ../VCF_NoMiss_noChr0_allSamples/PM.synonymous.2019.noChr0_noMiss.vcf --keep sample_files/PM_dadi.list --out PM_synonymous_MD --recode`
	#grep -v "#" *.recode.vcf: 76843


#	step2-3_____MD_BED
	Coordinates: within genes ; bim file: contiguous / REF allele=#2 / ALT allele=#1 / CHR column=gene name.

	Produce PLINK files with the original A1/A2 allele encoding, allowing extra chrom, converting spaces in sampleIDs (*.bed, *.bim, *.fam, *.log, *.nosex).
       `plink --vcf ../VCF_synonymous/LA_synonymous_MD.recode.vcf --keep-allele-order --vcf-idspace-to _ --make-bed --allow-extra-chr --out LA_synonymous_MD.recode.4P`
	#cat *.recode.4P.bim | wc -l: 372521
       `plink --vcf ../VCF_synonymous/MM_synonymous_MD.recode.vcf --keep-allele-order --vcf-idspace-to _ --make-bed --allow-extra-chr --out MM_synonymous_MD.recode.4P`
	#cat *.recode.4P.bim | wc -l: 46958
       `plink --vcf ../VCF_synonymous/PM_synonymous_MD.recode.vcf --keep-allele-order --vcf-idspace-to _ --make-bed --allow-extra-chr --out PM_synonymous_MD.recode.4P`
	#cat *.recode.4P.bim | wc -l: 76843

	Recode PLINK files into mostly tab-delimited files instead of mostly space-delimited files.
	`plink --bfile LA_synonymous_MD.recode.4P --recode --tab --out ../MD_PED/LA_synonymous_MD.recode. --allow-extra-chr`
	#cat *.recode..map | wc -l: 372521
	`plink --bfile MM_synonymous_MD.recode.4P --recode --tab --out ../MD_PED/MM_synonymous_MD.recode. --allow-extra-chr`
	#cat *.recode..map | wc -l: 46958
	`plink --bfile PM_synonymous_MD.recode.4P --recode --tab --out ../MD_PED/PM_synonymous_MD.recode. --allow-extra-chr`
	#cat *.recode..map | wc -l: 76843


#	step4_____MD_PED
	Coordinates: within genes ; map file: contiguous / CHR column="1" only.

	Change MAPs to get the good format.
	`mv LA_synonymous_MD.recode..map LA_synonymous_MD.recode..map.old`
	`gawk -vOFS="" '{print "1","\t", $1,"\t", $3,"\t", $4}' LA_synonymous_MD.recode..map.old > LA_synonymous_MD.recode..map`
	#cat *.recode..map | wc -l: 372521
	`mv MM_synonymous_MD.recode..map MM_synonymous_MD.recode..map.old`
	`gawk -vOFS="" '{print "1","\t", $1,"\t", $3,"\t", $4}' MM_synonymous_MD.recode..map.old > MM_synonymous_MD.recode..map`
	#cat *.recode..map | wc -l: 46958
	`mv PM_synonymous_MD.recode..map PM_synonymous_MD.recode..map.old`
	`gawk -vOFS="" '{print "1","\t", $1,"\t", $3,"\t", $4}' PM_synonymous_MD.recode..map.old > PM_synonymous_MD.recode..map`
	#cat *.recode..map | wc -l: 76843


#	step5_____MD_Pruning
	**THIS STEP WAS SKIPPED DUE TO REVIEWER SUGGESTION**
	LD-Pruning of the PED file.


#	step6_____MD_Filtered
	Coordinates: within genes ; bim file: NON-contiguous / REF allele=#2 / ALT allele=#1 (if monomorphic position: "0") / CHR column="1" only.

	Filter the list of genes under selection (dnds) and/or in low-recombining regions (linkage map). 
	`plink --file ../MD_PED/LA_synonymous_MD.recode. --exclude filter_files/LA_DNAsp_LD_ListToRemove.txt --make-bed --out LA_MD_Filtered_Syn --allow-extra-chr`
	#cat *.bim | wc -l: 285652
	`plink --file ../MD_PED/MM_synonymous_MD.recode. --exclude filter_files/MM_DNAsp_ListToRemove.txt --make-bed --out MM_MD_Filtered_Syn --allow-extra-chr`
	#cat *.bim | wc -l: 41330
	`plink --file ../MD_PED/PM_synonymous_MD.recode. --exclude filter_files/PM_DNAsp_LD_ListToRemove.txt --make-bed --out PM_MD_Filtered_Syn --allow-extra-chr`
	#cat *.bim | wc -l: 49034


#	step7_____4P-Files
	Coordinates: within genes ; map file: NON-contiguous / CHR column="1" only.

	Recode PLINK files into mostly tab-delimited files instead of mostly space-delimited files.
	`plink --bfile ../MD_Filtered/LA_MD_Filtered_Syn --recode --tab --out LA_MD_Filtered_Syn.4P --allow-extra-chr`
	#cat *.4P.map | wc -l: 285652
	`plink --bfile ../MD_Filtered/MM_MD_Filtered_Syn --recode --tab --out MM_MD_Filtered_Syn.4P --allow-extra-chr`
	#cat *.4P.map | wc -l: 41330
	`plink --bfile ../MD_Filtered/PM_MD_Filtered_Syn --recode --tab --out PM_MD_Filtered_Syn.4P --allow-extra-chr`
	#cat *.4P.map | wc -l: 49034

	Change PEDs to have tab instead of spaces in the genotype fields.
	`sed -i "s/ /\t/g" LA_MD_Filtered_Syn.4P.ped`
	`sed -i "s/ /\t/g" MM_MD_Filtered_Syn.4P.ped`
	`sed -i "s/ /\t/g" PM_MD_Filtered_Syn.4P.ped`

	Change MAPs to encode the within-gene coordinate system in the variant name with "_position".
	`mv LA_MD_Filtered_Syn.4P.map LA_MD_Filtered_Syn.4P.map.old`
	`gawk -vOFS="" '{print "1","\t", $2,"_",$4,"\t", $3,"\t", $4}' LA_MD_Filtered_Syn.4P.map.old > LA_MD_Filtered_Syn.4P.map`
	#cat *.4P.map | wc -l: 285652
	`mv MM_MD_Filtered_Syn.4P.map MM_MD_Filtered_Syn.4P.map.old`
	`gawk -vOFS="" '{print "1","\t", $2,"_",$4,"\t", $3,"\t", $4}' MM_MD_Filtered_Syn.4P.map.old > MM_MD_Filtered_Syn.4P.map`
	#cat *.4P.map | wc -l: 41330
	`mv PM_MD_Filtered_Syn.4P.map PM_MD_Filtered_Syn.4P.map.old`
	`gawk -vOFS="" '{print "1","\t", $2,"_",$4,"\t", $3,"\t", $4}' PM_MD_Filtered_Syn.4P.map.old > PM_MD_Filtered_Syn.4P.map`
	#cat *.4P.map | wc -l: 49034


#	step8_____RMD
	Coordinates: within chromosomes ; map file: NON-contiguous / CHR column="1" to "12". Only MAPs are modified.

	Change the coordinate system from "within genes" to "within chromosomes", and replace CHR column by the true chromosome id. 
	The temp file is created by the script (but removed at the end) and contains the START of the gene on the chromosome.
	`./11_gff_map_2_NewMapLA.py -i ~/ -mapin 4P-Files/LA_MD_Filtered_Syn.4P.map -temp RMD/LA_MD_Filtered_Syn.4P.map.temp -fin RMD/LA_MD_Filtered_Syn.4P.RMD.map -gff ref_genomes/ITAG3.2_gene_models.gff.gene -chrom Solyc`
	#cat *.4P.RMD.map | wc -l: 285652
	`./11_gff_map_2_NewMapMM.py -i ~/ -mapin 4P-Files/MM_MD_Filtered_Syn.4P.map -temp RMD/MM_MD_Filtered_Syn.4P.map.temp -fin RMD/MM_MD_Filtered_Syn.4P.RMD.map -gff ref_genomes/SMEL_V3.2016_11_01.gff.gene -chrom SMEL_0`
	#cat *.4P.RMD.map | wc -l: 41330
	`./11_gff_map_2_NewMapPM.py -i ~/ -mapin 4P-Files/PM_MD_Filtered_Syn.4P.map -temp RMD/PM_MD_Filtered_Syn.4P.map.temp -fin RMD/PM_MD_Filtered_Syn.4P.RMD.map -gff ref_genomes/Capsicum.annuum.L_Zunla-1_v2.0_genes.gff.gene -chrom Capana`
	#cat *.4P.RMD.map | wc -l: 49034


#	step9_____HET
	Coordinates: same as step7 ; map file: same as step7. Only PEDs and NOSEXs are modified.

	Produce PEDs by pooling two individuals rendered haploid by randomly sampling one of their allele. This is to take into account inbreeding.
	`./16_PED_INB_HETv2.py -i ~/ -ped 4P-Files/LA_MD_Filtered_Syn.4P.ped -InPop CROP -Nosex 4P-Files/LA_MD_Filtered_Syn.4P.nosex -outputped HET/LA_MD_Filtered_Syn.4P.CROP.HET.ped`
	`gawk -F "\t" 'BEGIN{OFS="\t"} {print $1,$2}' LA_MD_Filtered_Syn.4P.CROP.HET.ped > LA_MD_Filtered_Syn.4P.CROP.HET.nosex`
	#
	#
	`./16_PED_INB_HETv2.py -i ~/ -ped 4P-Files/MM_MD_Filtered_Syn.4P.ped -InPop CROP -Nosex 4P-Files/MM_MD_Filtered_Syn.4P.nosex -outputped HET/MM_MD_Filtered_Syn.4P.CROP.HET.ped`
	`gawk -F "\t" 'BEGIN{OFS="\t"} {print $1,$2}' MM_MD_Filtered_Syn.4P.CROP.HET.ped > MM_MD_Filtered_Syn.4P.CROP.HET.nosex`
	#
	`./16_PED_INB_HETv2.py -i ~/ -ped HET/MM_MD_Filtered_Syn.4P.CROP.HET.ped -InPop WILD -Nosex HET/MM_MD_Filtered_Syn.4P.CROP.HET.nosex -outputped HET/MM_MD_Filtered_Syn.4P.HET.ped`
	`gawk -F "\t" 'BEGIN{OFS="\t"} {print $1,$2}' MM_MD_Filtered_Syn.4P.HET.ped > MM_MD_Filtered_Syn.4P.HET.nosex`
	#
	#
	`./16_PED_INB_HETv2.py -i ~/ -ped 4P-Files/PM_MD_Filtered_Syn.4P.ped -InPop CROP -Nosex 4P-Files/PM_MD_Filtered_Syn.4P.nosex -outputped HET/PM_MD_Filtered_Syn.4P.CROP.HET.ped`
	`gawk -F "\t" 'BEGIN{OFS="\t"} {print $1,$2}' PM_MD_Filtered_Syn.4P.CROP.HET.ped > PM_MD_Filtered_Syn.4P.CROP.HET.nosex`
	#
	`./16_PED_INB_HETv2.py -i ~/ -ped HET/PM_MD_Filtered_Syn.4P.CROP.HET.ped -InPop WILD -Nosex HET/PM_MD_Filtered_Syn.4P.CROP.HET.nosex -outputped HET/PM_MD_Filtered_Syn.4P.HET.ped`
	`gawk -F "\t" 'BEGIN{OFS="\t"} {print $1,$2}' PM_MD_Filtered_Syn.4P.HET.ped > PM_MD_Filtered_Syn.4P.HET.nosex`


#	step10_____VCF
	Coordinates: within chromosomes ; vcf file: contiguous / REF allele=#1 / ALT allele=#2 (if monomorphic position: ".").

	Create TABs with the reference alleles.
	`grep -v "#" ../VCF_synonymous/LA_synonymous_MD.recode.vcf | gawk -vOFS="" '{print $1,"_", $2, "\t", $4}' > LA_ref.tab`
	`grep -v "#" ../VCF_synonymous/MM_synonymous_MD.recode.vcf | gawk -vOFS="" '{print $1,"_", $2, "\t", $4}' > MM_ref.tab`
	`grep -v "#" ../VCF_synonymous/PM_synonymous_MD.recode.vcf | gawk -vOFS="" '{print $1,"_", $2, "\t", $4}' > PM_ref.tab`

	Produce PEDs of homozygous diploids from each haplotype.
	`Manually run: 17a_ped_haploid.R / Infile: ../HET/LA_MD_Filtered_Syn.4P.CROP.HET / Outfile: LA_MD_Filtered_Syn.4P.CROP.HET-acc*.ped`
	`Manually run: 17a_ped_haploid.R / Infile: ../HET/MM_MD_Filtered_Syn.4P.HET / Outfile: MM_MD_Filtered_Syn.4P.HET-acc*.ped`
	`Manually run: 17a_ped_haploid.R / Infile: ../HET/PM_MD_Filtered_Syn.4P.HET / Outfile: PM_MD_Filtered_Syn.4P.HET-acc*.ped`

	Create VCFs from PEDs respecting the reference/alternate alleles and inbreeding solved.
	`plink --allow-extra-chr --ped LA_MD_Filtered_Syn.4P.CROP.HET-acc1.ped --map ../RMD/LA_MD_Filtered_Syn.4P.RMD.map --a2-allele LA_ref.tab 2 1 '#' --real-ref-alleles --recode vcf --out LA_synonymous_MD_HET-acc1`
	`plink --allow-extra-chr --ped LA_MD_Filtered_Syn.4P.CROP.HET-acc2.ped --map ../RMD/LA_MD_Filtered_Syn.4P.RMD.map --a2-allele LA_ref.tab 2 1 '#' --real-ref-alleles --recode vcf --out LA_synonymous_MD_HET-acc2`
	#grep -v "#" *-acc*.vcf | wc -l: 285652
	`plink --allow-extra-chr --ped MM_MD_Filtered_Syn.4P.HET-acc1.ped --map ../RMD/MM_MD_Filtered_Syn.4P.RMD.map --a2-allele MM_ref.tab 2 1 '#' --real-ref-alleles --recode vcf --out MM_synonymous_MD_HET-acc1`
	`plink --allow-extra-chr --ped MM_MD_Filtered_Syn.4P.HET-acc2.ped --map ../RMD/MM_MD_Filtered_Syn.4P.RMD.map --a2-allele MM_ref.tab 2 1 '#' --real-ref-alleles --recode vcf --out MM_synonymous_MD_HET-acc2`
	#grep -v "#" *acc*.vcf | wc -l: 41330
	`plink --allow-extra-chr --ped PM_MD_Filtered_Syn.4P.HET-acc1.ped --map ../RMD/PM_MD_Filtered_Syn.4P.RMD.map --a2-allele PM_ref.tab 2 1 '#' --real-ref-alleles --recode vcf --out PM_synonymous_MD_HET-acc1`
	`plink --allow-extra-chr --ped PM_MD_Filtered_Syn.4P.HET-acc2.ped --map ../RMD/PM_MD_Filtered_Syn.4P.RMD.map --a2-allele PM_ref.tab 2 1 '#' --real-ref-alleles --recode vcf --out PM_synonymous_MD_HET-acc2`
	#grep -v "#" *acc*.vcf | wc -l: 49034


#	step11_____toSFS
	Coordinates: within chromosomes.

	Produce SFS with 4P. It requires a config file: sumstat.par.
	`4P -f ../HET/LA_MD_Filtered_Syn.4P.CROP.HET.ped -i 0 -n 14 -s 285652 -m ../RMD/LA_MD_Filtered_Syn.4P.RMD.map -a LA_ALLConsensus_OUT.RMD.tab > LA.log`
	`cp AFS-U_CROP_WILD.dadi.txt AFS-LA-U_CROP_WILD.dadi.txt`
	#sum cells: 285652 [including 246317 with allele information]
	#crop=5 (pseudohaploids made diploids) ; wild=9 (diploids)
	`4P -f ../HET/MM_MD_Filtered_Syn.4P.HET.ped -i 0 -n 7 -s 41330 -m ../RMD/MM_MD_Filtered_Syn.4P.RMD.map -a MM_MD_Consensus_OUT.RMD.tab > MM.log`
	`cp AFS-U_WILD_CROP.dadi.txt AFS-MM-U_WILD_CROP.dadi.txt`
	#sum cells: 41330 [including 18861 with allele information]
	#crop=4 (pseudohaploids made diploids) ; wild=3 (pseudohaploids made diploids)
	`4P -f ../HET/PM_MD_Filtered_Syn.4P.HET.ped -i 0 -n 7 -s 49034 -m ../RMD/PM_MD_Filtered_Syn.4P.RMD.map -a PM_MD_Consensus_OUTALL.RMD.tab > PM.log`
	`cp AFS-U_WILD_CROP.dadi.txt AFS-PM-U_WILD_CROP.dadi.txt`
	#sum cells: 49034 [including 29563 with allele information]
	#crop=5 (pseudohaploids made diploids) ; wild=2 (pseudohaploids made diploids)

	Visualize the SFS, and eventually transform it.
	`Manually run: 17b_sfs_visualize.py / Infile: AFS-LA-U_CROP_WILD.dadi.txt / Outfile: AFS-LA-U_CROP_WILD.dadi.png`
	`Manually run: 17b_sfs_visualize.py / Infile: AFS-MM-U_WILD_CROP.dadi.txt / Outfile: AFS-MM-U_WILD_CROP.dadi.png`
	`Manually run: 17b_sfs_visualize.py / Infile: AFS-PM-U_WILD_CROP.dadi.txt / Outfile: AFS-PM-U_WILD_CROP.dadi.png`


#	step12_____toFasta
	Coordinates: within genes.

	Change VCFs to the within-gene coordinate system.
	`Manually run: 17c_vcf_coord.sh / Infile: LA_synonymous_MD_HET-acc*.vcf / Outfile: LA_synonymous_MD_HET-acc*-genes.vcf`
	#grep -v "#" *-acc*-genes.vcf | wc -l: 285652
	`Manually run: 17c_vcf_coord.sh / Infile: MM_synonymous_MD_HET-acc*.vcf / Outfile: MM_synonymous_MD_HET-acc*-genes.vcf`
	#grep -v "#" *-acc*-genes.vcf | wc -l: 41330
	`Manually run: 17c_vcf_coord.sh / Infile: PM_synonymous_MD_HET-acc*.vcf / Outfile: PM_synonymous_MD_HET-acc*-genes.vcf`
	#grep -v "#" *-acc*-genes.vcf | wc -l: 49034

	Produce Fasta from VCFs where individuals are homozygous diploids. Two identical sequences are created (:0 and :1), so a single one (:0) is conserved.
	`Manually run: 17d_vcf_2_fasta.sh / Infile: LA_synonymous_MD_HET-acc*-genes.vcf / Outfile: LA_synonymous_MD_HET-acc*-genes.fasta`
	#cat LA/tmp.pos | wc -l: 12717
	`Manually run: 17d_vcf_2_fasta.sh / Infile: MM_synonymous_MD_HET-acc*-genes.vcf / Outfile: MM_synonymous_MD_HET-acc*-genes.fasta`
	#cat MM/tmp.pos | wc -l: 11326
	`Manually run: 17d_vcf_2_fasta.sh / Infile: PM_synonymous_MD_HET-acc*-genes.vcf / Outfile: PM_synonymous_MD_HET-acc*-genes.fasta`
	#cat PM/tmp.pos | wc -l: 8814

	Rename headers according to PopPhyl.
	`SPE=LA ; awk -vOFS="" '/^>/{split($1,split1,"_") ; split(split1[3],split2,":") ; split(split1[1],split3,">") ; print ">" , split2[1] , "|" , "'"$SPE"'" , "_" , split3[2] , "|" , split1[2] , "|" , "allele_1" ; next}{print}' LA_synonymous_MD_HET-acc1-genes.temp > LA_synonymous_MD_HET-acc1-genes.fasta`
	`SPE=LA ; awk -vOFS="" '/^>/{split($1,split1,"_") ; split(split1[3],split2,":") ; split(split1[1],split3,">") ; print ">" , split2[1] , "|" , "'"$SPE"'" , "_" , split3[2] , "|" , split1[2] , "|" , "allele_1" ; next}{print}' LA_synonymous_MD_HET-acc2-genes.temp > LA_synonymous_MD_HET-acc2-genes.fasta`
	`cat LA_synonymous_MD_HET-acc1-genes.fasta LA_synonymous_MD_HET-acc2-genes.fasta > LA_synonymous_MD_HET-genes.fasta`
	`rm *-acc1-* ; rm *-acc2-*`
	#grep ">" *-genes.fasta | cut -d"|" -f1 | cut -d">" -f2 | sort -n | uniq | wc -l: 12717
	#grep -e '^$' *-genes.fasta -B1: ZERO BLANK LINE
	`SPE=MM ; awk -vOFS="" '/^>/{split($1,split1,"_") ; split(split1[4],split2,":") ; split(split1[1],split3,">") ; print ">" , split1[3] , "_" , split2[1] , "|" , "'"$SPE"'" , "_" , split3[2] , "|" , split1[2] , "|" , "allele_1" ; next}{print}' MM_synonymous_MD_HET-acc1-genes.temp > MM_synonymous_MD_HET-acc1-genes.fasta`
	`SPE=MM ; awk -vOFS="" '/^>/{split($1,split1,"_") ; split(split1[4],split2,":") ; split(split1[1],split3,">") ; print ">" , split1[3] , "_" , split2[1] , "|" , "'"$SPE"'" , "_" , split3[2] , "|" , split1[2] , "|" , "allele_1" ; next}{print}' MM_synonymous_MD_HET-acc2-genes.temp > MM_synonymous_MD_HET-acc2-genes.fasta`
	`cat MM_synonymous_MD_HET-acc1-genes.fasta MM_synonymous_MD_HET-acc2-genes.fasta > MM_synonymous_MD_HET-genes.fasta`
	`rm *-acc1-* ; rm *-acc2-*`
	#grep ">" *-genes.fasta | cut -d"|" -f1 | cut -d">" -f2 | sort -n | uniq | wc -l: 11326
	#grep -e '^$' *-genes.fasta -B1: ZERO BLANK LINE
	`SPE=PM ; awk -vOFS="" '/^>/{split($1,split1,"_") ; split(split1[3],split2,":") ; split(split1[1],split3,">") ; print ">" , split2[1] , "|" , "'"$SPE"'" , "_" , split3[2] , "|" , split1[2] , "|" , "allele_1" ; next}{print}' PM_synonymous_MD_HET-acc1-genes.temp > PM_synonymous_MD_HET-acc1-genes.fasta`
	`SPE=PM ; awk -vOFS="" '/^>/{split($1,split1,"_") ; split(split1[3],split2,":") ; split(split1[1],split3,">") ; print ">" , split2[1] , "|" , "'"$SPE"'" , "_" , split3[2] , "|" , split1[2] , "|" , "allele_1" ; next}{print}' PM_synonymous_MD_HET-acc2-genes.temp > PM_synonymous_MD_HET-acc2-genes.fasta`
	`cat PM_synonymous_MD_HET-acc1-genes.fasta PM_synonymous_MD_HET-acc2-genes.fasta > PM_synonymous_MD_HET-genes.fasta`
	`rm *-acc1-* ; rm *-acc2-*`
	#grep ">" *-genes.fasta | cut -d"|" -f1 | cut -d">" -f2 | sort -n | uniq | wc -l: 8814
	#grep -e '^$' *-genes.fasta -B1: ZERO BLANK LINE


#	step13_____toBootSFS
	Coordinates: within chromosomes.

	Sample among transcripts with replacement.
	`cp ../RMD/LA_MD_Filtered_Syn.4P.RMD.map LA_MD_Filtered_Syn.4P.boot.map ; cp ../HET/LA_MD_Filtered_Syn.4P.CROP.HET.ped LA_MD_Filtered_Syn.4P.boot.ped ; cp ../HET/LA_MD_Filtered_Syn.4P.CROP.HET.nosex LA_MD_Filtered_Syn.4P.boot.nosex`
	`cut -f2 LA_MD_Filtered_Syn.4P.boot.map | cut -d"_" -f1 | sort -n | uniq > tmp.pos`
	`Manually run: 17e_sample_replace.R / Infile: LA_MD_Filtered_Syn.4P.boot.map / Outfile: LA_MD_Filtered_Syn.4P.boot_*`
	#
	`cp ../RMD/MM_MD_Filtered_Syn.4P.RMD.map MM_MD_Filtered_Syn.4P.boot.map ; cp ../HET/MM_MD_Filtered_Syn.4P.HET.ped MM_MD_Filtered_Syn.4P.boot.ped ; cp ../HET/MM_MD_Filtered_Syn.4P.HET.nosex MM_MD_Filtered_Syn.4P.boot.nosex`
	`cut -f2 MM_MD_Filtered_Syn.4P.boot.map | cut -d"_" -f1,2 | sort -n | uniq > tmp.pos`
	`Manually run: 17e_sample_replace.R / Infile: MM_MD_Filtered_Syn.4P.boot.map / Outfile: MM_MD_Filtered_Syn.4P.boot_*`
	#
	`cp ../RMD/PM_MD_Filtered_Syn.4P.RMD.map PM_MD_Filtered_Syn.4P.boot.map ; cp ../HET/PM_MD_Filtered_Syn.4P.HET.ped PM_MD_Filtered_Syn.4P.boot.ped ; cp ../HET/PM_MD_Filtered_Syn.4P.HET.nosex PM_MD_Filtered_Syn.4P.boot.nosex`
	`cut -f2 PM_MD_Filtered_Syn.4P.boot.map | cut -d"_" -f1 | sort -n | uniq > tmp.pos`
	`Manually run: 17e_sample_replace.R / Infile: PM_MD_Filtered_Syn.4P.boot.map / Outfile: PM_MD_Filtered_Syn.4P.boot_*`

	Generate bootstrapped SFS with 4P. It requires a config file: sumstat.par.
	`Manually run: 17f_generate_boot.sh / Infile: LA_MD_Filtered_Syn.4P.boot_* / Outfile: boot_*_AFS-LA-U_CROP_WILD.dadi.txt`
	#
	`Manually run: 17f_generate_boot.sh / Infile: MM_MD_Filtered_Syn.4P.boot_* / Outfile: boot_*_AFS-MM-U_CROP_WILD.dadi.txt`
	#
	`Manually run: 17f_generate_boot.sh / Infile: PM_MD_Filtered_Syn.4P.boot_* / Outfile: boot_*_AFS-PM-U_CROP_WILD.dadi.txt`
	
	
## EXPLANATION previous version   
   
### 1. "Vcftools" VCF keep only the list of Crop and Wild, and remove missing data:  
In this step, we want to have a vcf file containing only the two populations of interest for the rest of the study. The input is a vcf file that has been filtered for paralogs.  
```bash    
/usr/local/public/bi/vcftools/bin/vcftools --vcf Input_Filtered.vcf --keep Crop_Wild_list.txt --out Name_TO_GO --recode   
grep -v '\.\/\.' Name_TO_GO.recode.vcf  > Name_TO_GO.noMiss.vcf  
```   
> ###### *vcftools v0.1.12b*    

### 2. "Plink" create bed and files formatted for 4P software:  
In this following line, we use plink to create a .bed file from the vcf input.  
```bash  
./plink --vcf Name_TO_GO.noMiss.vcf --keep-allele-order --vcf-idspace-to _ --make-bed --allow-extra-chr --out Name_TO_GO.noMiss.4P
```  
> ###### *plink v1.90p*    

### 3. Rename the .fam file before to create the ped files for 4P software:  
*WARNING:* It is really important to rename the file.fam with the population information for each individuals, as it will be used in the future step to seperate the crop population to the wild. Eventually if you use two populations wild, you can rename them wild1 and wild2 but keep it in mind for the following steps.    
```bash
./plink --bfile Name_TO_GO.noMiss.4P --recode --tab --out Name_TO_GO.noMiss. --allow-extra-chr
```  
> ###### *plink v1.90p*   

### 4. Change the .map file to get the good format'  
Once the previous step has been done, you obtain a .ped file that is paired with a .map file that needs to be slightly modified to be used for the 4P software. 
Here basically, we just need to pass from the left stage to the right one, the important part is to obtain a .map file as the right version following:  
> &nbsp; RefGenome00g005111 &nbsp; . &nbsp; 0 &nbsp; 792 &nbsp; &nbsp; &nbsp; &nbsp; 1 &nbsp; RefGenome00g005111 &nbsp; 0 &nbsp; 792  
> &nbsp; RefGenome00g005111 &nbsp; . &nbsp; 0 &nbsp; 808 &nbsp; ->  &nbsp; 1 &nbsp; RefGenome00g005111 &nbsp; 0 &nbsp; 808     
> &nbsp; RefGenome00g005111 &nbsp; . &nbsp; 0 &nbsp; 812 &nbsp; &nbsp; &nbsp; &nbsp; 1 &nbsp; RefGenome00g005111 &nbsp; 0 &nbsp; 8  
> &nbsp;  ...   
	
```bash
mv Name_TO_GO.noMiss..map  Name_TO_GO.noMiss.map.old  
awk -vOFS="" '{print "1","\t", $1,"_",$4,"\t", $3,"\t", $4}' Name_TO_GO.noMiss.map.old > Name_TO_GO.noMiss..map  
```  

### 5. pruning of the .ped file  
We are using a pruning method to remove linked sites, because dadi software does the assumption than all SNPs are independant and therefore not linked to each other.  
```bash
./plink --file Name_TO_GO.noMiss. --indep-pairwise 10 1 0.4 --allow-extra-chr --out LD_0.4_Pruned_files/Name_TO_GO.noMiss.
```  

### 6. filter the prune.in and prune.out  
Once the sites have been detected greater/lower than 0.4, you obtain two files, the prune.in and the prune.out. So the following step consist in keeping only the sites that were independant. Though the followind output will be a .bed file.
```bash  
./plink --file Name_TO_GO.noMiss. --extract LD_0.4_Pruned_files/Name_TO_GO.noMiss..prune.in --make-bed --out LD_0.4_Pruned_files/Name_TO_GO.noMiss.pruned  --allow-extra-chr
```   
> ###### *plink v1.90p* 

### 7. Create a Pruned ped file to get to 4p software
In this step, once again we change from the .bed file to the .ped file. 
```bash  
./plink --bfile LD_0.4_Pruned_files/Name_TO_GO.noMiss.pruned --recode --tab --out LD_0.4_Pruned_files/Name_TO_GO.noMiss.Pruned --allow-extra-chr  
```  
-> Need to get the .ped file as well with tab instead of spaces:  
```bash
sed -i 's/ /\t/g' Name_TO_GO.noMiss.Pruned.ped   
```  
> ###### *plink v1.90p; 4P software v1.0*  
And then move the files to the 4p-master folder or any folder containing your 4P software.

### 7b. Create a gene file for creating the reference    
In parallel to the previous a step, you need to prepare the reference file to obtain sites that have an ancestral state and therefore to unfold your SFS.
```bash  
grep 'mRNA' Ref_genome.gff3 > Ref_genome.gff.gene  
```  

### 8. Remake the .map file as it has to be for the 4P software !!!  FORMAT 4P MAP ENTRY    
In this step we will use the python script I created, I made one for each of the three species we were working on.   
MM: Eggplant  / PM: Pepper / LA: Tomato.    
This is mostly due to the ref_genome names that contain "_" , here again is what we got:  
> grep 'Solyc03g119960.3.1' LA.noMiss.pruned.map.temp LA.noMiss.pruned.map  
> &nbsp;LA.noMiss.pruned.map.temp: &nbsp; Solyc03g119960.3.1 &nbsp;&nbsp; 69962793    
> &nbsp;LA.noMiss.pruned.map: &nbsp; 1 &nbsp; Solyc03g119960.3.1_12 &nbsp;&nbsp;0&nbsp;&nbsp; 12    
> &nbsp;LA.noMiss.pruned.map: &nbsp; 1 &nbsp; Solyc03g119960.3.1_34 &nbsp;&nbsp;0&nbsp;&nbsp; 34    
> &nbsp;...  

What we want:   
> grep 'Solyc03g119960.3.1' LA.noMiss.pruned.map.temp LA.noMiss.pruned.ready4P.map  
> &nbsp;LA.noMiss.pruned.map.temp: &nbsp; Solyc03g119960.3.1 &nbsp;&nbsp; 69962793  
> &nbsp;LA.noMiss.pruned.map.fin: &nbsp; 1&nbsp;Solyc03g119960.3.1_12&nbsp;0&nbsp;69962805  
> &nbsp;LA.noMiss.pruned.map.fin: &nbsp; 1&nbsp;Solyc03g119960.3.1_34&nbsp;0&nbsp;69962827  
> &nbsp;...  

In all the scripts I made, you can ask for help and will see what is needed as input files and names for the command line.  
```bash   
./11_gff_map_2_NewMapPM,py --help
nohup ./11_gff_map_2_NewMapPM.py -mapin Name_TO_GO.noMiss.Pruned..map -temp RName_TO_GO.noMiss.Pruned..map.temp -fin Name_TO_GO.noMiss.Pruned.ready4P.map -gff Ref_genome.gff.gene -chrom RefGenome &   
```  
For the RefGenome as '-chrom' you only need to give the name used for the reference genome, for the tomato one gene is 'Solyc03g119960.3.1', in this case use: -chrom Solyc.  
> ###### *{Script '11_gff_map_2_NewMapXX.py'}*  

### 8b. FORMAT 4P ANC ENTRY  
In order to unfold the SFS, you need to get the .anc file. In the following script, you will first extract the outgroup individual from the original .vcf file. And then you will create the ancestral file from it.
```bash
/usr/local/public/bi/vcftools/bin/vcftools --vcf Input_Filtered.vcf --indv OUTGROUP1_ID --out OUTGROUP1_ID  --recode   
./12_VCF_to_ANC.py  -vcfin OUTGROUP1_ID.recode.vcf -tabout OUTGROUP1_ID.tab -chrom ref_genome  
```  
> ###### *plink v1.90p* 
> ###### *{Script '12_VCF_to_ANC.py'}*   
  
If you happen to have several outgroups, here is a way to make a consensus file that will be really conservative but will allow you to have a really clean unfolded file. Once again, you can run the --help to better know the input needed.    
```bash   
./13_Consensus_ANC.py --help
nohup ./13_Consensus_ANC.py -ancA OUTGROUP1_ID.tab -ancB OUTGROUP2_ID.tab -ancC OUTGROUPS.consensus.tab &
```
> ###### *{Script '13_Consensus_ANC.py'}*   
  
I developpped the script 14 then to format perfectly the ancestral file as following, here you will have to change the script according to the species (namely to the chromosomes name and separators) :    
```bash  
./14_ANC_2_NewANC_NchromMM.py --help
nohup ./14_ANC_2_NewANC_NchromMM.py -ANCin OUTGROUPS.consensus.tab -temp TEMP.Consensus_OUT.tab -fin OUTGROUPS.consensus.ready4P.tab -gff Ref_genome.gff.gene &
```  
> ###### *{Script '14_ANC_2_NewANC_NchromMM.py' & '14_ANC_2_NewANC_NchromLAandPM.py'}*  
  
  
### 9. 4P to get 2D SFS with high values of really good unfolding  
I will use the ancestral consensus file in the following step, to unfold the SFS. But before that, make sure that your .ped and .nosex file have both the order required (Crop first in the list) in order to obtain a 2D AFS with the Crop in y-axis and the Wild in the x-axis, this will be crucial for the later steps.    
For the software 4P you need to have the number of individuals you use -n X and the number of SNPs you input (an alternative is to run it with higher number and observe the limits that 4P detect).  
```bash    
./src/4P -f LA.noMiss.pruned.ped -i 0 -n 18 -s 31806 -m LA.noMiss.Pruned.ready4P.map -a OUTGROUPS.consensus.ready4P.tab 
```  
> ###### * 4P software v1.0 *  

### 10. Check visually the SFS with a script !      
> ###### *{Script '15_dadiscript.py'}   Written by Christopher Sauvage*


### Joker : INBRED LINES  
In that previous case I had inbreeding problems so I will use the following script to re-create dilpoids (two haploid for each of the fake diploid), this solution was suggested by Ryan Gutenkunst on the google group helping for dadi development.     

```bash  
./10_PED_INB_HETv2.py -ped LA.noMiss.pruned.ped -InPop CROP -Nosex 
LA.noMiss.pruned.nosex  -outputped     LA.noMiss.pruned.HET.ped  
awk -F  "\t" 'BEGIN{OFS="\t"} {print $1,$2}'  LA.noMiss.pruned.HET.ped > LA.noMiss.pruned.HET.nosex  
```  
> ###### *{Script '16_PED_INB_HETv2.py'}*  
