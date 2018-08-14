# B. From a clean vcf to a 2D Allele Frequency Spectrum :  
## EXPLANATION 

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
