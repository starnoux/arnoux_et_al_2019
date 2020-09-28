# Pipeline guideline {arnoux_et_al_2019}

### In this github, you will find all the scripts that were used for the data analyses in our paper.  

The github is devided in three sections: the first to detect and clean the SNPs, the second to produce the SFS unfolded and the last one to run the scripts from the dadi software.

### A. From RNA-seq data to a clean VCF: Sequences alignment, variant calling and paralogous site cleaning
**This first guideline is giving details on the steps to align the sequences and obtain clean and reliable variant call format files (vcf).**    

*WARNING:* Use the guidline **_A_RNAseq_to_vcf.md_** that is the detailed pipeline script.
  
All the scripts are named from *01_fastqc_quality.py* to *10_subVCFfromFilter.vbeta.py* and they are split in between 3 folders, representing the different alignment and cleaning steps:  
1.  Cleaning_raw_data
2.  Mapping_And_SNPs_Calling
3.  Paralog_Cleaning  
  
   
### B. From a VCF to a clean SFS
  
This second guideline is giving details on the statistical analyses:  
*WARNING:* Use the guidline **_B_SFS_production.md_** that is the detailed pipeline script.
   
4.  From_VCF_toPED_toSFS-&-FASTA

  
   
### C.Model Inferring with dadi software    
  
This third guideline is giving indication  on how we inferred our models, this part must be adapted to your own data if you want to use them, especially as the path are written for the author's hard paths on the cluster from genotool.  
 WARNING: Use the guidline **_C_Model_Inferring_dadi.md_** that is the detailed pipeline script.  
 
5.  dadi_inference  
  
And the scripts necessary to the run the inferences are described in the following folder.    
     
6. script_dadi
   
### Addditional scripts  
   
All the scripts necessary for the DNAsp analyses (PCA and PiN/PiS analyses) were described in the folder.    
   
7.  script_PCA   

All the scripts used for the modelling using DILS are described in the folder.   

8. DILS_scripts   

