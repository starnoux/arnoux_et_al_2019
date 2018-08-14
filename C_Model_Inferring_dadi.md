# C. Inferring evolutionnary model on the 2D Allele Frequency Spectrum :  
## EXPLANATION 

### 1. dadi installation // Virtual environment:
You'll need to install the modified version of dadi software. You can find the compressed folder in the 5.dadi_inference. 

```bash   
tar -xzvf dadi-1.7.0_modif.tar.gz  
```  
It would be greatly advised to use a virtual environement to have no problem with new versions of softwares.   
To use the virtual environment please have a look to the README that was developped by *Jacques Lagnel*.   
> ###### *dadi v1.7.0* 
> ###### *{Script 'README.virtualenv.txt'} Written by Jacques Lagnel*  

  
### 2. Use of slurm to run on a cluster the model inferrences.  
Here we explain the runninf script with slurm as most of the cluster are running with slurm nowadays.  
```bash  
while read line; do  
./run_dadi_MM.sh ${line}   
done<./toModels.txt  
```  
**Warning** You will need the following files:
- toModels.txt #List the models that have to be inferred
- dadi_slurm.template #Used as template to create jobs for each inferrence
- run_dadi_{SPECIES}.sh #To run for each species ** the proper command with the 2D_sfs
- dadi_define_models_arnoux18.py  #To define the models we use (See the Figure1 in our paper)
- dadi_inference_models_arnouxv2_{SPECIES}_2018.py #To give the parameter boundaries for dadi to infer the models  
  
**Warning** you need to keep the y = Crop and x = Wild, as the dadi inference are with asymetric models and that you want to keep sure you are interpreting your results in the right way.  
  
### 3. Detecting the best models and selecting the best parameters  
In this part we will recover all data and create a summary file giving all the necessary informations from inferred models, to compare them between each other.  

```bash  
for j in `ls -d Results_*`; do
  echo ${j}
  cd ${j}
  chmod a+wrx *
  for i in `ls -d *_2018_*`; do
    cd ${i}
    myFILE=`ls | grep "_2.txt"`
    cat ${myFILE}>> ~/work/${j}_sumstat.txt
    echo ${i}
    cd ..
  done
  cd ..
done
```  
In the following script, we will use the previously produced summary file in order to create lists of, so called, "Best_AIC" models. And once the models are retrieved, we will check on all the parameters that were produced with the dadi software.  
    
> ###### *{Script '17_PLot_Loglikelihood.R'}*     
   
Then we will retrieve all files from the best model list "Best_AIC_Models_{species}.txt" and create a folder "BestAIC_Models_{species}{newname}" containing all the output files related to the best likelihood model of each scenario.  

```bash  
for j in `ls Best_AIC_Models_{species}*.txt`; do  
  echo $j  
  Species=`echo {species}`  
  newname=`echo $(echo $j | sed 's/Best_AIC_Models_{species}_//' | sed 's/.txt//')`  
  #num=`echo $(echo $j | sed 's/Best_AIC_Models_{species}_/' | sed 's/.txt//')`  
  while read line; do cp Results_${newname}_${Species}/${line}/${line}*.png Best_AIC_Models_{species}_${newname}/; done< ${j}  
  while read line; do cp Results_${newname}_${Species}/${line}/${line}.txt Best_AIC_Models_{species}_${newname}/; done< ${j}   
  while read line; do cp Results_${newname}_${Species}/${line}/${line}_2.txt Best_AIC_Models_{species}_${newname}/; done< ${j}  
  while read line; do tail -5 Results_${newname}_${Species}/${line}/${line}.txt > Best_AIC_Models_{species}_${newname}/${line}.opt ; done< ${j}   
  cd Best_AIC_Models_{species}_${newname}/  
  more *.opt >> ../${Species}_${newname}_opt_sumstat.txt  
  cd ..  
done  
```


### 4. Bootstrap production
" The general idea is to resample the data several times, and to find the best model (with 1,000 replicates) for each bootstrapped dataset. This gives you a distribution of parameter values, from which you get the confidence intervals. Since version 1.7.0, a new set of functions has been added (Godambe.py) which includes methods for computationally-efficient parameter uncertainty estimation; it avoids doing the fit of each bootstrapped dataset.

Importantly, this method (Fisher Information Matrix) relies on a Gaussian distribution of the parameter uncertainties (which breaks down is the uncertainties are too large); and independant data. If the confidence intervals are larger than half of the parameter value (large uncertainties), the method does not apply, and one should do "conventional" bootstrapping. Moreover, the method struggles if the estimates are very low (but in that case, a nested model should work better !). "
*AJ Coffman, P Hsieh, S Gravel, RN Gutenkunst "Computationally efficient composite likelihood statistics for demographic inference" Molecular Biology and Evolution 33:591-593 (2016)*

**Warning** Any time you want to use dadi, remember to activate your virtual environment.   
  
```bash    
python 18_dadi_estim_params_arnoux_{SPECIES}_v2.py  
```  
> ###### *{Script '18_dadi_estim_params_arnoux_{SPECIES}_v2.py'}*  
