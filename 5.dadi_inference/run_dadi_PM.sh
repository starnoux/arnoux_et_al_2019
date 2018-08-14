#!/bin/bash

MOD=$1

toScript=/home/sarnoux/work
toFold=/home/sarnoux/work/Results_July_2018_PM
toFS=/home/sarnoux/work/LD_0.4_2D_SFS
jobsscripts=/home/sarnoux/work/dadi_jobsscripts_PM


mkdir -p $jobsscripts
mkdir -p $toFold


#use at genotoul with container
#active='module load system/singularity-2.5.1'
#otherwise use python virtualenv 
#source ~/virtual_env_projects/dadi/bin/activate
active='source ~/save/virt_env_projects/dadi/bin/activate'
#use at genotoul
#deact='module unload system/singularity-2.5.1'
#otherwise use python virtualenv 
#deactivate
deact='deactivate'

#for container:
#inter='/mycontainerpath/dadi'
#for virtualenv
inter='python2'


for i in {1..50}; do
	cd ${toFold}
	nue="${MOD}_${i}"
	echo $nue
	mycmd="$inter ${toScript}/dadi_inference_models_arnouxv2_PM_2018.py -o ${MOD} -y Crop -x Wild -f ${toFS}/PMgl_order_AFS-U_CROP_WILD.dadi.txt -m ${MOD} -l -p 10 20 30"

#use the template SLURM to generate submission file prefix=nue
sed 's/MYJOBNAME/'"$nue"'/g;s/LOGFILENAME/'"$nue"'/g' ${toScript}/dadi_slurm.template >${jobsscripts}/${nue}.dadi.slurm
sed -i 's|MYWPATH|'"$toFold"'|g' ${jobsscripts}/${nue}.dadi.slurm
sed -i 's|MYNUE|'"$nue"'|g' ${jobsscripts}/${nue}.dadi.slurm
sed -i 's|MYACTIVATE|'"$active"'|g' ${jobsscripts}/${nue}.dadi.slurm
sed -i 's|MYDEACTIVATE|'"$deact"'|g' ${jobsscripts}/${nue}.dadi.slurm
sed -i 's|MYCMDFULL|'"$mycmd"'|g' ${jobsscripts}/${nue}.dadi.slurm

#submit the job
sbatch ${jobsscripts}/${nue}.dadi.slurm

#sleep 0.5

#exit 0

done

exit 0

