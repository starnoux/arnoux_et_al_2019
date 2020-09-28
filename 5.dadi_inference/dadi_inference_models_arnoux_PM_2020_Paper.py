#!/usr/local/bin/python
# -*- coding: utf-8 -*-

######################################################################################################################################################
#								IMPORT PACKAGES
######################################################################################################################################################
import os
import sys
import time
import getopt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
import numpy
from numpy import array
from math import log
sys.path.append('/home/fchristelle/.local/lib/python2.7/site-packages') #local path to python packages where dadi is installed
from scipy import stats
import dadi
import dadi_define_models_arnoux_2020_Paper


######################################################################################################################################################
#								DEFINE FUNCTIONS
######################################################################################################################################################
##### Help function
def usage():
	print("# Minimal command to run 20 replicates: for i in $(seq 1 1 20); do python dadi_inference_anneal.py -i path_to_input_name > log-${i}.txt &; done\n"+
		"# Output files are: 1) full output 2) easy-to-handle output 3) model sfs. The full output contains the inference setting, the summary statistics, starting and final logL, AIC, theta and parameter array. Note that optimization is performed with the classic annealing algorithm followed by BFGS (i.e. quasi-Newton). This requires scipy < 0.14.\n"+
	 	"# -h --help: Display this help.\n"+
		"# -v --verbose: If set, print intermediate optimization steps on stdout (best logL, current step logL, parameter array). [False]\n"+
		"# -z --masked: If set, mask the singletons in the sfs: [1,0] and [0,1]. [False]\n"+
		"# -d --folded: If set, the sfs is folded (i.e. no outgroup to polarize the data). [False]\n"+
		"# -j --projected: If set, the sfs is projected down to n1,n2 (values to be set with option '-n'). If your data have missing calls for some individuals, projecting down to a smaller sample size will increase the number of SNPs you can use. It will average over all possible re-samplings of the larger sample size data (assuming random mating). [False]\n"+
		"# -o --nameoutput: Prefix for output files. [dadi]\n"+
		"# -i --fs_file_name: Path-name to the input.\n"+
		"# -c --pop_file_name: Path-name to the file describing how individuals map to populations. This file is ignored if datatype is not 'vcf'. [None]\n"+
		"# -t --datatype: Type of input: 1) simul: ms simulation. 2) sfs: sfs data. 3) snp-dadi: snp data with dadi format; the sample size is required (values to be set with options '-j' and '-n'). 4) vcf: vcf data; it requires '-c' to be provided and it will only include sites without missing data. [sfs]\n"+
		"# -y --namepop1: Name of population 1 in the sfs (y-axis). [Pop1]\n"+
		"# -x --namepop2: Name of population 2 in the sfs (x-axis). [Pop2]\n"+
		"# -n --proj_sample: Sample size of population 1 and 2 in the sfs. Recall that for a diploid organism we get two samples from each individual. These values are ignored if '-j' is not set. [10,10]\n"+
		"# -p --grid_points: Size of the grids for extrapolation. dadi solves a partial differential equation, approximating the solution using a grid of points in population frequency space (diffusion approximation). The grid takes 3 numbers separated by a coma. Good results are obtained by setting the smallest grid size slightly larger than the largest sample size, but this can be altered depending on usage. If you are fitting a complex model, it may speed up the analysis considerably to first run an optimization at small grid sizes. They can then be refined by running another optimization with a finer grid. [nMax,nMax*2,nMax*3]\n"+
		"# -m --model_list: List the models to include in the inference. All available models are defined in 'dadi_models_2pop.py'. [SI,IM,SC,AM]\n\n\n")
	return()


##### Function that records the arguments from the command line
def takearg(argv):
	# default values
	verbose = False
	masked = False
	folded = False
	projected = False
	checkfile = False
	nameoutput = "dadi"
	pop_file_name = None
	datatype = "sfs"
	namepop1 = "Pop1"
	namepop2 = "Pop2"
	proj_sample = [10,10]
	pts_l = None
	model_list = ["SI", "IM", "SC", "AM"]

	# check values
	if len(argv) < 2:
		print("You should give, at least, the name of the sfs file !")
		sys.exit(1)
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hvzdjo:i:c:t:y:x:n:p:m:", ["help", "verbose", "masked", "folded", "projected", "nameoutput=", "fs_file_name=", "pop_file_name=", "datatype=", "namepop1=", "namepop2=", "proj_sample=", "grid_points=", "model_list="])
	except getopt.GetoptError as err:
		print(err)
		usage()
		sys.exit(2)

	# extract values from command line, or set to "True" if boolean is absent from command line.
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			usage()
			sys.exit()
		elif opt in ("-v", "--verbose"):
			verbose = True
		elif opt in ("-z", "--masked"):
			masked = True
		elif opt in ("-d", "--folded"):
			folded = True
		elif opt in ("-j", "--projected"):
			projected = True
		elif opt in ("-o", "--nameoutput"):
			nameoutput = arg
		elif opt in ("-i", "--fs_file_name"):
			fs_file_name = arg
			checkfile = True
		elif opt in ("-c", "--pop_file_name"):
			pop_file_name = arg
		elif opt in ("-t", "--datatype"):
			datatype = arg
		elif opt in ("-y", "--namepop1"):
			namepop1 = arg
		elif opt in ("-x", "--namepop2"):
			namepop2 = arg
		elif opt in ("-n", "--proj_sample"):
			proj_sample = arg.split(",")
		elif opt in ("-p", "--grid_points"):
			pts_l = arg.split(",")
		elif opt in ("-m", "--model_list"):
			model_list = arg.split(",")
		else:
			print("Option {} inconnue".format(opt))
			sys.exit(2)
	if not checkfile:
		print("You should give, at least, the name of the sfs file !")
		sys.exit(1)

	# return values
	return(verbose, masked, folded, projected, nameoutput, fs_file_name, pop_file_name, datatype, namepop1, namepop2, proj_sample, pts_l, model_list)


##### Saturation function
def lnLmax(fs,ns):
	lnL = 0
	imax = ns[0]-1
	jmax = ns[1]-1
	for i in range(0,imax):
		for j in range(0,jmax):
			k = fs[i,j]
			if k > 0:
				#fs[i,j] is Poisson distributed (Sylvain Glemin)
				lnL += -k + k*log(k) - numpy.math.lgamma(k+1)
	return(lnL)


##### Inference function
def callmodel(func, data, pts_l=None, ns=None, modeldemo=None, namepop1=None, namepop2=None,
		params=None, fixed_params=None, lower_bound=None, upper_bound=None,
		maxiterGlobal=100, Tini=50, Tfin=0, learn_rate=0.005, schedule= "cauchy", maxiterLocal=100,
		output_file="output-1", output_file_2="output-2", output_file_3="output-3",
		nameoutput=None, nameoutput2=None, verbose=True, full_output=True):

	# Extrapolate the model function.
	func_ex = dadi.Numerics.make_extrap_log_func(func)

	# Calculate the expected model SFS.
	model = func_ex(params, ns, pts_l)

	# Likelihood of the data given the model SFS.
	ll_model = dadi.Inference.ll_multinom(model, data)

	# The optimal value of theta (4*Nref*mu*L) given the model.
	theta = dadi.Inference.optimal_sfs_scaling(model, data)

	# Do the optimization.
	# Global searchs. Both 'full_output=True' and 'full_output=False' return the maxlogL.
	if optimizationstate == "anneal":

		# Perturb the parameter array before optimization by taking each parameter up to a <fold> factor of 2 up or down. 
		# This is required to generate a new initial point for the optimization between replicates.
		p0 = dadi.Misc.perturb_params(params, fold=1, lower_bound=lower_bound, upper_bound=upper_bound)

		popt = dadi.Inference.optimize_anneal(p0=p0, data=data, model_func=func_ex, pts=pts_l,
							lower_bound=lower_bound, upper_bound=upper_bound,
							maxiter=maxiterGlobal, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, schedule=schedule,
							verbose=verbose, full_output=True)
	# Local searchs (BFGS). It requires 'full_output=True' to return the maxlogL.
	else:
		popt = dadi.Inference.optimize_log(p0=params, data=data, model_func=func_ex, pts=pts_l, 
							lower_bound=lower_bound, upper_bound=upper_bound,
							maxiter=maxiterLocal, verbose=verbose, full_output=True)

	# Extract results.
	model = func_ex(popt, ns, pts_l)
	ll_opt = dadi.Inference.ll_multinom(model, data)
	theta_opt = dadi.Inference.optimal_sfs_scaling(model, data)
	AIC = 2*len(params)-2*ll_opt
	
	# Write results.
	if optimizationstate == "anneal":
		line = (str(modeldemo) + "\n" + "Starting log-likelihood: " + repr(ll_model) + "\n" + "Starting theta: " + repr(theta) + "\n" + "Starting parameters: " + ';'.join(map(str, p0)) + "\n\n" "Optimization: " + repr(optimizationstate) + "\n" + "AIC: " + repr(AIC) + "\n" + "Optimized log-likelihood: " + repr(ll_opt) + "\n" + "Optimized theta: " + repr(theta_opt) + "\n" + "Optimized parameters: " + ';'.join(map(str, popt)) + "\n")
		output_file.write(line)
	else:
		line = ("\n" "Optimization: " + repr(optimizationstate) + "\n" + "AIC: " + repr(AIC) + "\n" + "Optimized log-likelihood: " + repr(ll_opt) + "\n" + "Optimized theta: " + repr(theta_opt) + "\n" + "Optimized parameters: " + ';'.join(map(str, popt)) + "\n")
		output_file.write(line)

	# Write results at the end of the optimization.
	if optimizationstate == "BFGS":
		#Plot SFS: data, model, residuals.
		fig = pylab.figure(4,figsize=(10,8))
		fig.clear()
		dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, pop_ids =(namepop1,namepop2), residual='Anscombe', show=False)
		fig.savefig(nameoutput + "_" + modeldemo + ".png")

		#Extract best logL and AIC for each run.
		line = (str(nameoutput2) + "\t" + repr(ll_opt)+ "\t" + repr(AIC) + "\n")
		output_file_2.write(line)

		#Extract model SFS in a file.
		model_sfs = dadi.Inference.optimally_scaled_sfs(model, data)
		for i in range(1,len(model_sfs)-1):
			output_file_3.write(str(model_sfs[i]) + '\n')

 	done=True
	return(done, popt)


######################################################################################################################################################
#								LOADINGS, SETTINGS, STATISTICS
######################################################################################################################################################
##### Load parameters from the command line
verbose, masked, folded, projected, nameoutput, fs_file_name, pop_file_name, datatype, namepop1, namepop2, proj_sample, pts_l, model_list = takearg(sys.argv)

##### Load data
datastate = "not_folded"
if datatype == "simul":
	if folded == True:
		print(("\nnPolarization is knwown in ms simulations !\n"))
		datastate = "folded"
	else:
		data = dadi.Spectrum.from_ms_file(fs_file_name, average=False)
		ns = data.sample_sizes
elif datatype == "sfs":
	data = dadi.Spectrum.from_file(fs_file_name)
	if folded == True:
		data = data.fold()
		datastate = "folded"
	ns = data.sample_sizes
elif datatype == "snp-dadi":
	dd = dadi.Misc.make_data_dict(fs_file_name)
	if projected == True:
		if folded == True:
			data = dadi.Spectrum.from_data_dict(dd, pop_id=[namepop1,namepop2], projections=[int(proj_sample[0]),int(proj_sample[1])], polarized=False)
			datastate = "folded"
		else:
			data = dadi.Spectrum.from_data_dict(dd, pop_id=[namepop1,namepop2], projections=[int(proj_sample[0]),int(proj_sample[1])], polarized=True)
	else:
		print(("\nYou should set '-p' and specify the sample sizes with '-n' !\n"))
	ns = data.sample_sizes
else:
	dd = dadi.Misc.make_data_dict_vcf(fs_file_name, pop_file_name)
	if projected == True:
		if folded == True:
			data = dadi.Spectrum.from_data_dict(dd, pop_id=[namepop1,namepop2], projections=[int(proj_sample[0]),int(proj_sample[1])], polarized=False)
			datastate = "folded"
		else:
			print(("\nPolarized data is not implemented in the 'import from VCF' module !\n"))
	else:
		print(("\nYou should set '-p' and specify the sample sizes with '-n' !\n"))
	ns = data.sample_sizes

##### Set singleton masking
datastate2 = "not_masked"
if masked == True:
	data.mask[1,0] = True
	data.mask[0,1] = True
	datastate2 = "masked"

##### Set down-projection
datastate3 = "not_projdown"
if projected == True:
	data = data.project([int(proj_sample[0]),int(proj_sample[1])])
        ns = data.sample_sizes
	datastate3 = "projdown"

##### Set grid sizes
if pts_l != None:
	for i in range(len(pts_l)):
		pts_l[i] = int(pts_l[i])
else:
	pts_l = [max(ns),max(ns)*2,max(ns)*3]

##### Set optimization algorithm
opt_list = ["anneal","BFGS"]

##### Set prefix
nameoutput = nameoutput + "_" + repr(time.localtime()[0]) + repr(time.localtime()[1]) + repr(time.localtime()[2]) + repr(time.localtime()[3]) + repr(time.localtime()[4]) + repr(time.localtime()[5])

##### Set outputs
os.mkdir(nameoutput)
output_file = open((nameoutput + "/" + nameoutput + "-full.txt"), "w")
output_file_2 = open((nameoutput + "/" + nameoutput + "-easy.txt"), "w")
output_file_3 = open((nameoutput + "/" + nameoutput + "-modelsfs.txt"), "w")


##### Calculate statistics from the data
# Likelihood of the saturated model
ll_max = lnLmax(data, ns)

# Total number of segregating sites S
S = data.S()

# Wright's FST by the method of Weir and Cockerham.
Fst = round(data.Fst(),3)

# One-pop statistics
data1 = data.marginalize([1])
data2 = data.marginalize([0])
S1 = data1.S()
S2 = data2.S()

# Watterson's theta
thetaW1 = round(data1.Watterson_theta(),0)
thetaW2 = round(data2.Watterson_theta(),0)

# Expected heterozygosity pi assuming random mating
pi1 = round(data1.pi(),0)
pi2 = round(data2.pi(),0)

# Tajima's D
D1 = round(data1.Tajima_D(),3)
D2 = round(data2.Tajima_D(),3)

# Write-down inference setting
line = ("Model(s): " + repr(model_list) + "; Optimization(s): " + repr(opt_list) + "\n" + "Data type: " + repr(datatype) + "; Pop1: " + repr(namepop1) + " of sample size: " + repr(ns[0]) + "; Pop2: " + repr(namepop2) + " of sample size: " + repr(ns[1]) + "\n" + "SFS: " + repr(datastate) + "; Singletons: " + repr(datastate2) + "; Down-projection: " + repr(datastate3) + "; Grid size: " + repr(pts_l) + "\n" + "Number of seg. sites S1: " + repr(S1) + "; Watterson's Theta1: " + repr(thetaW1) + "; Expected heterozygosity pi1: " + repr(pi1) + "; Tajima's D1: " + repr(D1) + "\n" + "Number of seg. sites S2: " + repr(S2) + "; Watterson's Theta2: " + repr(thetaW2) + "; Expected heterozygosity pi2: " + repr(pi2) + "; Tajima's D2: " + repr(D2) + "\n" + "Weir & Cockerham's Fst: " + repr(Fst) + "; Total number of seg. sites S: " + repr(S) + "\n" + "Likelihood saturated model: " + repr(ll_max) + "\n\n\n")
output_file.write(line)


######################################################################################################################################################
#								PERFORM THE INFERENCE
######################################################################################################################################################
##### Set optimization values
maxiterGlobal=50 #maximum global search iterations
Tini=50 #initial temperature
Tfin=0 #final temperature
learn_rate=0.005 #scale constant for adjusting guesses
schedule= "cauchy" #annealing schedule
maxiterLocal=20 #maximum local search iterations

##### Set bounds and starting values
# If fits often push the bounds of the parameter space, this indicates 1) that bounds are too conservative. 2) that the model is misspecified.
# It is often useful to optimize only a subset of model parameters, and so to fix others with 'fixed_params'.
nuC_max = 100 ; nuC_min = 1e-4 ; nuC_start=1 #current size of Crop
nuW_max = 100 ; nuW_min = 1e-4 ; nuW_start=1 #current size of Wild
nuCe_max = 100 ; nuCe_min = 1e-4 ; nuCe_start=1 #size of Crop after exponential change
nuWe_max = 100 ; nuWe_min = 1e-4 ; nuWe_start=1 #size of Wild after exponential change
nuCb_max = 1 ; nuCb_min = 1e-4 ; nuCb_start=0.1 #size of Crop after bottleneck
Ts_max = 20 ; Ts_min = 1e-10 ; Ts_start=0.5 #duration of the first period (after split)
Tb_max = 20 ; Tb_min = 1e-10 ; Tb_start=0.5 #duration of the intermediate period
Te_max = 20 ; Te_min = 1e-10 ; Te_start=0.5 #duration of the last period
mCW_max = 10 ; mCW_min = 1e-10 ; mCW_start=1 #migration from Wild to Crop
mWC_max = 10 ; mWC_min = 1e-10 ; mWC_start=1 #migration from Crop to Wild
mCW2_max = 10 ; mCW2_min = 1e-10 ; mCW2_start=1 #migration from Wild to Crop during bottleneck in 2-migration models
mWC2_max = 10 ; mWC2_min = 1e-10 ; mWC2_start=1 #migration from Crop to Wild during bottleneck in 2-migration models
O_max = 1 ; O_min = 1e-10 ; O_start=0.8 #fraction of SNPs accurately oriented relative to the outroup

##### Inference for each model in the list
for namemodel in model_list:
	print(namemodel)
	time.sleep(1.0)

	# SI_C: Ts, O
	if namemodel == "SI_C":
		func = dadi_define_models_arnoux_2020_Paper.SI_C

		for optimizationstate in opt_list:
			print(optimizationstate)

			if optimizationstate == "anneal":
				params = [Ts_start, O_start]
			else:
				params = [popt[0], popt[1]]

			upper_bound = [Ts_max, O_max]
			lower_bound = [Ts_min, O_min]

			done, popt = callmodel(func=func, data=data, pts_l=pts_l, ns=ns, modeldemo=namemodel, namepop1=namepop1, namepop2=namepop2,
					 params=params, fixed_params=None, lower_bound=lower_bound, upper_bound=upper_bound,
					 maxiterGlobal=maxiterGlobal, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, schedule=schedule, maxiterLocal=maxiterLocal, 
					 output_file=output_file, output_file_2=output_file_2, output_file_3=output_file_3, 
					 nameoutput=nameoutput + "/" + nameoutput, nameoutput2=nameoutput, verbose=verbose, full_output=True)

		if done: print(("\n" + namemodel + " : done\n"))

	# IM_C: mCW, mWC, Ts, O
	if namemodel == "IM_C":
		func = dadi_define_models_arnoux_2020_Paper.IM_C

		for optimizationstate in opt_list:
			print(optimizationstate)

			if optimizationstate == "anneal":
				params = [mCW_start, mWC_start, Ts_start, O_start]
			else :
				params = [popt[0], popt[1], popt[2], popt[3]]
		
			upper_bound = [mCW_max, mWC_max, Ts_max, O_max]
			lower_bound = [mCW_min, mWC_min, Ts_min, O_min]

			done, popt = callmodel(func=func, data=data, pts_l=pts_l, ns=ns, modeldemo=namemodel, namepop1=namepop1, namepop2=namepop2,
					 params=params, fixed_params=None, lower_bound=lower_bound, upper_bound=upper_bound,
					 maxiterGlobal=maxiterGlobal, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, schedule=schedule, maxiterLocal=maxiterLocal, 
					 output_file=output_file, output_file_2=output_file_2, output_file_3=output_file_3, 
					 nameoutput=nameoutput + "/" + nameoutput, nameoutput2=nameoutput, verbose=verbose, full_output=True)

		if done: print(("\n" + namemodel + " : done\n"))

	# IM_BcCw: nuCb, mCW, mWC, Ts, O
	if namemodel == "IM_BcCw":
		func = dadi_define_models_arnoux_2020_Paper.IM_BcCw

		for optimizationstate in opt_list:
			print(optimizationstate)

			if optimizationstate == "anneal":
				params = [nuCb_start, mCW_start, mWC_start, Ts_start, O_start]
			else :
				params = [popt[0], popt[1], popt[2], popt[3], popt[4]]
		
			upper_bound = [nuCb_max, mCW_max, mWC_max, Ts_max, O_max]
			lower_bound = [nuCb_min, mCW_min, mWC_min, Ts_min, O_min]

			done, popt = callmodel(func=func, data=data, pts_l=pts_l, ns=ns, modeldemo=namemodel, namepop1=namepop1, namepop2=namepop2,
					 params=params, fixed_params=None, lower_bound=lower_bound, upper_bound=upper_bound,
					 maxiterGlobal=maxiterGlobal, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, schedule=schedule, maxiterLocal=maxiterLocal, 
					 output_file=output_file, output_file_2=output_file_2, output_file_3=output_file_3, 
					 nameoutput=nameoutput + "/" + nameoutput, nameoutput2=nameoutput, verbose=verbose, full_output=True)

		if done: print(("\n" + namemodel + " : done\n"))

	# IM_BcCw_E: nuCb, nuC, nuW, mCW, mWC, Ts, Te, O
	if namemodel == "IM_BcCw_E":
		func = dadi_define_models_arnoux_2020_Paper.IM_BcCw_E
		
		for optimizationstate in opt_list:
			print(optimizationstate)

			if optimizationstate == "anneal":
				params = [nuCb_start, nuC_start, nuW_start, mCW_start, mWC_start, Ts_start, Te_start, O_start]
			else :
				params = [popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7]]
		
			upper_bound = [nuCb_max, nuC_max, nuW_max, mCW_max, mWC_max, Ts_max, Te_max, O_max]
			lower_bound = [nuCb_min, nuC_min, nuW_min, mCW_min, mWC_min, Ts_min, Te_min, O_min]

			done, popt = callmodel(func=func, data=data, pts_l=pts_l, ns=ns, modeldemo=namemodel, namepop1=namepop1, namepop2=namepop2,
					 params=params, fixed_params=None, lower_bound=lower_bound, upper_bound=upper_bound,
					 maxiterGlobal=maxiterGlobal, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, schedule=schedule, maxiterLocal=maxiterLocal, 
					 output_file=output_file, output_file_2=output_file_2, output_file_3=output_file_3, 
					 nameoutput=nameoutput + "/" + nameoutput, nameoutput2=nameoutput, verbose=verbose, full_output=True)

		if done: print(("\n" + namemodel + " : done\n"))

	# IM2_BcCw_E: nuCb, nuC, nuW, mCW, mWC, mWC2, mCW2, Ts, Te, O
	if namemodel == "IM2_BcCw_E":
		func = dadi_define_models_arnoux_2020_Paper.IM2_BcCw_E
		
		for optimizationstate in opt_list:
			print(optimizationstate)

			if optimizationstate == "anneal":
				params = [nuCb_start, nuC_start, nuW_start, mCW_start, mWC_start, mCW2_start, mWC2_start, Ts_start, Te_start, O_start]
			else :
				params = [popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9]]
		
			upper_bound = [nuCb_max, nuC_max, nuW_max, mCW_max, mWC_max, mCW2_max, mWC2_max, Ts_max, Te_max, O_max]
			lower_bound = [nuCb_min, nuC_min, nuW_min, mCW_min, mWC_min, mCW2_min, mWC2_min, Ts_min, Te_min, O_min]

			done, popt = callmodel(func=func, data=data, pts_l=pts_l, ns=ns, modeldemo=namemodel, namepop1=namepop1, namepop2=namepop2,
					 params=params, fixed_params=None, lower_bound=lower_bound, upper_bound=upper_bound,
					 maxiterGlobal=maxiterGlobal, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, schedule=schedule, maxiterLocal=maxiterLocal, 
					 output_file=output_file, output_file_2=output_file_2, output_file_3=output_file_3, 
					 nameoutput=nameoutput + "/" + nameoutput, nameoutput2=nameoutput, verbose=verbose, full_output=True)

		if done: print(("\n" + namemodel + " : done\n"))

	# IM_C_BcCw_E: nuCb, nuC, nuW, mCW, mWC, Ts, Tb, Te, O
	if namemodel == "IM_C_BcCw_E":
		func = dadi_define_models_arnoux_2020_Paper.IM_C_BcCw_E

		for optimizationstate in opt_list:
			print(optimizationstate)

			if optimizationstate == "anneal":
				params = [nuCb_start, nuC_start, nuW_start, mCW_start, mWC_start, Ts_start, Tb_start, Te_start, O_start]
			else :
				params = [popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8]]
		
			upper_bound = [nuCb_max, nuC_max, nuW_max, mCW_max, mWC_max, Ts_max, Tb_max, Te_max, O_max]
			lower_bound = [nuCb_min, nuC_min, nuW_min, mCW_min, mWC_min, Ts_min, Tb_min, Te_min, O_min]

			done, popt = callmodel(func=func, data=data, pts_l=pts_l, ns=ns, modeldemo=namemodel, namepop1=namepop1, namepop2=namepop2,
					 params=params, fixed_params=None, lower_bound=lower_bound, upper_bound=upper_bound,
					 maxiterGlobal=maxiterGlobal, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, schedule=schedule, maxiterLocal=maxiterLocal, 
					 output_file=output_file, output_file_2=output_file_2, output_file_3=output_file_3, 
					 nameoutput=nameoutput + "/" + nameoutput, nameoutput2=nameoutput, verbose=verbose, full_output=True)

		if done: print(("\n" + namemodel + " : done\n"))

	# IM2_C_BcCw_E: nuCb, nuC, nuW, mCW, mWC, mWC2, mCW2, Ts, Tb, Te, O
	if namemodel == "IM2_C_BcCw_E":
		func = dadi_define_models_arnoux_2020_Paper.IM2_C_BcCw_E

		for optimizationstate in opt_list:
			print(optimizationstate)

			if optimizationstate == "anneal":
				params = [nuCb_start, nuC_start, nuW_start, mCW_start, mWC_start, mCW2_start, mWC2_start, Ts_start, Tb_start, Te_start, O_start]
			else :
				params = [popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10]]
		
			upper_bound = [nuCb_max, nuC_max, nuW_max, mCW_max, mWC_max, mCW2_max, mWC2_max, Ts_max, Tb_max, Te_max, O_max]
			lower_bound = [nuCb_min, nuC_min, nuW_min, mCW_min, mWC_min, mCW2_min, mWC2_min, Ts_min, Tb_min, Te_min, O_min]

			done, popt = callmodel(func=func, data=data, pts_l=pts_l, ns=ns, modeldemo=namemodel, namepop1=namepop1, namepop2=namepop2,
					 params=params, fixed_params=None, lower_bound=lower_bound, upper_bound=upper_bound,
					 maxiterGlobal=maxiterGlobal, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, schedule=schedule, maxiterLocal=maxiterLocal, 
					 output_file=output_file, output_file_2=output_file_2, output_file_3=output_file_3, 
					 nameoutput=nameoutput + "/" + nameoutput, nameoutput2=nameoutput, verbose=verbose, full_output=True)

		if done: print(("\n" + namemodel + " : done\n"))

	# IM_E: nuC, nuW, mCW, mWC, Ts, O
	if namemodel == "IM_E":
		func = dadi_define_models_arnoux_2020_Paper.IM_E

		for optimizationstate in opt_list:
			print(optimizationstate)

			if optimizationstate == "anneal":
					params = [nuC_start, nuW_start, mCW_start, mWC_start, Ts_start, O_start]
			else :
					params = [popt[0], popt[1], popt[2], popt[3], popt[4], popt[5]]

			upper_bound = [nuC_max, nuW_max, mCW_max, mWC_max, Ts_max, O_max]
			lower_bound = [nuC_min, nuW_min, mCW_min, mWC_min, Ts_min, O_min]

			done, popt = callmodel(func=func, data=data, pts_l=pts_l, ns=ns, modeldemo=namemodel, namepop1=namepop1, namepop2=namepop2,
					 params=params, fixed_params=None, lower_bound=lower_bound, upper_bound=upper_bound,
					 maxiterGlobal=maxiterGlobal, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, schedule=schedule, maxiterLocal=maxiterLocal, 
					 output_file=output_file, output_file_2=output_file_2, output_file_3=output_file_3, 
					 nameoutput=nameoutput + "/" + nameoutput, nameoutput2=nameoutput, verbose=verbose, full_output=True)

		if done: print(("\n" + namemodel + " : done\n"))	 

	# IM_E_E: nuCe, nuC, nuWe, nuW, mCW, mWC, Ts, Te, O
	if namemodel == "IM_E_E":
		func = dadi_define_models_arnoux_2020_Paper.IM_E_E

		for optimizationstate in opt_list:
			print(optimizationstate)

			if optimizationstate == "anneal":
				params = [nuCe_start, nuC_start, nuWe_start, nuW_start, mCW_start, mWC_start, Ts_start, Te_start, O_start]
			else :
				params = [popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8]]
		
			upper_bound = [nuCe_max, nuC_max, nuWe_max, nuW_max, mCW_max, mWC_max, Ts_max, Te_max, O_max]
			lower_bound = [nuCe_min, nuC_min, nuWe_min, nuW_min, mCW_min, mWC_min, Ts_min, Te_min, O_min]

			done, popt = callmodel(func=func, data=data, pts_l=pts_l, ns=ns, modeldemo=namemodel, namepop1=namepop1, namepop2=namepop2,
					 params=params, fixed_params=None, lower_bound=lower_bound, upper_bound=upper_bound,
					 maxiterGlobal=maxiterGlobal, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, schedule=schedule, maxiterLocal=maxiterLocal, 
					 output_file=output_file, output_file_2=output_file_2, output_file_3=output_file_3, 
					 nameoutput=nameoutput + "/" + nameoutput, nameoutput2=nameoutput, verbose=verbose, full_output=True)

		if done: print(("\n" + namemodel + " : done\n"))

	# IM2_E_E: nuCe, nuC, nuWe, nuW, mCW, mWC, mCW2, mWC2, Ts, Te, O
	if namemodel == "IM2_E_E":
		func = dadi_define_models_arnoux_2020_Paper.IM2_E_E

		for optimizationstate in opt_list:
			print(optimizationstate)

			if optimizationstate == "anneal":
				params = [nuCe_start, nuC_start, nuWe_start, nuW_start, mCW_start, mWC_start, mCW2_start, mWC2_start, Ts_start, Te_start, O_start]
			else :
				params = [popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10]]
		
			upper_bound = [nuCe_max, nuC_max, nuWe_max, nuW_max, mCW_max, mWC_max, mCW2_max, mWC2_max, Ts_max, Te_max, O_max]
			lower_bound = [nuCe_min, nuC_min, nuWe_min, nuW_min, mCW_min, mWC_min, mCW2_min, mWC2_min, Ts_min, Te_min, O_min]

			done, popt = callmodel(func=func, data=data, pts_l=pts_l, ns=ns, modeldemo=namemodel, namepop1=namepop1, namepop2=namepop2,
					 params=params, fixed_params=None, lower_bound=lower_bound, upper_bound=upper_bound,
					 maxiterGlobal=maxiterGlobal, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, schedule=schedule, maxiterLocal=maxiterLocal, 
					 output_file=output_file, output_file_2=output_file_2, output_file_3=output_file_3, 
					 nameoutput=nameoutput + "/" + nameoutput, nameoutput2=nameoutput, verbose=verbose, full_output=True)

		if done: print(("\n" + namemodel + " : done\n"))

	# IM_C_E_E: nuCe, nuC, nuWe, nuW, mCW, mWC, Ts, Tb, Te, O
	if namemodel == "IM_C_E_E":
		func = dadi_define_models_arnoux_2020_Paper.IM_C_E_E

		for optimizationstate in opt_list:
			print(optimizationstate)

			if optimizationstate == "anneal":
				params = [nuCe_start, nuC_start, nuWe_start, nuW_start, mCW_start, mWC_start, Ts_start, Tb_start, Te_start, O_start]
			else :
				params = [popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9]]
		
			upper_bound = [nuCe_max, nuC_max, nuWe_max, nuW_max, mCW_max, mWC_max, Ts_max, Tb_max, Te_max, O_max]
			lower_bound = [nuCe_min, nuC_min, nuWe_min, nuW_min, mCW_min, mWC_min, Ts_min, Tb_min, Te_min, O_min]

			done, popt = callmodel(func=func, data=data, pts_l=pts_l, ns=ns, modeldemo=namemodel, namepop1=namepop1, namepop2=namepop2,
					 params=params, fixed_params=None, lower_bound=lower_bound, upper_bound=upper_bound,
					 maxiterGlobal=maxiterGlobal, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, schedule=schedule, maxiterLocal=maxiterLocal, 
					 output_file=output_file, output_file_2=output_file_2, output_file_3=output_file_3, 
					 nameoutput=nameoutput + "/" + nameoutput, nameoutput2=nameoutput, verbose=verbose, full_output=True)

		if done: print(("\n" + namemodel + " : done\n"))

	# IM2_C_E_E: nuCe, nuC, nuWe, nuW, mCW, mWC, mCW2, mWC2, Ts, Tb, Te, O
	if namemodel == "IM2_C_E_E":
		func = dadi_define_models_arnoux_2020_Paper.IM2_C_E_E

		for optimizationstate in opt_list:
			print(optimizationstate)

			if optimizationstate == "anneal":
				params = [nuCe_start, nuC_start, nuWe_start, nuW_start, mCW_start, mWC_start, mCW2_start, mWC2_start, Ts_start, Tb_start, Te_start, O_start]
			else :
				params = [popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10], popt[11]]
		
			upper_bound = [nuCe_max, nuC_max, nuWe_max, nuW_max, mCW_max, mWC_max, mCW2_max, mWC2_max, Ts_max, Tb_max, Te_max, O_max]
			lower_bound = [nuCe_min, nuC_min, nuWe_min, nuW_min, mCW_min, mWC_min, mCW2_min, mWC2_min, Ts_min, Tb_min, Te_min, O_min]

			done, popt = callmodel(func=func, data=data, pts_l=pts_l, ns=ns, modeldemo=namemodel, namepop1=namepop1, namepop2=namepop2,
					 params=params, fixed_params=None, lower_bound=lower_bound, upper_bound=upper_bound,
					 maxiterGlobal=maxiterGlobal, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, schedule=schedule, maxiterLocal=maxiterLocal, 
					 output_file=output_file, output_file_2=output_file_2, output_file_3=output_file_3, 
					 nameoutput=nameoutput + "/" + nameoutput, nameoutput2=nameoutput, verbose=verbose, full_output=True)

		if done: print(("\n" + namemodel + " : done\n"))

output_file.close()
output_file_2.close()
output_file_3.close()
