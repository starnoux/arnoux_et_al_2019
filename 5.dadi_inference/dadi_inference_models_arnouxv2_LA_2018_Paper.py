#!/usr/local/bin/python
# -*- coding: utf-8 -*-

#########################################
########## Library importation ##########
######################################### 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import getopt
import pylab
import time
from scipy import stats

import dadi
import dadi_define_models_arnoux18_Paper

#sys.path.append('python')
import numpy
from numpy import array
import random


###################################
########## Help function ########## 
###################################
def usage():
	""" Function for help """
	print("# This script allow you to test different demographic models on your genomic data\n"+
		  "# and will give you which one is the best fitted.\n\n"+
		  "# This is an exemple of the most complete command line :\n"+
		  "# -o pathoutput -y population1 -x population2 -p 10,20,30 -f pathfsfile -m SI_C,IM_C -l -a -h -v\n\n"+
		  "# where:\n"+
		  "# -f pathfsfile\n\n"+
		  "# -h --help : Display the help you are looking at.\n"+
		  "# -v --verbose : Print steps while the code is running\n"+
		  "# -y --population1 : Take the name of the first population in the sfs (y-axis)\n"+
		  "# -x --population2 : Take the name of the second population in the sfs (x-axis)\n"+
		  "# -o --outputname : Take the path of output file.\n"+
		  "# -f --fs_file_name : Take the path of the fs file from thr parent directory.\n"+
		  "# -p --grid_points : Take 3 numbers separated by a coma, for the size of grids for extrapolation.\n"+
		  "# -m --model_list : Take names of model separated by a coma.\n"+
		  "# -z : mask the singletons.\n"+
		  "# -l : record the final parameters in the output file.\n\n\n"
		  "########################## Enjoy ###########################")
	return()
		  
		  
#######################################
########## Argument function ##########
#######################################
def takearg(argv):
	""" Function which record arguments from the command line."""
	# Default values
	masked = False # Freq 0,1 and 1,0 masked if masked = 1
	pts_l = None  # Grids sizes for extrapolation
	outputname = "fs_2d_optlog"
	model_list = ["SI_C","IM_C","IM_BcCw","IM_BcCw_E","IM2_BcCw_E","IM_C_BcCw_E","IM2_C_BcCw_E","IM_E","IM_E_E","IM_C_E_E"]
	verbose = False
	logparam = False
	nompopC = "PonbC"
	nompopW = "PonbW"
	checkfile = False #Initilization. if True fs file needed exists, if False it doesn't

	if len(argv) < 2:
		usage()
		sys.exit(1)
	try:
		opts, args = getopt.getopt(argv[1:], "hvo:y:x:azf:p:m:l", ["help", "verbose", "outputname=", "population1=", "population2=", "masked", "fs_file_name=", "grid_points=", "model_list=", "log"])
	except getopt.GetoptError as err:
		# Print help, and exit
		print(err)
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			usage()
			sys.exit()
		elif opt in ("-v", "--verbose"):
			verbose = True
		elif opt in ("-o", "--outputname"):
			outputname = arg
		elif opt in ("-y", "--population1"):
			nompopC = arg
		elif opt in ("-x", "--population2"):
			nompopW = arg
		elif opt in ("-z", "--masked"):
			masked = True
		elif opt in ("-f", "--fs_file_name"):
			fs_file_name = arg
			checkfile = True
		elif opt in ("-p", "--grid_points"):
			pts_l = arg.split(",")
		elif opt in ("-m", "--model_list"):
			model_list = arg.split(",")
		elif opt in ("-l", "--log"):
			logparam = True
		else:
			print("Option {} unknown".format(opt))
			sys.exit(2)
	if not checkfile:
		print("You should give, at least, the name of the fs file !")
		sys.exit(1)
	return(masked, pts_l, outputname, nompopC, nompopW, fs_file_name, model_list, verbose, logparam)

######################################## PARAMETERS BEGIN
# PLOT 
#tmp = data.copy();tmp[1][0]=data.min();tmp[0][1]=data.min();vMin = int(0.8*numpy.floor(tmp.min()));vMax = int(1.2*numpy.ceil(tmp.max()))
vMin = int(1)
vMax = int(10000)
residRange = int(3)

# OPTIMIZATION "ANNEAL"
# maxiter = restricts how long the optimizer will run. You may want to set this value higher, to encourage better convergence.
maxiter=50 
# Tini = initial temperature of the chain.
Tini=50
# Tfin = final temperature of the chain.
Tfin=0
# Learn rate = decreasing rate in the probability of accepting worse solutions as it explores the solution space. 
learn_rate=0.005

# PRIORS: upper bound; lower bound; starting value
nuC_max = 12 ; nuC_min = 1e-4 ; nuC_start=1.0
nuCe_max = 20 ; nuCe_min = 1e-4 ; nuCe_start=1.0 #CF_Nov17
nuW_max = 12 ; nuW_min = 1e-4 ; nuW_start=1.0
nuWe_max = 20 ; nuWe_min = 1e-4 ; nuWe_start=1.0 #CF_Nov17
nuCb_max = 1 ; nuCb_min = 1e-6 ; nuCb_start=0.1
Ts_max = 12.0 ; Ts_min = 0.0 ; Ts_start=0.5
Tb_max = 12.0 ; Tb_min = 0 ; Tb_start=0.5
Te_max = 12.0 ; Te_min = 0.0 ; Te_start=0.5
mCW_max = 4 ; mCW_min = 0.0 ; mCW_start=1.0
mWC_max = 4 ; mWC_min = 0.0 ; mWC_start=1.0
mCW2_max = 4 ; mCW2_min = 0.0 ; mCW2_start=1.0
mWC2_max = 4 ; mWC2_min = 0.0 ; mWC2_start=1.0
O_max = 1.0 ; O_min = 0.0 ; O_start=0.8
######################################## PARAMETERS END


########################################
########## Inference function ########## 
########################################
def callmodel(func, data, output_file, output_file_2, output_file_3, modeldemo, ll_opt_dic, nbparam_dic,
	      nompopC="PonbC", nompopW="PonbW", params=None, fixed_params=None, lower_bound=None, upper_bound=None,
	      pts_l=None, ns=None,outputname=None, outputname2=None, verbose=False, maxiter=3,
	      Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, schedule= "cauchy"):

	# Make the extrapolating version of our demographic model function.
	func_ex = dadi.Numerics.make_extrap_log_func(func)
	# Calculate the model AFS.
	model = func_ex(params, ns, pts_l)
	# Likelihood of the data given the model AFS.
	ll_model = dadi.Inference.ll_multinom(model, data)
	print 'Model log-likelihood:', ll_model
	# The optimal value of theta (4*No*u) given the model.
	theta = dadi.Inference.optimal_sfs_scaling(model, data)
	print 'theta:', theta
	
	# Do the optimization. By default we assume that theta is a free parameter, since it's trivial to find given the other parameters.
	# If you want to fix theta, add a multinom=False to the call. (This is commented out by default, since it takes several minutes.)
	if optimizationstate == "anneal_hot" :
		# Perturb our parameter array before optimization. This does so by taking each
		# parameter a up to a factor of two up or down.
		p0 = dadi.Misc.perturb_params(params, fold=1, lower_bound=lower_bound, upper_bound=upper_bound)
		print 'Starting parameters', repr(p0)
		
		popt = dadi.Inference.optimize_anneal(p0, data, func_ex, pts_l, 
						      lower_bound=lower_bound,
						      upper_bound=upper_bound,
						      verbose=verbose,
						      maxiter=maxiter, Tini=Tini, Tfin=Tfin, 
						      learn_rate=learn_rate, schedule=schedule)
	else :
		popt = dadi.Inference.optimize_log(params, data, func_ex, pts_l, 
						   lower_bound=lower_bound,
						   upper_bound=upper_bound,
						   verbose=verbose,
						   maxiter=maxiter/2)
	
	# Computation of statistics
	model = func_ex(popt, ns, pts_l)
	ll_opt = dadi.Inference.ll_multinom(model, data)
	theta = dadi.Inference.optimal_sfs_scaling(model, data)
	AIC = 2*len(params)-2*ll_opt
        ll_opt_dic[modeldemo] = ll_opt
        nbparam_dic[modeldemo] = len(params)

	# Print results
	print 'Optimized parameters', repr(popt)
	print 'Optimized log-likelihood:', ll_opt
	print 'theta:', theta
	
	# Write results
	line = ("\n" + str(modeldemo) + "\n" + "Model log-likelihood: " + repr(ll_model) + "\n" "Optimization : " + repr(optimizationstate) + "\n"  "Optimized parameters: " + "\t".join([ str(param) for param in popt]) + "\n" + "Optimized log-likelihood: " + repr(ll_opt) + "\n" + "theta: " + repr(theta) + "\n" + "AIC: " + repr(AIC) + "\n")
	output_file.write(line)

	# Plot a comparison of the resulting fs with the data.
	if optimizationstate == "BFGS" :
		import pylab
		#pylab.figure()
		#dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=vMin, vmax=vMax, resid_range=residRange,
		#					    pop_ids =(nompopC,nompopW),
		#					    saveplot=True, nomplot=(outputname + "_" + modeldemo), showplot=False)
		#matplotlib.use('Agg')
		fig = plt.figure()
		dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=vMin, vmax=vMax, resid_range=residRange, pop_ids =(nompopC,nompopW))
		fig.savefig(outputname + "_" + modeldemo)
		plt.close(fig)
		line = (str(outputname2) + "\t" + repr(ll_opt)+ "\t" + repr(AIC) + "\n")
		output_file_2.write(line)
		#Extract data, model, residuals
		model = dadi.Inference.optimally_scaled_sfs(model, data)
		model = model
		masked_model, masked_data = dadi.Numerics.intersect_masks(model, data)
		max_toplot = max(masked_model.max(), masked_data.max())
		min_toplot = min(masked_model.min(), masked_data.min())
		resid = dadi.Inference.Anscombe_Poisson_residual(masked_model, masked_data, mask=min_toplot)
		line2 = (repr(data) + "\n" + repr(model) + "\n" + repr(resid) + "\n")		
		output_file_3.write(line2)
 	done=True
	return(done, ll_opt_dic, nbparam_dic, popt)



##################################################################################################
########################### Main function: actually make the inference ###########################
##################################################################################################



# Load the parameters
masked, pts_l, outputname, nompopC, nompopW, fs_file_name, model_list, verbose, logparam = takearg(sys.argv)

if pts_l != None:
	for i in range(len(pts_l)):
		pts_l[i] = int(pts_l[i])

# Load the data
data = dadi.Spectrum.from_file(fs_file_name)
ns = data.sample_sizes


# Creation of outputname and setting default params if they are not in the args
datastate = "not_masked"
opt_list = ["anneal_hot", "BFGS"]

if pts_l == None:
	pts_l = [ns[0]+10,ns[0]+20,ns[0]+30]

if masked:
	data.mask[1,0] = True
	data.mask[0,1] = True
	outputname = outputname + "_masked"
	datastate = "masked"

rn=random.randrange(0, 101, 2)
tt=time.time()
millis=int(round((tt - int(tt))*10000000))+rn;
outputname = outputname + "_" + repr(time.localtime()[0]) + "_" + repr(time.localtime()[1]) + "_" + repr(time.localtime()[2]) + "_" + repr(time.localtime()[3]) + repr(time.localtime()[4]) + repr(time.localtime()[5]) + "_" + str(millis)

poub=1
while os.path.exists(outputname) and poub<20:
	rn=random.randrange(0, 101, 2)
	tt=time.time()
	millis=int(round((tt - int(tt))*10000000)+rn+poub);
	outputname = outputname + "_" + repr(time.localtime()[0]) + "_" + repr(time.localtime()[1]) + "_" + repr(time.localtime()[2]) + "_" + repr(time.localtime()[3]) + repr(time.localtime()[4]) + repr(time.localtime()[5])+"_"+str(millis)
	poub=poub+1

# Create output dir and file
os.mkdir(outputname)
output_file = open((outputname + "/" + outputname + ".txt"), "w")


# Save the parameters
if logparam :
	line = ("Model(s) : " + repr(model_list) + "\n" + "Data state : " + repr(datastate) + "\n" + "Grid points : " + repr(pts_l) + "\n\n\n")
	output_file.write(line)

output_file_2 = open((outputname + "/" + outputname + "_2.txt"), "w")
output_file_3 = open((outputname + "/" + outputname + "_3.txt"), "w")
	
# Create dic for ll to make lrt
ll_opt_dic = {}
nbparam_dic = {}


# ML inference for each model
for namemodel in model_list:
	print namemodel
	time.sleep(1.0)

	if namemodel == "SI_C":

		# SI_C: Ts, O
		func = dadi_define_models_arnoux18_Paper.SI_C

		for optimizationstate in opt_list:
			print optimizationstate

			if optimizationstate == "anneal_hot":
				params = (Ts_start, O_start)
			else:
				params = (popt[0], popt[1])

			# The upper_bound array is for use in optimization. Occasionally the optimizer
			# will try wacky parameter values. We in particular want to exclude values with
			# very long times, as they will take a long time to evaluate.
			upper_bound = [Ts_max, O_max]
			lower_bound = [Ts_min, O_min]

			done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, output_file_2, output_file_3, namemodel, ll_opt_dic, nbparam_dic,
								  nompopC=nompopC, nompopW=nompopW, params=params, fixed_params=None, lower_bound=lower_bound, 
								  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
								  outputname= outputname + "/" + outputname,  outputname2=outputname,
								  verbose=verbose, maxiter=maxiter, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, 
								  schedule= "cauchy")
		if done: print(("\n" + namemodel + " : done\n"))


	if namemodel == "IM_C":

		# IM_C: mCW, mWC, Ts, O
		func = dadi_define_models_arnoux18_Paper.IM_C

		for optimizationstate in opt_list:
			print optimizationstate

			if optimizationstate == "anneal_hot":
				params = (mCW_start, mWC_start, Ts_start, O_start)
			else :
				params = (popt[0], popt[1], popt[2], popt[3])
		
			# The upper_bound array is for use in optimization. Occasionally the optimizer
			# will try wacky parameter values. We in particular want to exclude values with
			# very long times, as they will take a long time to evaluate.
			upper_bound = [mCW_max, mWC_max, Ts_max, O_max]
			lower_bound = [mCW_min, mWC_min, Ts_min, O_min]

			done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, output_file_2, output_file_3, namemodel, ll_opt_dic, nbparam_dic,
								  nompopC=nompopC, nompopW=nompopW, params=params, fixed_params=None, lower_bound=lower_bound, 
								  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
								  outputname= outputname + "/" + outputname,  outputname2=outputname,
								  verbose=verbose, maxiter=maxiter, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, 
								  schedule= "cauchy")
		if done: print(("\n" + namemodel + " : done\n"))


	if namemodel == "IM_BcCw":

		# IM_BcCw: nuCb, mCW, mWC, Ts, O
		func = dadi_define_models_arnoux18_Paper.IM_BcCw

		for optimizationstate in opt_list:
			print optimizationstate

			if optimizationstate == "anneal_hot":
				params = (nuCb_start, mCW_start, mWC_start, Ts_start, O_start)
			else :
				params = (popt[0], popt[1], popt[2], popt[3], popt[4])
		
			# The upper_bound array is for use in optimization. Occasionally the optimizer
			# will try wacky parameter values. We in particular want to exclude values with
			# very long times, as they will take a long time to evaluate.
			upper_bound = [nuCb_max, mCW_max, mWC_max, Ts_max, O_max]
			lower_bound = [nuCb_min, mCW_min, mWC_min, Ts_min, O_min]

			done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, output_file_2, output_file_3, namemodel, ll_opt_dic, nbparam_dic,
								  nompopC=nompopC, nompopW=nompopW, params=params, fixed_params=None, lower_bound=lower_bound, 
								  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
								  outputname= outputname + "/" + outputname,  outputname2=outputname,
								  verbose=verbose, maxiter=maxiter, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, 
								  schedule= "cauchy")
		if done: print(("\n" + namemodel + " : done\n"))


	if namemodel == "IM_BcCw_E":

		# IM_BcCw_E: nuCb, nuC, nuW, mCW, mWC, Ts, Te, O
		func = dadi_define_models_arnoux18_Paper.IM_BcCw_E
		
		for optimizationstate in opt_list:
			print optimizationstate

			if optimizationstate == "anneal_hot":
				params = (nuCb_start, nuC_start, nuW_start, mCW_start, mWC_start, Ts_start, Te_start, O_start)
			else :
				params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7])
		
			# The upper_bound array is for use in optimization. Occasionally the optimizer
			# will try wacky parameter values. We in particular want to exclude values with
			# very long times, as they will take a long time to evaluate.
			upper_bound = [nuCb_max, nuC_max, nuW_max, mCW_max, mWC_max, Ts_max, Te_max, O_max]
			lower_bound = [nuCb_min, nuC_min, nuW_min, mCW_min, mWC_min, Ts_min, Te_min, O_min]

			done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, output_file_2, output_file_3, namemodel, ll_opt_dic, nbparam_dic,
								  nompopC=nompopC, nompopW=nompopW, params=params, fixed_params=None, lower_bound=lower_bound, 
								  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
								  outputname= outputname + "/" + outputname,  outputname2=outputname,
								  verbose=verbose, maxiter=maxiter, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, 
								  schedule= "cauchy")
		if done: print(("\n" + namemodel + " : done\n"))


	if namemodel == "IM2_BcCw_E":

		# IM2_BcCw_E: nuCb, nuC, nuW, mCW, mWC, mWC2, mCW2, Ts, Te, O
		func = dadi_define_models_arnoux18_Paper.IM2_BcCw_E
		
		for optimizationstate in opt_list:
			print optimizationstate

			if optimizationstate == "anneal_hot":
				params = (nuCb_start, nuC_start, nuW_start, mCW_start, mWC_start, mCW2_start, mWC2_start, Ts_start, Te_start, O_start)
			else :
				params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9])
		
			# The upper_bound array is for use in optimization. Occasionally the optimizer
			# will try wacky parameter values. We in particular want to exclude values with
			# very long times, as they will take a long time to evaluate.
			upper_bound = [nuCb_max, nuC_max, nuW_max, mCW_max, mWC_max, mCW2_max, mWC2_max, Ts_max, Te_max, O_max]
			lower_bound = [nuCb_min, nuC_min, nuW_min, mCW_min, mWC_min, mCW2_min, mWC2_min, Ts_min, Te_min, O_min]

			done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, output_file_2, output_file_3, namemodel, ll_opt_dic, nbparam_dic,
								  nompopC=nompopC, nompopW=nompopW, params=params, fixed_params=None, lower_bound=lower_bound, 
								  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
								  outputname= outputname + "/" + outputname,  outputname2=outputname,
								  verbose=verbose, maxiter=maxiter, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, 
								  schedule= "cauchy")
		if done: print(("\n" + namemodel + " : done\n"))


	if namemodel == "IM_C_BcCw_E":

		# IM_C_BcCw_E: nuCb, nuC, nuW, mCW, mWC, Ts, Tb, Te, O
		func = dadi_define_models_arnoux18_Paper.IM_C_BcCw_E

		for optimizationstate in opt_list:
			print optimizationstate

			if optimizationstate == "anneal_hot":
				params = (nuCb_start, nuC_start, nuW_start, mCW_start, mWC_start, Ts_start, Tb_start, Te_start, O_start)
			else :
				params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8])
		
			# The upper_bound array is for use in optimization. Occasionally the optimizer
			# will try wacky parameter values. We in particular want to exclude values with
			# very long times, as they will take a long time to evaluate.
			upper_bound = [nuCb_max, nuC_max, nuW_max, mCW_max, mWC_max, Ts_max, Tb_max, Te_max, O_max]
			lower_bound = [nuCb_min, nuC_min, nuW_min, mCW_min, mWC_min, Ts_min, Tb_min, Te_min, O_min]

			done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, output_file_2, output_file_3, namemodel, ll_opt_dic, nbparam_dic,
								  nompopC=nompopC, nompopW=nompopW, params=params, fixed_params=None, lower_bound=lower_bound, 
								  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
								  outputname= outputname + "/" + outputname,  outputname2=outputname,
								  verbose=verbose, maxiter=maxiter, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, 
								  schedule= "cauchy")
		if done: print(("\n" + namemodel + " : done\n"))


	if namemodel == "IM2_C_BcCw_E":

		# IM2_C_BcCw_E: nuCb, nuC, nuW, mCW, mWC, mWC2, mCW2, Ts, Tb, Te, O
		func = dadi_define_models_arnoux18_Paper.IM2_C_BcCw_E

		for optimizationstate in opt_list:
			print optimizationstate

			if optimizationstate == "anneal_hot":
				params = (nuCb_start, nuC_start, nuW_start, mCW_start, mWC_start, mCW2_start, mWC2_start, Ts_start, Tb_start, Te_start, O_start)
			else :
				params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10])
		
			# The upper_bound array is for use in optimization. Occasionally the optimizer
			# will try wacky parameter values. We in particular want to exclude values with
			# very long times, as they will take a long time to evaluate.
			upper_bound = [nuCb_max, nuC_max, nuW_max, mCW_max, mWC_max, mCW2_max, mWC2_max, Ts_max, Tb_max, Te_max, O_max]
			lower_bound = [nuCb_min, nuC_min, nuW_min, mCW_min, mWC_min, mCW2_min, mWC2_min, Ts_min, Tb_min, Te_min, O_min]

			done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, output_file_2, output_file_3, namemodel, ll_opt_dic, nbparam_dic,
								  nompopC=nompopC, nompopW=nompopW, params=params, fixed_params=None, lower_bound=lower_bound, 
								  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
								  outputname= outputname + "/" + outputname,  outputname2=outputname,
								  verbose=verbose, maxiter=maxiter, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, 
								  schedule= "cauchy")
		if done: print(("\n" + namemodel + " : done\n"))


	if namemodel == "IM_E":

		# IM_E: nuC, nuW, mCW, mWC, Ts, O
		func = dadi_define_models_arnoux18_Paper.IM_E

		for optimizationstate in opt_list:
			print optimizationstate

			if optimizationstate == "anneal_hot":
					params = (nuC_start, nuW_start, mCW_start, mWC_start, Ts_start, O_start)
			else :
					params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])

			# The upper_bound array is for use in optimization. Occasionally the optimizer
			# will try wacky parameter values. We in particular want to exclude values with
			# very long times, as they will take a long time to evaluate.
			upper_bound = [nuC_max, nuW_max, mCW_max, mWC_max, Ts_max, O_max]
			lower_bound = [nuC_min, nuW_min, mCW_min, mWC_min, Ts_min, O_min]

			done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, output_file_2, output_file_3, namemodel, ll_opt_dic, nbparam_dic,
								  nompopC=nompopC, nompopW=nompopW, params=params, fixed_params=None, lower_bound=lower_bound, 
								  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
								  outputname= outputname + "/" + outputname,  outputname2=outputname,
								  verbose=verbose, maxiter=maxiter, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, 
								  schedule= "cauchy")
		if done: print(("\n" + namemodel + " : done\n"))	   


	if namemodel == "IM_E_E":

		# IM_E_E: nuCe, nuC, nuWe, nuW, mCW, mWC, Ts, Te, O
		func = dadi_define_models_arnoux18_Paper.IM_E_E

		for optimizationstate in opt_list:
			print optimizationstate

			if optimizationstate == "anneal_hot":
				params = (nuCe_start, nuC_start, nuWe_start, nuW_start, mCW_start, mWC_start, Ts_start, Te_start, O_start)
			else :
				params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8])
		
			# The upper_bound array is for use in optimization. Occasionally the optimizer
			# will try wacky parameter values. We in particular want to exclude values with
			# very long times, as they will take a long time to evaluate.
			upper_bound = [nuCe_max, nuC_max, nuWe_max, nuW_max, mCW_max, mWC_max, Ts_max, Te_max, O_max]
			lower_bound = [nuCe_min, nuC_min, nuWe_min, nuW_min, mCW_min, mWC_min, Ts_min, Te_min, O_min]

			done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, output_file_2, output_file_3, namemodel, ll_opt_dic, nbparam_dic,
								  nompopC=nompopC, nompopW=nompopW, params=params, fixed_params=None, lower_bound=lower_bound, 
								  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
								  outputname= outputname + "/" + outputname,  outputname2=outputname,
								  verbose=verbose, maxiter=maxiter, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, 
								  schedule= "cauchy")
		if done: print(("\n" + namemodel + " : done\n"))


	if namemodel == "IM_C_E_E":

		# IM_C_E_E: nuCe, nuC, nuWe, nuW, mCW, mWC, Ts, Tb, Te, O
		func = dadi_define_models_arnoux18_Paper.IM_C_E_E

		for optimizationstate in opt_list:
			print optimizationstate

			if optimizationstate == "anneal_hot":
				params = (nuCe_start, nuC_start, nuWe_start, nuW_start, mCW_start, mWC_start, Ts_start, Te_start, Tb_start, O_start)
			else :
				params = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9])
		
			# The upper_bound array is for use in optimization. Occasionally the optimizer
			# will try wacky parameter values. We in particular want to exclude values with
			# very long times, as they will take a long time to evaluate.
			upper_bound = [nuCe_max, nuC_max, nuWe_max, nuW_max, mCW_max, mWC_max, Ts_max, Te_max, Tb_max, O_max]
			lower_bound = [nuCe_min, nuC_min, nuWe_min, nuW_min, mCW_min, mWC_min, Ts_min, Te_min, Tb_min, O_min]

			done, ll_opt_dic, nbparam_dic, popt = callmodel(func, data, output_file, output_file_2, output_file_3, namemodel, ll_opt_dic, nbparam_dic,
								  nompopC=nompopC, nompopW=nompopW, params=params, fixed_params=None, lower_bound=lower_bound, 
								  upper_bound=upper_bound,  pts_l=pts_l, ns=ns,
								  outputname= outputname + "/" + outputname,  outputname2=outputname,
								  verbose=verbose, maxiter=maxiter, Tini=Tini, Tfin=Tfin, learn_rate=learn_rate, 
								  schedule= "cauchy")
		if done: print(("\n" + namemodel + " : done\n"))


output_file.close()
output_file_2.close()
output_file_3.close()

