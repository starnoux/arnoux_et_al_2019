# -*- coding: utf-8 -*-

# Import des différentes librairie

import os
import pylab
import numpy
from numpy import array
import dadi

os.chdir('/Users/christophersauvage/Desktop/ARCAD/analysis/_SNP/dadi/spectres/_2pops_6-9')

import new_model

#Import des données
os.chdir('/Users/christophersauvage/Desktop/ARCAD/analysis/_SNP/dadi/spectres/_2pops_6-9')

fs_data = dadi.Spectrum.from_file('Tomato_unfolded.fs')

#Variables générales
ns = fs_data.sample_sizes
npop = fs_data.Npop
xmin = max(ns)+10 #taille de la plus petite grille
pts = [xmin,xmin+10,xmin+20]

# Plot 1D des données
# projection 1D du spectre sauvage
fs_1d_wild=fs_data.marginalize([0])
dadi.Plotting.plot_1d_fs(fs_1d_wild)

# projection 1D du spectre cultivé
pylab.figure()
fs_1d_cult=fs_data.marginalize([1])
dadi.Plotting.plot_1d_fs(fs_1d_cult)

# Plot 2D des données
pylab.figure()
dadi.Plotting.plot_single_2d_sfs(fs_data,vmin=0.1)


# Statsitiques de bases: utile pour calibrer les valeurs initiales
pi_C = fs_data.marginalize([1]).pi()
pi_W = fs_data.marginalize([0]).pi()
thetaW_C = fs_data.marginalize([1]).Watterson_theta()
thetaW_W = fs_data.marginalize([0]).Watterson_theta()
Dtaj_C = fs_data.marginalize([1]).Tajima_D()
Dtaj_W = fs_data.marginalize([0]).Tajima_D()
Fst = fs_data.Fst()
SC = fs_data.marginalize([1]).S
SW = fs_data.marginalize([0]).S

### Valeurs initiales possibles. On peut aussi les fixer à la main

nuCi = pi_C/pi_W # ici l'idée est que le sauvage est proche de l'ancetre
nuWi = 1-Dtaj_W # ici l'idée est qu'un D < 0 correspond à une expension et que c'est donc un moyen (très grossier) de déterminer la taille du sauvage par rapport à l'ancetre
mi = (1/Fst - 1)/2 # pour le modèle avec migration symétrique. Attendu dans le modèle en iles (encore une fois très grossier)
Ti = 0.5 # le temps de coalescence moyen de deux séquences dans le modèle standard (ici T est en 4Ne générations)
erri = 0.05 # erreur
Pi = 0.1 # fraction de gènes sous sélection. Rq: ne pas utiliser pi = 3.14..!


### Fonction pour calculer la vraisemblance du modèle saturé

def lnLmax(fs):
    ns = fs_data.sample_sizes
    lnL = 0
    imax = ns[0]-1
    jmax = ns[1]-1
    for i in range(0,imax):
        for j in range(0,jmax):
            k = fs[i,j]
            if k > 0:
                lnL += -k + k*log(k) -numpy.math.lgamma(k+1)
    return lnL
 
lnLmax(fs_data)


###################################
###      MODELS IN A ROW        ###
###################################

#### SPLIT BOTTLENECK MIGRATION ASYMETRIC ERROR ####

##################################
### ONLY ONE CATEGORY OF GENES ###
##################################

###############
### ROUND 1 ###
###############

func = new_model.split_bottleneck_migration

# nuC1, nuC2, nuW, mCW, mWC, T1, T2, e
params = (0.02141713, 4.26487927, 1.76771081, 2.50919151, 0.04575117, 9.1988774 , 0.08353503, 0.15498357)

# define the upper and lower bound of the parameters
up = [10,10,10,20,20,10,10,0.5]
low = [0.00001,0.00001,0.00001,0.00001,0.00001,0.001,0.001,0]

# creation de la version extrapolee du modele demographique
func_ex = dadi.Numerics.make_extrap_log_func(func)
model = func_ex(params, ns, pts)

# vraisemblance initiale
ll_init = dadi.Inference.ll_multinom(model, fs_data)
print 'Model_init log-likelihood:', ll_init

# perturbation initiale
p0 = dadi.Misc.perturb_params(params, fold=1, lower_bound=low, upper_bound=up)

# optimisation

# using the anneal method (large range of parameters variation)
p_migration_err = dadi.Inference.optimize_anneal(p0, fs_data, func_ex, pts,lower_bound=low,upper_bound=up, maxiter=1, verbose=1, Tini=50, Tfin=0, learn_rate=0.01, schedule="cauchy")

print 'Optimzed parameters', repr(p_migration_err)
model_migration_err = func_ex(p_migration_err, ns, pts)

# paramètres au ML
p_migration_err
# vraisemblance max
ll_migration_err = dadi.Inference.ll_multinom(model_migration_err, fs_data)
print 'Model_optimized log-likelihood:', ll_migration_err

# valeur optimale de theta
theta_migration_err = dadi.Inference.optimal_sfs_scaling(model_migration_err, fs_data)
print 'Theta:', theta_migration_err

# plot
pylab.figure()
dadi.Plotting.plot_2d_comp_multinom(model_migration_err, fs_data, vmin=10, resid_range=5)

###############
### ROUND 2 ###
###############

func = new_model.split_bottleneck_migration

# nuC1, nuC2, nuW, mCW, mWC, T1, T2, e
params = array([p_migration_err[0],p_migration_err[1],p_migration_err[2],p_migration_err[3],p_migration_err[4],p_migration_err[5],p_migration_err[6],p_migration_err[7]])
# params = (4.57771448e-03, 4.54648931e-01, 5.20974260e-01, 5.07082797e+00, 2.78374778e-01, 4.61302922e-01, 4.56521474e-02, 1.35222389e-01)

# define the upper and lower bound of the parameters
up = [10,10,10,20,20,10,10,0.5]
low = [0.00001,0.00001,0.00001,0.00001,0.00001,0.001,0.001,0]

# creation de la version extrapolee du modele demographique
func_ex = dadi.Numerics.make_extrap_log_func(func)
model = func_ex(params, ns, pts)

# vraisemblance initiale
ll_init = dadi.Inference.ll_multinom(model, fs_data)
print 'Model_init log-likelihood:', ll_init

# perturbation initiale
p0 = dadi.Misc.perturb_params(params, fold=1, lower_bound=low, upper_bound=up)

# optimisation

# using the anneal method (large range of parameters variation)
p_migration_err_2 = dadi.Inference.optimize_anneal(p0, fs_data, func_ex, pts,lower_bound=low,upper_bound=up, maxiter=5, verbose=1, Tini=50, Tfin=0, learn_rate=0.01, schedule="cauchy")

print 'Optimzed parameters', repr(p_migration_err_2)
model_migration_err_2 = func_ex(p_migration_err_2, ns, pts)

# paramètres au ML
p_migration_err_2
# vraisemblance max
ll_migration_err_2 = dadi.Inference.ll_multinom(model_migration_err_2, fs_data)
print 'Model_optimized log-likelihood:', ll_migration_err_2

# valeur optimale de theta
theta_migration_err_2 = dadi.Inference.optimal_sfs_scaling(model_migration_err_2, fs_data)
print 'Theta:', theta_migration_err_2

# plot
pylab.figure()
dadi.Plotting.plot_2d_comp_multinom(model_migration_err_2, fs_data, vmin=10, resid_range=5)

###############
### ROUND 3 ###
###############

func = new_model.split_bottleneck_migration

# nuC1, nuC2, nuW, mCW, mWC, T1, T2, e
params = array([p_migration_err_2[0],p_migration_err_2[1],p_migration_err_2[2],p_migration_err_2[3],p_migration_err_2[4],p_migration_err_2[5],p_migration_err_2[6],p_migration_err_2[7]])

# define the upper and lower bound of the parameters
up = [10,10,10,20,20,10,10,0.5]
low = [0.00001,0.00001,0.00001,0.00001,0.00001,0.001,0.001,0]

# creation de la version extrapolee du modele demographique
func_ex = dadi.Numerics.make_extrap_log_func(func)
model = func_ex(params, ns, pts)

# vraisemblance initiale
ll_init = dadi.Inference.ll_multinom(model, fs_data)
print 'Model_init log-likelihood:', ll_init

# perturbation initiale
p0 = dadi.Misc.perturb_params(params, fold=1, lower_bound=low, upper_bound=up)

# optimisation

# using the anneal method (large range of parameters variation)
p_migration_err_3 = dadi.Inference.optimize_log_fmin(p0, fs_data, func_ex, pts,lower_bound=low,upper_bound=up, maxiter=5, verbose=1)

print 'Optimzed parameters', repr(p_migration_err_3)
model_migration_err_3 = func_ex(p_migration_err_3, ns, pts)

# paramètres au ML
p_migration_err_3
# vraisemblance max
ll_migration_err_3 = dadi.Inference.ll_multinom(model_migration_err_3, fs_data)
print 'Model_optimized log-likelihood:', ll_migration_err_3

# valeur optimale de theta
theta_migration_err_3 = dadi.Inference.optimal_sfs_scaling(model_migration_err_3, fs_data)
print 'Theta:', theta_migration_err_3

# plot
pylab.figure()
dadi.Plotting.plot_2d_comp_multinom(model_migration_err_3, fs_data, vmin=10, resid_range=5)




###############################
### TWO CATEGORIES OF GENES ###
###############################

###############
### ROUND 1 ###
###############

func = new_model.split_bottleneck_migration_2cat_error

# nuC1x,nuC1y,nuC2,nuW,mCW,mWC,T1,T2,p,e
params = (1.82573289e-05, 5.21611269e-03, 4.32263364e-01, 5.61241520e-01, 5.68756922e+00, 4.05973583e-01, 2.93795323e-01, 4.54987283e-02, 6.47386758e-02, 1.39968519e-01)

up = [10,10,10,10,20,20,10,10,0.5,0.5]
low = [0,0,0,0.0001,0.0001,0.0001,0.001,0.001,0,0]

# creation de la version extrapolee du modele demographique
func_ex = dadi.Numerics.make_extrap_log_func(func)
model = func_ex(params, ns, pts)

# vraisemblance initiale
ll_init = dadi.Inference.ll_multinom(model, fs_data)
print 'Model_init log-likelihood:', ll_init

p0 = dadi.Misc.perturb_params(params, fold=1, lower_bound=low, upper_bound=up)

# optimisation

# using the anneal method (large range of parameters variation)
p_migration_err_2cat_error = dadi.Inference.optimize_anneal(p0, fs_data, func_ex, pts,lower_bound=low,upper_bound=up, maxiter=3, verbose=1, Tini=50, Tfin=0, learn_rate=0.01, schedule="cauchy")

print 'Optimzed parameters', repr(p_migration_err)
model_migration_err_2cat_error = func_ex(p_migration_err_2cat_error, ns, pts)

# paramètres au ML
p_migration_err_2cat_error
# vraisemblance max
ll_migration_err_2cat_error = dadi.Inference.ll_multinom(model_migration_err_2cat_error, fs_data)
print 'Model_optimized log-likelihood:', ll_migration_err

# valeur optimale de theta
theta_migration_err_2cat_error = dadi.Inference.optimal_sfs_scaling(model_migration_err_2cat_error, fs_data)
print 'Theta:', theta_migration_err_2cat_error

# plot
pylab.figure()
dadi.Plotting.plot_2d_comp_multinom(model_migration_err_2cat_error, fs_data, vmin=10, resid_range=5)

###############
### ROUND 2 ###
###############

func = new_model.split_bottleneck_migration_2cat_error

# nuC1x,nuC1y,nuC2,nuW,mCW,mWC,T1,T2,p,e
params = array([p_migration_err_2cat_error[0],p_migration_err_2cat_error[1],p_migration_err_2cat_error[2],p_migration_err_2cat_error[3],p_migration_err_2cat_error[4],p_migration_err_2cat_error[5],p_migration_err_2cat_error[6],p_migration_err_2cat_error[7],p_migration_err_2cat_error[8],p_migration_err_2cat_error[9]])

up = [10,10,10,10,20,20,10,10,0.5,0.5]
low = [0,0,0,0.0001,0.0001,0.0001,0.001,0.001,0,0]

# creation de la version extrapolee du modele demographique
func_ex = dadi.Numerics.make_extrap_log_func(func)
model = func_ex(params, ns, pts)

# vraisemblance initiale
ll_init = dadi.Inference.ll_multinom(model, fs_data)
print 'Model_init log-likelihood:', ll_init

p0 = dadi.Misc.perturb_params(params, fold=1, lower_bound=low, upper_bound=up)

# optimisation

# using the anneal method (large range of parameters variation)
p_migration_err_2cat_error_2 = dadi.Inference.optimize_anneal(p0, fs_data, func_ex, pts,lower_bound=low,upper_bound=up, maxiter=3, verbose=1, Tini=50, Tfin=0, learn_rate=0.01, schedule="cauchy")

print 'Optimzed parameters', repr(p_migration_err)
model_migration_err_2cat_error_2 = func_ex(p_migration_err_2cat_error_2, ns, pts)

# paramètres au ML
p_migration_err_2cat_error_2
# vraisemblance max
ll_migration_err_2cat_error_2 = dadi.Inference.ll_multinom(model_migration_err_2cat_error_2, fs_data)
print 'Model_optimized log-likelihood:', ll_migration_err_2

# valeur optimale de theta
theta_migration_err_2cat_error_2 = dadi.Inference.optimal_sfs_scaling(model_migration_err_2cat_error_2, fs_data)
print 'Theta:', theta_migration_err_2cat_error_2

# plot
pylab.figure()
dadi.Plotting.plot_2d_comp_multinom(model_migration_err_2cat_error_2, fs_data, vmin=10, resid_range=5)

###############
### ROUND 3 ###
###############

func = new_model.split_bottleneck_migration_2cat_error

# nuC1x,nuC1y,nuC2,nuW,mCW,mWC,T1,T2,p,e
params = array([p_migration_err_2cat_error_2[0],p_migration_err_2cat_error_2[1],p_migration_err_2cat_error_2[2],p_migration_err_2cat_error_2[3],p_migration_err_2cat_error_2[4],p_migration_err_2cat_error_2[5],p_migration_err_2cat_error_2[6],p_migration_err_2cat_error_2[7],p_migration_err_2cat_error_2[8],p_migration_err_2cat_error_2[9]])

up = [10,10,10,10,20,20,10,10,0.5,0.5]
low = [0,0,0,0.0001,0.0001,0.0001,0.001,0.001,0,0]

# creation de la version extrapolee du modele demographique
func_ex = dadi.Numerics.make_extrap_log_func(func)
model = func_ex(params, ns, pts)

# vraisemblance initiale
ll_init = dadi.Inference.ll_multinom(model, fs_data)
print 'Model_init log-likelihood:', ll_init

p0 = dadi.Misc.perturb_params(params, fold=1, lower_bound=low, upper_bound=up)

# optimisation

# using the anneal method (large range of parameters variation)
p_migration_err_2cat_error_3 = dadi.Inference.optimize_log_fmin(p0, fs_data, func_ex, pts,lower_bound=low,upper_bound=up, maxiter=3, verbose=1)

print 'Optimzed parameters', repr(p_migration_err)
model_migration_err_2cat_error_3 = func_ex(p_migration_err_2cat_error_3, ns, pts)

# paramètres au ML
p_migration_err_2cat_error_3
# vraisemblance max
ll_migration_err_2cat_error_3 = dadi.Inference.ll_multinom(model_migration_err_2cat_error_3, fs_data)
print 'Model_optimized log-likelihood:', ll_migration_err_3

# valeur optimale de theta
theta_migration_err_2cat_error_3 = dadi.Inference.optimal_sfs_scaling(model_migration_err_2cat_error_3, fs_data)
print 'Theta:', theta_migration_err_2cat_error_3

# plot
pylab.figure()
dadi.Plotting.plot_2d_comp_multinom(model_migration_err_2cat_error_3, fs_data, vmin=10, resid_range=5)
