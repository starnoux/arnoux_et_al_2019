# -*- coding: utf-8 -*-
# 19.07.2017

# Import des differentes librairies
import matplotlib
matplotlib.use('Agg')
import os
import matplotlib.pyplot as plt
import pylab
import numpy
from numpy import array
import scipy
import dadi
from math import log

os.chdir('/path/to/dadi/')

import models2D
import models2D_error
import models2D_2cat
import models2D_2cat_error
import new_model

# Import des donnees
os.chdir('/path/to/sfs/4p-master/')

fs_data = dadi.Spectrum.from_file('AFS-U_WILD_CROP.dadi.txt')

#Variables generales

ns = fs_data.sample_sizes
npop = fs_data.Npop
xmin = max(ns)+10 #taille de la plus petite grille
pts = [xmin,xmin+10,xmin+20]

# Plot 1D des donnees
# projection 1D du spectre sauvage
fs_1d_wild=fs_data.marginalize([0])
fig = plt.figure()
dadi.Plotting.plot_1d_fs(fs_1d_wild)
fig.savefig('/data_nas1/wgs/solanaceae/Solution/starnoux/Tools/4p-master/MM_HETtest_2018_1D_Wild.png')
plt.close(fig)

# projection 1D du spectre cultive
#pylab.figure()

fs_1d_cult=fs_data.marginalize([1])

fig = plt.figure()
plt.figure(dadi.Plotting.plot_1d_fs(fs_1d_cult))
fig.savefig('/data_nas1/wgs/solanaceae/Solution/starnoux/Tools/4p-master/MM_HETtest_2018_1D_CROP.png')
plt.close(fig)


# Plot 2D des donnees
fig = plt.figure()
dadi.Plotting.plot_single_2d_sfs(fs_data,vmin=0.1)
fig.savefig('/data_nas1/wgs/solanaceae/Solution/starnoux/Tools/4p-master/MM_HETtest_2018_2D.png')
plt.close(fig)

...