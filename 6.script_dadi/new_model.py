# -*- coding: utf-8 -*-
import numpy
from numpy import array
import dadi

### NEW MODELS ###

def split_bottleneck_migration(params,(n1,n2),pts):
    nuC1, nuC2, nuW, mCW, mWC, T1, T2, e = params
    """
    nuC1 : Size of the cultivated population during bottleneck
    nuC2 : Size of the cultivated population after bottleneck
    nuW  : Size of the wild population
    mCW	 : Migration rate from wild to cultivated
    mWC  : Migration rate from cultivated to wild
    T1	 : Duration of the bottleneck
    T2	 : Time since bottleneck 
    e	 : Error rate
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # split
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # phase 1
    phi = dadi.Integration.two_pops(phi, xx, T1, nuC1, nuW, m12=mCW, m21=mWC)
    # phase 2
    phi = dadi.Integration.two_pops(phi, xx, T2, nuC2, nuW, m12=mCW, m21=mWC)
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    # adding error to the spectrum
    fs_rev = dadi.Numerics.reverse_array(fs)
    fs_tot = (1-e)*fs + e*fs_rev
    return fs_tot
    
    
def split_bottleneck_migration_2cat_error((nuC1x,nuC1y,nuC2,nuW,mCW,mWC,T1,T2,p,e),(n1,n2),pts):
    """
    nuC1x : Size of the cultivated population during bottleneck, gene category X
    nuC1y : Size of the cultivated population during bottleneck, gene category Y
    nuC2  : Size of the cultivated population after bottleneck
    nuW   : Size of the wild population
    mCW	  : Migration rate from wild to cultivated
    mWC   : Migration rate from cultivated to wild
    T1	  : Duration of the bottleneck
    T2	  : Time since bottleneck 
    p	  : Proportion of genes under category Y
    e    : Error rate
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # CAT X
    # split
    phiX = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # phase 1
    phiX = dadi.Integration.two_pops(phiX, xx, T1, nu1=nuC1x, nu2=nuW, m12=mCW, m21=mWC)
    # phase 2
    phiX = dadi.Integration.two_pops(phiX, xx, T2, nu1=nuC2, nu2=nuW, m12=mCW, m21=mWC)
    # computing fsX
    fsX = dadi.Spectrum.from_phi(phiX, (n1,n2), (xx,xx))
    # CAT Y
    # split
    phiY = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # phase 1
    phiY = dadi.Integration.two_pops(phiY, xx, T1, nu1=nuC1y, nu2=nuW, m12=mCW, m21=mWC)
    # phase 2
    phiY = dadi.Integration.two_pops(phiY, xx, T2, nu1=nuC2, nu2=nuW, m12=mCW, m21=mWC)
    # computing fsY
    fsY = dadi.Spectrum.from_phi(phiY, (n1,n2), (xx,xx))
    # complete spectrum
    fs = (1-p)*fsX + p*fsY
    # adding error to the spectrum
    fs_rev = dadi.Numerics.reverse_array(fs)
    fs_tot = (1-e)*fs + e*fs_rev
    return fs_tot
    
    
    #### General model with two time periods ####

def general_2T((nuC1,nuW1,nuC2,nuW2,mCW1,mWC1,mCW2,mWC2,T1,T2,e),ns,pts):
    """
    nuC1 : Size of the cultivated population during phase 1
    nuC2 : Size of the cultivated population during phase 2
    nuW1 : Size of the wild population during phase 1
    nuW2 : Size of the wild population during phase 2
    mCW1 : Migration rate from wild to cultivated during phase 1
    mWC1 : Migration rate from cultivated to wild during phase 1
    mCW2 : Migration rate from wild to cultivated during phase 2
    mWC2 : Migration rate from cultivated to wild during phase 2
    T1	 : Duration of the bottleneck
    T2	 : Time since bottleneck 
    e    : error rate
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Split
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # phase 1
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1=nuC1, nu2=nuW1, m12=mCW1, m21=mWC1)
    # phase 2
    phi = dadi.Integration.two_pops(phi, xx, T2, nu1=nuC2, nu2=nuW2, m12=mCW2, m21=mWC2)
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    # adding error to the spectrum
    fs_rev = dadi.Numerics.reverse_array(fs)
    fs_tot = (1-e)*fs + e*fs_rev
    return fs_tot
    