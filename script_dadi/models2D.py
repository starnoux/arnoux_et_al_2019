# -*- coding: utf-8 -*-
import numpy
from numpy import array
import dadi

# Preliminary remark: The terms used are cultivated (pop0) and wild (pop1) populations but they can apply to any kind of populations
# nu = population size relative to the ancestral reference population (Nref)
# T: time in 4Nref population


#### MAIN MODELS FOR ARCAD PROJECT ####


#### Split of the ancestral population into two new populations of different sizes. No migration. ####

def split((nuC,nuW,T),(n1,n2),pts):
    """
    nuC : Size of the cultivated population
    nuW : Size of the wild population 
    T	: Total time since the split
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # split
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # divergence
    phi = dadi.Integration.two_pops(phi, xx, T, nu1=nuC, nu2=nuW, m12=0, m21=0)
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs



#### Split of the ancestral population into two new populations of different sizes. Symetrical migration. ####

def split_migration_sym((nuC,nuW,m,T),(n1,n2),pts):
    """
    nuC : Size of the cultivated population
    nuW : Size of the wild population 
    m   : Migration rate
    T	: Total time since the split
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # split
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # divergence
    phi = dadi.Integration.two_pops(phi, xx, T, nu1=nuC, nu2=nuW, m12=m, m21=m)
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs
    
    
    
#### Split of the ancestral population into two new populations of different sizes. Asymetrical migration. ####

def split_migration((nuC,nuW,mCW,mWC,T),(n1,n2),pts):
    """
    nuC : Size of the cultivated population
    nuW : Size of the wild population 
    mCW	 : Migration rate from wild to cultivated
    mWC  : Migration rate from cultivated to wild
    T	: Total time since the split
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # split
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # divergence
    phi = dadi.Integration.two_pops(phi, xx, T, nu1=nuC, nu2=nuW, m12=mCW, m21=mWC)
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs
    


#### Split of the ancestral population into two new populations of different sizes. Two sizes in the cultivated population. No migration ####

def split_bottleneck((nuC1,nuC2,nuW,T1,T2),(n1,n2),pts):
    """
    nuC1 : Size of the cultivated population during bottleneck
	nuC2 : Size of the cultivated population after bottleneck
    nuW  : Size of the wild population 
    T1	 : Duration of the bottleneck
	T2	 : Time since bottleneck 
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # split
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # phase 1
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1=nuC1, nu2=nuW, m12=0, m21=0)
    # phase 2
    phi = dadi.Integration.two_pops(phi, xx, T2, nu1=nuC2, nu2=nuW, m12=0, m21=0)
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs


#### Split of the ancestral population into two new populations of different sizes. Two sizes in the cultivated population. Symetrical migration ####

def split_bottleneck_migration_sym((nuC1,nuC2,nuW,m,T1,T2),(n1,n2),pts):
    """
    nuC1 : Size of the cultivated population during bottleneck
    nuC2 : Size of the cultivated population after bottleneck
    nuW  : Size of the wild population
    m	 : Migration rate
    T1	 : Duration of the bottleneck
    T2	 : Time since bottleneck 
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # split
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # phase 1
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1=nuC1, nu2=nuW, m12=m, m21=m)
    # phase 2
    phi = dadi.Integration.two_pops(phi, xx, T2, nu1=nuC2, nu2=nuW, m12=m, m21=m)
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs



#### Split of the ancestral population into two new populations of different sizes. Two sizes in the cultivated population. Asymetrical migration ####

def split_bottleneck_migration((nuC1,nuC2,nuW,mCW,mWC,T1,T2),(n1,n2),pts):
    """
    nuC1 : Size of the cultivated population during bottleneck
    nuC2 : Size of the cultivated population after bottleneck
    nuW  : Size of the wild population
    mCW	 : Migration rate from wild to cultivated
    mWC  : Migration rate from cultivated to wild
    T1	 : Duration of the bottleneck
    T2	 : Time since bottleneck 
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # split
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # phase 1
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1=nuC1, nu2=nuW, m12=mCW, m21=mWC)
    # phase 2
    phi = dadi.Integration.two_pops(phi, xx, T2, nu1=nuC2, nu2=nuW, m12=mCW, m21=mWC)
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs





#### ADDITIONAL MODELS ####



#### General model with two time periods ####

def general_2T((nuC1,nuW1,nuC2,nuW2,mCW1,mWC1,mCW2,mWC2,T1,T2),ns,pts):
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
    return fs


#### Split and population growth in the cultivated population. No migration ####
	
def growth((nuCi,nuCf,nuW,T),(n1,n2),pts):
    """
    nuCi : Initial size of the cultivated population
    nuCf : Final size of the cultivated population
    nuW  : Size of the wild population
    T	 : Total time since the split
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # split
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # bottleneck and growth
    nuC_func = lambda t: nuCi*(nuCf/nuCi)**(t/T)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1=nuC_func, nu2=nuW, m12=0, m21=0)
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs


#### Split and population growth in the cultivated population. Symetrical migration ####
	
def growth_migration_sym((nuCi,nuCf,nuW,m,T),(n1,n2),pts):
    """
    nuCi : Initial size of the cultivated population
	nuCf : Final size of the cultivated population
    nuW  : Size of the wild population
    m	 : Migration rate
    T	 : Total time since the split
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # split
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # bottleneck and growth
    nuC_func = lambda t: nuCi*(nuCf/nuCi)**(t/T)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1=nuC_func, nu2=nuW, m12=m, m21=m)
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs


#### Split and population growth in the cultivated population. Asymetrical migration ####
	
def growth_migration((nuCi,nuCf,nuW,mCW,mWC,T),(n1,n2),pts):
    """
    nuCi : Initial size of the cultivated population
    nuCf : Final size of the cultivated population
    nuW  : Size of the wild population
    mCW	 : Migration rate from wild to cultivated
    mWC  : Migration rate from cultivated to wild
    T	 : Total time since the split
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # split
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # bottleneck and growth
    nuC_func = lambda t: nuCi*(nuCf/nuCi)**(t/T)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1=nuC_func, nu2=nuW, m12=mCW, m21=mWC)
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs



#### Split. No migration during phase 1. Migration during phase 2 ####


def admixture((nuC,nuW,mCW,mWC,T1,T2),(n1,n2),pts):
    """
    nuC : Size of the cultivated population
    nuW : Size of the wild population
    mCW	 : Migration rate from wild to cultivated
    mWC  : Migration rate from cultivated to wild
    T1	 : Duration of phase 1
    T2	 : Duration of phase 2
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # split
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # divergence without migration
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1=nuC, nu2=nuW, m12=0, m21=0)
    # Secondary contact
    phi = dadi.Integration.two_pops(phi, xx, T2, nu1=nuC, nu2=nuW, m12=mCW, m21=mWC)
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs