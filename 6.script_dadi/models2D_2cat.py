# -*- coding: utf-8 -*-
import numpy
from numpy import array
import dadi

# Same models as in models2D with two categories of genes
# The difference between the two categories is the intensity of the bottleneck in the cultivated population
# All other parameters are equal between the two categories


# Preliminary remark: The terms used are cultivated (pop0) and wild (pop1) populations but they can apply to any kind of populations
# nu = population size relative to the ancestral reference population (Nref)
# T: time in 4Nref population



#### MAIN MODELS FOR ARCAD PROJECT ####


#### Split of the ancestral population into two new populations of different sizes. No migration. ####

def split((nuCx,nuCy,nuW,T,p),(n1,n2),pts):
    """
    nuCx : Size of the cultivated population, gene category X
    nuCy : Size of the cultivated population, gene category Y
    nuW  : Size of the wild population 
    T	 : Total time since the split
    p	 : Proportion of genes under category Y
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # CAT X
    # split
    phiX = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # divergence
    phiX = dadi.Integration.two_pops(phiX, xx, T, nu1=nuCx, nu2=nuW, m12=0, m21=0)
    # computing fsX
    fsX = dadi.Spectrum.from_phi(phiX, (n1,n2), (xx,xx))
    # CAT Y
    # split
    phiY = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # divergence
    phiY = dadi.Integration.two_pops(phiY, xx, T, nu1=nuCy, nu2=nuW, m12=0, m21=0)
    # computing fsY
    fsY = dadi.Spectrum.from_phi(phiY, (n1,n2), (xx,xx))
    # complete spectrum
    fs = (1-p)*fsX + p*fsY
    return fs



#### Split of the ancestral population into two new populations of different sizes. Symetrical migration. ####

def split_migration_sym((nuCx,nuCy,nuW,m,T,p),(n1,n2),pts):
    """
    nuCx : Size of the cultivated population, gene category X
    nuCy : Size of the cultivated population, gene category Y
    nuW  : Size of the wild population 
    m    : Migration rate
    T	 : Total time since the split
    p	 : Proportion of genes under category Y
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # CAT X
    # split
    phiX = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # divergence
    phiX = dadi.Integration.two_pops(phiX, xx, T, nu1=nuCx, nu2=nuW, m12=m, m21=m)
    # computing fsX
    fsX = dadi.Spectrum.from_phi(phiX, (n1,n2), (xx,xx))
    # CAT Y
    # split
    phiY = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # divergence
    phiY = dadi.Integration.two_pops(phiY, xx, T, nu1=nuCy, nu2=nuW, m12=m, m21=m)
    # computing fsY
    fsY = dadi.Spectrum.from_phi(phiY, (n1,n2), (xx,xx))
    # complete spectrum
    fs = (1-p)*fsX + p*fsY
    return fs



#### Split of the ancestral population into two new populations of different sizes. Asymetrical migration. ####

def split_migration((nuCx,nuCy,nuW,mCW,mWC,T,p),(n1,n2),pts):
    """
    nuCx : Size of the cultivated population, gene category X
    nuCy : Size of the cultivated population, gene category Y
    nuW  : Size of the wild population 
    mCW	  : Migration rate from wild to cultivated
    mWC   : Migration rate from cultivated to wild
    T	 : Total time since the split
    p	 : Proportion of genes under category Y
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # CAT X
    # split
    phiX = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # divergence
    phiX = dadi.Integration.two_pops(phiX, xx, T, nu1=nuCx, nu2=nuW, m12=mCW, m21=mWC)
    # computing fsX
    fsX = dadi.Spectrum.from_phi(phiX, (n1,n2), (xx,xx))
    # CAT Y
    # split
    phiY = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # divergence
    phiY = dadi.Integration.two_pops(phiY, xx, T, nu1=nuCy, nu2=nuW, m12=mCW, m21=mWC)
    # computing fsY
    fsY = dadi.Spectrum.from_phi(phiY, (n1,n2), (xx,xx))
    # complete spectrum
    fs = (1-p)*fsX + p*fsY
    return fs




def split_migration_sym2((nuCx,nuCy,nuW,mx,my,T,p),(n1,n2),pts):
    """
    nuCx : Size of the cultivated population, gene category X
    nuCy : Size of the cultivated population, gene category Y
    nuW  : Size of the wild population 
    m    : Migration rate
    T	 : Total time since the split
    p	 : Proportion of genes under category Y
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # CAT X
    # split
    phiX = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # divergence
    phiX = dadi.Integration.two_pops(phiX, xx, T, nu1=nuCx, nu2=nuW, m12=mx, m21=mx)
    # computing fsX
    fsX = dadi.Spectrum.from_phi(phiX, (n1,n2), (xx,xx))
    # CAT Y
    # split
    phiY = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # divergence
    phiY = dadi.Integration.two_pops(phiY, xx, T, nu1=nuCy, nu2=nuW, m12=my, m21=my)
    # computing fsY
    fsY = dadi.Spectrum.from_phi(phiY, (n1,n2), (xx,xx))
    # complete spectrum
    fs = (1-p)*fsX + p*fsY
    return fs





#### Split of the ancestral population into two new populations of different sizes. Two sizes in the cultivated population. No migration ####

def split_bottleneck((nuC1x,nuC1y,nuC2,nuW,T1,T2,p),(n1,n2),pts):
    """
    nuC1x : Size of the cultivated population during bottleneck, gene category X
    nuC1y : Size of the cultivated population during bottleneck, gene category Y
    nuC2  : Size of the cultivated population after bottleneck
    nuW   : Size of the wild population 
    T1	  : Duration of the bottleneck
    T2	  : Time since bottleneck
    p	  : Proportion of genes under category Y
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # CAT X
    # split
    phiX = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # phase 1
    phiX = dadi.Integration.two_pops(phiX, xx, T1, nu1=nuC1x, nu2=nuW, m12=0, m21=0)
    # phase 2
    phiX = dadi.Integration.two_pops(phiX, xx, T2, nu1=nuC2, nu2=nuW, m12=0, m21=0)
    # computing fsX
    fsX = dadi.Spectrum.from_phi(phiX, (n1,n2), (xx,xx))
    # CAT Y
    # split
    phiY = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # phase 1
    phiY = dadi.Integration.two_pops(phiY, xx, T1, nu1=nuC1x, nu2=nuW, m12=0, m21=0)
    # phase 2
    phiY = dadi.Integration.two_pops(phiY, xx, T2, nu1=nuC2, nu2=nuW, m12=0, m21=0)
    # computing fsX
    fsY = dadi.Spectrum.from_phi(phiY, (n1,n2), (xx,xx))
    # complete spectrum
    fs = (1-p)*fsX + p*fsY
    return fs


#### Split of the ancestral population into two new populations of different sizes. Two sizes in the cultivated population. Symetrical migration ####

def split_bottleneck_migration_sym((nuC1x,nuC1y,nuC2,nuW,m,T1,T2,p),(n1,n2),pts):
    """
    nuC1x : Size of the cultivated population during bottleneck, gene category X
    nuC1y : Size of the cultivated population during bottleneck, gene category Y
    nuW   : Size of the wild population
    m	  : Migration rate
    T1	  : Duration of the bottleneck
    T2	  : Time since bottleneck 
    p	  : Proportion of genes under category Y
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # CAT X
    # split
    phiX = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # phase 1
    phiX = dadi.Integration.two_pops(phiX, xx, T1, nu1=nuC1x, nu2=nuW, m12=m, m21=m)
    # phase 2
    phiX = dadi.Integration.two_pops(phiX, xx, T2, nu1=nuC2, nu2=nuW, m12=m, m21=m)
    # computing fsX
    fsX = dadi.Spectrum.from_phi(phiX, (n1,n2), (xx,xx))
    # CAT Y
    # split
    phiY = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # phase 1
    phiY = dadi.Integration.two_pops(phiY, xx, T1, nu1=nuC1y, nu2=nuW, m12=m, m21=m)
    # phase 2
    phiY = dadi.Integration.two_pops(phiY, xx, T2, nu1=nuC2, nu2=nuW, m12=m, m21=m)
    # computing fsY
    fsY = dadi.Spectrum.from_phi(phiY, (n1,n2), (xx,xx))
    # complete spectrum
    fs = (1-p)*fsX + p*fsY
    return fs



#### Split of the ancestral population into two new populations of different sizes. Two sizes in the cultivated population. Asymetrical migration ####

def split_bottleneck_migration((nuC1x,nuC1y,nuC2,nuW,mCW,mWC,T1,T2,p),(n1,n2),pts):
    """
    nuC1x : Size of the cultivated population during bottleneck, gene category X
    nuC1y : Size of the cultivated population during bottleneck, gene category Y
    nuW   : Size of the wild population
    mCW	  : Migration rate from wild to cultivated
    mWC   : Migration rate from cultivated to wild
    T1	  : Duration of the bottleneck
    T2	  : Time since bottleneck 
    p	  : Proportion of genes under category Y
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
    return fs
    
    
    
    
#### Split and population growth in the cultivated population. No migration ####
	
def growth((nuCi,nuCf,nuW,T,y,p),(n1,n2),pts):
    """
    nuCi : Initial size of the cultivated population 
    nuCf : Final size of the cultivated population
    nuW  : Size of the wild population
    T	 : Total time since the split
    y    : Reduction of Ne of selected genes, category Y
    p    : Proportion of selected genes
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # CATX
    # split
    phiX = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # bottleneck and growth
    nuC_func = lambda t: nuCi*(nuCf/nuCi)**(t/T)
    phiX = dadi.Integration.two_pops(phiX, xx, T, nu1=nuC_func, nu2=nuW, m12=0, m21=0)
    # computing fsX
    fsX = dadi.Spectrum.from_phi(phiX, (n1,n2), (xx,xx))
    # CATY
    phiY = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # bottleneck and growth
    nuCiY = nuCi*y
    nuCfY = nuCf*y
    nuC_func = lambda t: nuCiY*(nuCfY/nuCiY)**(t/T)
    phiY = dadi.Integration.two_pops(phiY, xx, T, nu1=nuC_func, nu2=nuW, m12=0, m21=0)
    # computing fsX
    fsY = dadi.Spectrum.from_phi(phiY, (n1,n2), (xx,xx))
    # complete spectrum
    fs = (1-p)*fsX + p*fsY
    return fs


#### Split and population growth in the cultivated population. Symetrical migration ####
	
def growth_migration_sym((nuCi,nuCf,nuW,m,T,y,p),(n1,n2),pts):
    """
    nuCi : Initial size of the cultivated population 
    nuCf : Final size of the cultivated population
    nuW  : Size of the wild population
    m    : Migration rate
    T	 : Total time since the split
    y    : Reduction of Ne of selected genes, category Y
    p    : Proportion of selected genes
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # CATX
    # split
    phiX = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # bottleneck and growth
    nuC_func = lambda t: nuCi*(nuCf/nuCi)**(t/T)
    phiX = dadi.Integration.two_pops(phiX, xx, T, nu1=nuC_func, nu2=nuW, m12=m, m21=m)
    # computing fsX
    fsX = dadi.Spectrum.from_phi(phiX, (n1,n2), (xx,xx))
    # CATY
    phiY = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # bottleneck and growth
    nuCiY = nuCi*y
    nuCfY = nuCf*y
    nuC_func = lambda t: nuCiY*(nuCfY/nuCiY)**(t/T)
    phiY = dadi.Integration.two_pops(phiY, xx, T, nu1=nuC_func, nu2=nuW, m12=m, m21=m)
    # computing fsX
    fsY = dadi.Spectrum.from_phi(phiY, (n1,n2), (xx,xx))
    # complete spectrum
    fs = (1-p)*fsX + p*fsY
    return fs

#### Split and population growth in the cultivated population. Asymetrical migration ####


def growth_migration((nuCi,nuCf,nuW,mCW,mWC,T,y,p),(n1,n2),pts):
    """
    nuCi : Initial size of the cultivated population 
    nuCf : Final size of the cultivated population
    nuW  : Size of the wild population
    m    : Migration rate
    T	 : Total time since the split
    y    : Reduction of Ne of selected genes, category Y
    p    : Proportion of selected genes
    """
    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # CATX
    # split
    phiX = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # bottleneck and growth
    nuC_func = lambda t: nuCi*(nuCf/nuCi)**(t/T)
    phiX = dadi.Integration.two_pops(phiX, xx, T, nu1=nuC_func, nu2=nuW, m12=mCW, m21=mWC)
    # computing fsX
    fsX = dadi.Spectrum.from_phi(phiX, (n1,n2), (xx,xx))
    # CATY
    phiY = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # bottleneck and growth
    nuCiY = nuCi*y
    nuCfY = nuCf*y
    nuC_func = lambda t: nuCiY*(nuCfY/nuCiY)**(t/T)
    phiY = dadi.Integration.two_pops(phiY, xx, T, nu1=nuC_func, nu2=nuW, m12=mCW, m21=mWC)
    # computing fsX
    fsY = dadi.Spectrum.from_phi(phiY, (n1,n2), (xx,xx))
    # complete spectrum
    fs = (1-p)*fsX + p*fsY
    return fs
