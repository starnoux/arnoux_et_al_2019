#!/usr/local/bin/python
# -*- coding: utf-8 -*-
import sys
import dadi
import numpy

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
@date: 2020/05

@authors: ste.arnoux@gmail.com, christelle.fraisse.rios@gmail.com

@readme: a python module to define the demographic models tested with dadi.


In all models, an initial ancestral population is split in two sister populations (Crop and Wild).

(1) Timing of gene flow
	"SI": Strict Isolation: divergence occurs without gene flow between Crop and Wild.
	"IM": Isolation with Migration: divergence occurs with continuous gene flow between Crop and Wild.
	"IM2": Isolation with Migration 2: divergence occurs with continuous gene flow between Crop and Wild at different rates before and after the bottleneck.

(2) Demography
	"C": One-time period: after the split, Crop and Wild have a constant population size equals to the ancestral size.
	"BcCw": One-time period: after the split, Crop undergoes a bottleneck and stays at that size, Wild has a constant population size.
	"BcCw_E": Two-time period: after the split, Crop goes through a bottleneck for some time; then Crop and Wild change size exponentially.
	"C_BcCw_E": Three-time period: after the split, Crop and Wild have a constant population size; then Crop goes through a bottleneck for some time;
		    then Crop and Wild change size exponentially.
	"E": One-time period: after the split, Crop and Wild change size exponentially.
	"E_E": Two-time period: after the split, Crop and Wild change size exponentially twice.
	"C_E_E": Three-time period: after the split, Crop and Wild have a constant population size; then Crop and Wild change size exponentially twice.

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
									Demographic models
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def SI_C(params, (n1,n2), pts):

    # Define the parameters of the model
    Ts, O = params
    """""""""""""""""""""""""""""""""""""""
    Ts: Duration of divergence in isolation with constant sizes (phase 1).
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """""""""""""""""""""""""""""""""""""""

    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # Equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)

    # Split event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    # Integrate by setting divergence time to Ts, the population sizes to 1*nref and 1*nref and the migration rates to zero
    phi = dadi.Integration.two_pops(phi, xx, T=Ts, nu1=1, nu2=1, m12=0, m21=0)

    # Correctly oriented spectrum
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))

    # Misoriented spectrum
    fsM = dadi.Numerics.reverse_array(fsO)

    # Final spectrum
    fs = O*fsO + (1-O)*fsM
    return fs


def IM_C(params, (n1,n2), pts):
    mCW, mWC, Ts, O = params
    """""""""""""""""""""""""""""""""""""""
    mCW: Migration from Wild to Crop.
    mWC: Migration from Crop to Wild.
    Ts: Duration of divergence in continuous migration with constant sizes (phase 1).
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """""""""""""""""""""""""""""""""""""""
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T=Ts, nu1=1, nu2=1, m12=mCW, m21=mWC)
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    fsM = dadi.Numerics.reverse_array(fsO)
    fs = O*fsO + (1-O)*fsM
    return fs


def IM_BcCw(params, (n1,n2), pts):
    nuCb, mCW, mWC, Ts, O = params
    """""""""""""""""""""""""""""""""""""""
    nuCb: Size of Crop after split & bottleneck.
    mCW: Migration from Wild to Crop.
    mWC: Migration from Crop to Wild.
    Ts: Duration of divergence in continuous migration during which Crop goes through a bottleneck and Wild stays constant (phase 1).
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """""""""""""""""""""""""""""""""""""""
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T=Ts, nu1=nuCb, nu2=1, m12=mCW, m21=mWC)
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    fsM = dadi.Numerics.reverse_array(fsO)
    fs = O*fsO + (1-O)*fsM
    return fs


def IM_BcCw_E(params, (n1,n2), pts):
    nuCb, nuC, nuW, mCW, mWC, Ts, Te, O = params
    """""""""""""""""""""""""""""""""""""""
    nuCb: Size of Crop after split & bottleneck.
    nuC: Size of Crop at the end of exponential change.
    nuW: Size of Wild at the end of exponential change.
    mCW: Migration from Wild to Crop.
    mWC: Migration from Crop to Wild.
    Ts: Duration of divergence in continuous migration during which Crop goes through a bottleneck and Wild stays constant (phase 1).
    Te: Duration of divergence in continuous migration during which Crop and Wild change size exponentially (phase 2).
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """""""""""""""""""""""""""""""""""""""
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T=Ts, nu1=nuCb, nu2=1, m12=mCW, m21=mWC)
    nuC_func = lambda t: nuCb*(nuC/nuCb) ** (t/Te)
    nuW_func = lambda t: 1*(nuW/1) ** (t/Te)
    phi = dadi.Integration.two_pops(phi, xx, T=Te, nu1=nuC_func, nu2=nuW_func, m12=mCW, m21=mWC)
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    fsM = dadi.Numerics.reverse_array(fsO)
    fs = O*fsO + (1-O)*fsM
    return fs


def IM2_BcCw_E(params, (n1,n2), pts):
    nuCb, nuC, nuW, mCW, mWC, mCW2, mWC2, Ts, Te, O = params
    """""""""""""""""""""""""""""""""""""""
    nuCb: Size of Crop after split & bottleneck.
    nuC: Size of Crop at the end of exponential change.
    nuW: Size of Wild at the end of exponential change.
    mCW: Migration from Wild to Crop before bottleneck.
    mWC: Migration from Crop to Wild before bottleneck.
    mCW2: Migration from Wild to Crop after bottleneck.
    mWC2: Migration from Crop to Wild after bottleneck.
    Ts: Duration of divergence in continuous migration during which Crop goes through a bottleneck and Wild stays constant (phase 1).
    Te: Duration of divergence in continuous migration during which Crop and Wild change size exponentially (phase 2).
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """""""""""""""""""""""""""""""""""""""
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T=Ts, nu1=nuCb, nu2=1, m12=mCW, m21=mWC)
    nuC_func = lambda t: nuCb*(nuC/nuCb) ** (t/Te)
    nuW_func = lambda t: 1*(nuW/1) ** (t/Te)
    phi = dadi.Integration.two_pops(phi, xx, T=Te, nu1=nuC_func, nu2=nuW_func, m12=mCW2, m21=mWC2)
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    fsM = dadi.Numerics.reverse_array(fsO)
    fs = O*fsO + (1-O)*fsM
    return fs


def IM_C_BcCw_E(params, (n1,n2), pts):
    nuCb, nuC, nuW, mCW, mWC, Ts, Tb, Te, O = params
    """""""""""""""""""""""""""""""""""""""
    nuCb: Size of Crop after split & bottleneck.
    nuC: Size of Crop at the end of exponential change.
    nuW: Size of Wild at the end of exponential change.
    mCW: Migration from Wild to Crop.
    mWC: Migration from Crop to Wild.
    Ts: Duration of divergence in continuous migration with constant sizes (phase 1).
    Tb: Duration of divergence in continuous migration during which Crop goes through a bottleneck and Wild stays constant (phase 2).
    Te: Duration of divergence in continuous migration during which Crop and Wild change size exponentially (phase 3).
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """""""""""""""""""""""""""""""""""""""
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T=Ts, nu1=1, nu2=1, m12=mCW, m21=mWC)
    phi = dadi.Integration.two_pops(phi, xx, T=Tb, nu1=nuCb, nu2=1, m12=mCW, m21=mWC)
    nuC_func = lambda t: nuCb*(nuC/nuCb) ** (t/Te)
    nuW_func = lambda t: 1*(nuW/1) ** (t/Te)
    phi = dadi.Integration.two_pops(phi, xx, T=Te, nu1=nuC_func, nu2=nuW_func, m12=mCW, m21=mWC) 
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    fsM = dadi.Numerics.reverse_array(fsO)
    fs = O*fsO + (1-O)*fsM
    return fs


def IM2_C_BcCw_E(params, (n1,n2), pts):
    nuCb, nuC, nuW, mCW, mWC, mCW2, mWC2, Ts, Tb, Te, O = params
    """""""""""""""""""""""""""""""""""""""
    nuCb: Size of Crop after split & bottleneck.
    nuC: Size of Crop at the end of exponential change.
    nuW: Size of Wild at the end of exponential change.
    mCW: Migration from Wild to Crop before bottleneck.
    mWC: Migration from Crop to Wild before bottleneck.
    mCW2: Migration from Wild to Crop after bottleneck.
    mWC2: Migration from Crop to Wild after bottleneck.
    Ts: Duration of divergence in continuous migration with constant sizes (phase 1).
    Tb: Duration of divergence in continuous migration during which Crop goes through a bottleneck and Wild stays constant (phase 2).
    Te: Duration of divergence in continuous migration during which Crop and Wild change size exponentially (phase 3).
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """""""""""""""""""""""""""""""""""""""
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T=Ts, nu1=1, nu2=1, m12=mCW, m21=mWC)
    phi = dadi.Integration.two_pops(phi, xx, T=Tb, nu1=nuCb, nu2=1, m12=mCW, m21=mWC)
    nuC_func = lambda t: nuCb*(nuC/nuCb) ** (t/Te)
    nuW_func = lambda t: 1*(nuW/1) ** (t/Te)
    phi = dadi.Integration.two_pops(phi, xx, T=Te, nu1=nuC_func, nu2=nuW_func, m12=mCW2, m21=mWC2) 
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    fsM = dadi.Numerics.reverse_array(fsO)
    fs = O*fsO + (1-O)*fsM
    return fs


def IM_E(params, (n1,n2), pts):
    nuC, nuW, mCW, mWC, Ts, O = params
    """""""""""""""""""""""""""""""""""""""
    nuC: Size of Crop at the end of exponential change.
    nuW: Size of Wild at the end of exponential change.
    mCW: Migration from Wild to Crop.
    mWC: Migration from Crop to Wild.
    Ts: Duration of divergence in continuous migration during which Crop and Wild change size exponentially (phase 1).
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """""""""""""""""""""""""""""""""""""""
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    nuC_func = lambda t: 1*(nuC/1) ** (t/Ts)
    nuW_func = lambda t: 1*(nuW/1) ** (t/Ts)
    phi = dadi.Integration.two_pops(phi, xx, T=Ts, nu1=nuC_func, nu2=nuW_func, m12=mCW, m21=mWC)
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    fsM = dadi.Numerics.reverse_array(fsO)
    fs = O*fsO + (1-O)*fsM
    return fs


def IM_E_E(params, (n1,n2), pts):
    nuCe, nuC, nuWe, nuW, mCW, mWC, Ts, Te, O = params
    """""""""""""""""""""""""""""""""""""""
    nuCe: Size of Crop at the end of exponential change (phase 1).
    nuC: Size of Crop at the end of exponential change (phase 2).
    nuWe: Size of Wild at the end of exponential change (phase 1).
    nuW: Size of Wild at the end of exponential change (phase 2).
    mCW: Migration from Wild to Crop.
    mWC: Migration from Crop to Wild.
    Ts: Duration of divergence in continuous migration during which Crop and Wild change size exponentially (phase 1).
    Te: Duration of divergence in continuous migration during which Crop and Wild change size exponentially (phase 2).
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """""""""""""""""""""""""""""""""""""""
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    nuC_func = lambda t: 1*(nuCe/1) ** (t/Ts)
    nuW_func = lambda t: 1*(nuWe/1) ** (t/Ts)
    phi = dadi.Integration.two_pops(phi, xx, T=Ts, nu1=nuC_func, nu2=nuW_func, m12=mCW, m21=mWC)
    nuC_func = lambda t: nuCe*(nuC/nuCe) ** (t/Te)
    nuW_func = lambda t: nuWe*(nuW/nuWe) ** (t/Te)
    phi = dadi.Integration.two_pops(phi, xx, T=Te, nu1=nuC_func, nu2=nuW_func, m12=mCW, m21=mWC) 
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    fsM = dadi.Numerics.reverse_array(fsO)
    fs = O*fsO + (1-O)*fsM
    return fs


def IM2_E_E(params, (n1,n2), pts):
    nuCe, nuC, nuWe, nuW, mCW, mWC, mCW2, mWC2, Ts, Te, O = params
    """""""""""""""""""""""""""""""""""""""
    nuCe: Size of Crop at the end of exponential change (phase 1).
    nuC: Size of Crop at the end of exponential change (phase 2).
    nuWe: Size of Wild at the end of exponential change (phase 1).
    nuW: Size of Wild at the end of exponential change (phase 2).
    mCW: Migration from Wild to Crop.
    mWC: Migration from Crop to Wild.
    mCW2: Migration from Wild to Crop after phase 1.
    mWC2: Migration from Crop to Wild after phase 1.
    Ts: Duration of divergence in continuous migration during which Crop and Wild change size exponentially (phase 1).
    Te: Duration of divergence in continuous migration during which Crop and Wild change size exponentially (phase 2).
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """""""""""""""""""""""""""""""""""""""
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    nuC_func = lambda t: 1*(nuCe/1) ** (t/Ts)
    nuW_func = lambda t: 1*(nuWe/1) ** (t/Ts)
    phi = dadi.Integration.two_pops(phi, xx, T=Ts, nu1=nuC_func, nu2=nuW_func, m12=mCW, m21=mWC)
    nuC_func = lambda t: nuCe*(nuC/nuCe) ** (t/Te)
    nuW_func = lambda t: nuWe*(nuW/nuWe) ** (t/Te)
    phi = dadi.Integration.two_pops(phi, xx, T=Te, nu1=nuC_func, nu2=nuW_func, m12=mCW2, m21=mWC2) 
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    fsM = dadi.Numerics.reverse_array(fsO)
    fs = O*fsO + (1-O)*fsM
    return fs


def IM_C_E_E(params, (n1,n2), pts):
    nuCe, nuC, nuWe, nuW, mCW, mWC, Ts, Tb, Te, O = params
    """""""""""""""""""""""""""""""""""""""
    nuCe: Size of Crop at the end of exponential change (phase 2).
    nuC: Size of Crop at the end of exponential change (phase 3).
    nuWe: Size of Wild at the end of exponential change (phase 2).
    nuW: Size of Wild at the end of exponential change (phase 3).
    mCW: Migration from Wild to Crop.
    mWC: Migration from Crop to Wild.
    Ts: Duration of divergence in continuous migration with constant sizes (phase 1).
    Tb: Duration of divergence in continuous migration during which Crop and Wild change size exponentially (phase 2).
    Te: Duration of divergence in continuous migration during which Crop and Wild change size exponentially (phase 3).
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """""""""""""""""""""""""""""""""""""""
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T=Ts, nu1=1, nu2=1, m12=mCW, m21=mWC)
    nuC_func = lambda t: 1*(nuCe/1) ** (t/Tb)
    nuW_func = lambda t: 1*(nuWe/1) ** (t/Tb)
    phi = dadi.Integration.two_pops(phi, xx, T=Tb, nu1=nuC_func, nu2=nuW_func, m12=mCW, m21=mWC)
    nuC_func = lambda t: nuCe*(nuC/nuCe) ** (t/Te)
    nuW_func = lambda t: nuWe*(nuW/nuWe) ** (t/Te)
    phi = dadi.Integration.two_pops(phi, xx, T=Te, nu1=nuC_func, nu2=nuW_func, m12=mCW, m21=mWC) 
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    fsM = dadi.Numerics.reverse_array(fsO)
    fs = O*fsO + (1-O)*fsM
    return fs


def IM2_C_E_E(params, (n1,n2), pts):
    nuCe, nuC, nuWe, nuW, mCW, mWC, mCW2, mWC2, Ts, Tb, Te, O = params
    """""""""""""""""""""""""""""""""""""""
    nuCe: Size of Crop at the end of exponential change (phase 2).
    nuC: Size of Crop at the end of exponential change (phase 3).
    nuWe: Size of Wild at the end of exponential change (phase 2).
    nuW: Size of Wild at the end of exponential change (phase 3).
    mCW: Migration from Wild to Crop.
    mWC: Migration from Crop to Wild.
    mCW2: Migration from Wild to Crop after phase 2.
    mWC2: Migration from Crop to Wild after phase 2.
    Ts: Duration of divergence in continuous migration with constant sizes (phase 1).
    Tb: Duration of divergence in continuous migration during which Crop and Wild change size exponentially (phase 2).
    Te: Duration of divergence in continuous migration during which Crop and Wild change size exponentially (phase 3).
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """""""""""""""""""""""""""""""""""""""
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T=Ts, nu1=1, nu2=1, m12=mCW, m21=mWC)
    nuC_func = lambda t: 1*(nuCe/1) ** (t/Tb)
    nuW_func = lambda t: 1*(nuWe/1) ** (t/Tb)
    phi = dadi.Integration.two_pops(phi, xx, T=Tb, nu1=nuC_func, nu2=nuW_func, m12=mCW, m21=mWC)
    nuC_func = lambda t: nuCe*(nuC/nuCe) ** (t/Te)
    nuW_func = lambda t: nuWe*(nuW/nuWe) ** (t/Te)
    phi = dadi.Integration.two_pops(phi, xx, T=Te, nu1=nuC_func, nu2=nuW_func, m12=mCW2, m21=mWC2) 
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    fsM = dadi.Numerics.reverse_array(fsO)
    fs = O*fsO + (1-O)*fsM
    return fs

