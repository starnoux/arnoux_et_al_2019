#!/usr/local/bin/python
# -*- coding: utf-8 -*-



import sys

import dadi



sys.path.append('python')

import numpy



"""

Defines the different models of divergence tested with DaDi. In all models, an initial ancestral population is split in two sister populations.



(1) Timing of gene flow 

	"SI": Strict Isolation: divergence occurs without gene flow between the two sister species

	"IM": Isolation with Migration: divergence occurs with continuous gene flow between the two sister species

        "IM2": Isolation with Migration: divergence occurs with two seperated continuous gene flow between the two sister species before and after bottleneck



(2) Demography

	"C": At the time of split, Crop and Wild have a constant population size equals to the ancestral size

	"BcCw": At the time of split, Crop undergoes a bottleneck and stays at that size, Wild has a constant population size

	"BcCw_E": At the time of split, Crop goes through a bottleneck for some time; then Crop and Wild grow/decline exponentially

	"C_BcCw_E": At the time of split, Crop and Wild have a constant population size; then Crop goes through a bottleneck for some time; then Crop and Wild grow/decline exponentially

	"E": At the time of split, Crop and Wild grows/declines exponentially

	"E_E": At the time of split, Crop and Wild grow/decline exponentially twice

	"C_E_E": At the time of split, Crop and Wild have a constant population size; then Crop and Wild grow/decline exponentially twice


"""









"""""""""""""""""""""""""""""""""""""""

Strict Isolation models (SI)

"""""""""""""""""""""""""""""""""""""""



def SI_C(params, (n1,n2), pts):

    Ts, O = params

    """""""""""""""""""""""""""""""""""""""

    Ts: Time of split in strict isolation with constant population size (in units of 2*Na generations).

    O: Proportion of accurate SNP orientation.



    n1,n2: Size of fs to generate.

    pts: Number of points to use in grid for evaluation.

    """""""""""""""""""""""""""""""""""""""

    # Define the grid we'll use

    xx = dadi.Numerics.default_grid(pts)

    # Ancestral population

    phi = dadi.PhiManip.phi_1D(xx)

    # Split event

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    # Divergence in strict isolation: crop and wild stay constant

    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1=1, nu2=1, m12=0, m21=0)

    # Correctly oriented spectrum

    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))

    # Misoriented spectrum

    fsM = dadi.Numerics.reverse_array(fsO)



    # Final spectrum

    fs = O*fsO+(1-O)*fsM

    return fs



""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""



"""""""""""""""""""""""""""""""""""""""

Isolation with Migration models (IM)

"""""""""""""""""""""""""""""""""""""""



def IM_C(params, (n1,n2), pts):

    mCW, mWC, Ts, O = params

    """""""""""""""""""""""""""""""""""""""

    Ts: Time of split in continuous migration (in units of 2*Na generations).

    mCW: Migration from Wild population to Crop population.

    mWC: Migration from Crop population to Wild population.

    O: Proportion of accurate SNP orientation.



    n1,n2: Size of fs to generate.

    pts: Number of points to use in grid for evaluation.

    """""""""""""""""""""""""""""""""""""""

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1=1, nu2=1, m12=mCW, m21=mWC)

    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))

    fsM = dadi.Numerics.reverse_array(fsO)



    fs = O*fsO+(1-O)*fsM

    return fs





def IM_BcCw(params, (n1,n2), pts):

    nuCb, mCW, mWC, Ts, O = params

    """""""""""""""""""""""""""""""""""""""

    nuCb: Size of Crop bottlenecked population after split (relative to Na).

    Ts: Time of split in continuous migration during which Crop goes through a bottleneck and Wild stays constant (in units of 2*Na generations).

    mCW: Migration from Wild population to Crop population.

    mWC: Migration from Crop population to Wild population.

    O: Proportion of accurate SNP orientation.



    n1,n2: Size of fs to generate.

    pts: Number of points to use in grid for evaluation.

    """""""""""""""""""""""""""""""""""""""

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, Ts, nuCb, 1, mCW, mWC)

    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))

    fsM = dadi.Numerics.reverse_array(fsO)



    fs = O*fsO+(1-O)*fsM

    return fs




def IM_BcCw_E(params, (n1,n2), pts):

    nuCb, nuC, nuW, mCW, mWC, Ts, Te, O = params

    """""""""""""""""""""""""""""""""""""""

    nuCb: Size of Crop bottlenecked population after split (relative to Na).

    nuC: Size of Crop population after exponential size change (relative to Na).

    nuW: Size of Wild population after exponential size change (relative to Na).

    Ts: Time of split in continuous migration during which Crop undergoes a bottleneck (in units of 2*Na generations).

    Te: Time in continuous migration during which Crop and Wild grow/decline exponentially (in units of 2*Na generations).

    mCW: Migration from Wild population to Crop population.

    mWC: Migration from Crop population to Wild population.

    O: Proportion of accurate SNP orientation.



    n1,n2: Size of fs to generate.

    pts: Number of points to use in grid for evaluation.

    """""""""""""""""""""""""""""""""""""""

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, Ts, nuCb, 1, mCW, mWC)

    nuC_func = lambda t: nuCb*(nuC/nuCb) ** (t/Te)

    nuW_func = lambda t: 1*(nuW/1) ** (t/Te)

    phi = dadi.Integration.two_pops(phi, xx, Te, nuC_func, nuW_func, mCW, mWC)

    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))

    fsM = dadi.Numerics.reverse_array(fsO)



    fs = O*fsO+(1-O)*fsM

    return fs




def IM2_BcCw_E(params, (n1,n2), pts):

    nuCb, nuC, nuW, mCW, mWC, mCW2, mWC2, Ts, Te, O = params

    """""""""""""""""""""""""""""""""""""""

    nuCb: Size of Crop bottlenecked population after split (relative to Na).

    nuC: Size of Crop population after exponential size change (relative to Na).

    nuW: Size of Wild population after exponential size change (relative to Na).

    Ts: Time of split in continuous migration during which Crop undergoes a bottleneck (in units of 2*Na generations).

    Te: Time in continuous migration during which Crop and Wild grow/decline exponentially (in units of 2*Na generations).

    mCW: Migration from Wild population to Crop population.

    mWC: Migration from Crop population to Wild population.

    mCW2: Migration from Wild population to Crop population After bottleneck.

    mWC2: Migration from Crop population to Wild population After bottleneck.

    O: Proportion of accurate SNP orientation.



    n1,n2: Size of fs to generate.

    pts: Number of points to use in grid for evaluation.

    """""""""""""""""""""""""""""""""""""""

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, Ts, nuCb, 1, mCW, mWC)

    nuC_func = lambda t: nuCb*(nuC/nuCb) ** (t/Te)

    nuW_func = lambda t: 1*(nuW/1) ** (t/Te)

    phi = dadi.Integration.two_pops(phi, xx, Te, nuC_func, nuW_func, mCW2, mWC2)

    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))

    fsM = dadi.Numerics.reverse_array(fsO)



    fs = O*fsO+(1-O)*fsM

    return fs




def IM_C_BcCw_E(params, (n1,n2), pts):

    nuCb, nuC, nuW, mCW, mWC, Ts, Tb, Te, O = params

    """""""""""""""""""""""""""""""""""""""

    nuCb: Size of Crop bottlenecked population after split (relative to Na).

    nuC: Size of Crop population after exponential size change (relative to Na).

    nuW: Size of Wild population after exponential size change (relative to Na).

    Ts: Time of split in continuous migration with constant population size (in units of 2*Na generations).

    Tb: Time in continuous migration during which Crop undergoes a bottleneck (in units of 2*Na generations).

    Te: Time in continuous migration during which Crop and Wild grow/decline exponentially (in units of 2*Na generations).

    mCW: Migration from Wild population to Crop population.

    mWC: Migration from Crop population to Wild population.

    O: Proportion of accurate SNP orientation.



    n1,n2: Size of fs to generate.

    pts: Number of points to use in grid for evaluation.

    """""""""""""""""""""""""""""""""""""""

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, Ts, 1, 1, mCW, mWC)

    phi = dadi.Integration.two_pops(phi, xx, Tb, nuCb, 1, mCW, mWC)

    nuC_func = lambda t: nuCb*(nuC/nuCb) ** (t/Te)

    nuW_func = lambda t: 1*(nuW/1) ** (t/Te)

    phi = dadi.Integration.two_pops(phi, xx, Te, nuC_func, nuW_func, mCW, mWC) 

    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))

    fsM = dadi.Numerics.reverse_array(fsO)



    fs = O*fsO+(1-O)*fsM

    return fs




def IM2_C_BcCw_E(params, (n1,n2), pts):

    nuCb, nuC, nuW, mCW, mWC, mCW2, mWC2, Ts, Tb, Te, O = params

    """""""""""""""""""""""""""""""""""""""

    nuCb: Size of Crop bottlenecked population after split (relative to Na).

    nuC: Size of Crop population after exponential size change (relative to Na).

    nuW: Size of Wild population after exponential size change (relative to Na).

    Ts: Time of split in continuous migration with constant population size (in units of 2*Na generations).

    Tb: Time in continuous migration during which Crop undergoes a bottleneck (in units of 2*Na generations).

    Te: Time in continuous migration during which Crop and Wild grow/decline exponentially (in units of 2*Na generations).

    mCW: Migration from Wild population to Crop population.

    mWC: Migration from Crop population to Wild population.

    mCW2: Migration from Wild population to Crop population After bottleneck.

    mWC2: Migration from Crop population to Wild population After bottleneck.

    O: Proportion of accurate SNP orientation.



    n1,n2: Size of fs to generate.

    pts: Number of points to use in grid for evaluation.

    """""""""""""""""""""""""""""""""""""""

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, Ts, 1, 1, mCW, mWC )

    phi = dadi.Integration.two_pops(phi, xx, Tb, nuCb, 1, mCW, mWC )

    nuC_func = lambda t: nuCb*(nuC/nuCb) ** (t/Te)

    nuW_func = lambda t: 1*(nuW/1) ** (t/Te)

    phi = dadi.Integration.two_pops(phi, xx, Te, nuC_func, nuW_func, mCW2, mWC2 ) 

    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))

    fsM = dadi.Numerics.reverse_array(fsO)



    fs = O*fsO+(1-O)*fsM

    return fs




def IM_E(params, (n1,n2), pts):

    nuC, nuW, mCW, mWC, Ts, O = params

    """""""""""""""""""""""""""""""""""""""

    nuC: Size of Crop population after exponential size change (relative to Na).

    nuW: Size of Wild population after exponential size change (relative to Na).

    Ts: Time of split in continuous migration during which Crop and Wild grow/decline exponentially (in units of 2*Na generations).

    mCW: Migration from Wild population to Crop population.

    mWC: Migration from Crop population to Wild population.

    O: Proportion of accurate SNP orientation.


    n1,n2: Size of fs to generate.

    pts: Number of points to use in grid for evaluation.

    """""""""""""""""""""""""""""""""""""""

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nuC_func = lambda t: 1*(nuC/1) ** (t/Ts)

    nuW_func = lambda t: 1*(nuW/1) ** (t/Ts)

    phi = dadi.Integration.two_pops(phi, xx, Ts, nuC_func, nuW_func, mCW, mWC )

    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))

    fsM = dadi.Numerics.reverse_array(fsO)


    fs = O*fsO+(1-O)*fsM

    return fs




def IM_E_E(params, (n1,n2), pts):

    nuCe, nuC, nuWe, nuW, mCW, mWC, Ts, Te, O = params

    """""""""""""""""""""""""""""""""""""""

    nuCe: Size of Crop population after split (relative to Na).

    nuC: Size of Crop population after exponential size change (relative to Na).

    nuWe: Size of Crop population after split (relative to Na).

    nuW: Size of Wild population after exponential size change (relative to Na).

    Ts: Time of split in continuous migration during which Crop and Wild grow/decline exponentially (in units of 2*Na generations).

    Te: Time in continuous migration during which Crop and Wild grow/decline exponentially (in units of 2*Na generations).

    mCW: Migration from Wild population to Crop population.

    mWC: Migration from Crop population to Wild population.

    O: Proportion of accurate SNP orientation.


    n1,n2: Size of fs to generate.

    pts: Number of points to use in grid for evaluation.

    """""""""""""""""""""""""""""""""""""""

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nuC_func = lambda t: 1*(nuCe/1) ** (t/Ts)

    nuW_func = lambda t: 1*(nuWe/1) ** (t/Ts)

    phi = dadi.Integration.two_pops(phi, xx, Ts, nuC_func, nuW_func, mCW, mWC )

    nuC_func = lambda t: nuCe*(nuC/nuCe) ** (t/Te)

    nuW_func = lambda t: nuWe*(nuW/nuWe) ** (t/Te)

    phi = dadi.Integration.two_pops(phi, xx, Te, nuC_func, nuW_func, mCW, mWC ) 

    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))

    fsM = dadi.Numerics.reverse_array(fsO)


    fs = O*fsO+(1-O)*fsM

    return fs




def IM_C_E_E(params, (n1,n2), pts):

    nuCe, nuC, nuWe, nuW, mCW, mWC, Ts, Tb, Te, O = params

    """""""""""""""""""""""""""""""""""""""

    nuCe: Size of Crop population after exponential size change (relative to Na).

    nuC: Size of Crop population after exponential size change (relative to Na).

    nuWe: Size of Crop population after exponential size change (relative to Na).

    nuW: Size of Wild population after exponential size change (relative to Na).

    Ts: Time of split in continuous migration during which Crop and Wild stay constant (in units of 2*Na generations).

    Tb: Time in continuous migration during which Crop grows/declines exponentially (in units of 2*Na generations).

    Te: Time in continuous migration during which Crop grows/declines exponentially (in units of 2*Na generations).

    mCW: Migration from Wild population to Crop population.

    mWC: Migration from Crop population to Wild population.

    O: Proportion of accurate SNP orientation.


    n1,n2: Size of fs to generate.

    pts: Number of points to use in grid for evaluation.

    """""""""""""""""""""""""""""""""""""""

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, Ts, 1, 1, mCW, mWC )

    nuC_func = lambda t: 1*(nuCe/1) ** (t/Tb)

    nuW_func = lambda t: 1*(nuWe/1) ** (t/Tb)

    phi = dadi.Integration.two_pops(phi, xx, Tb, nuC_func, nuW_func, mCW, mWC )

    nuC_func = lambda t: nuCe*(nuC/nuCe) ** (t/Te)

    nuW_func = lambda t: nuWe*(nuW/nuWe) ** (t/Te)

    phi = dadi.Integration.two_pops(phi, xx, Te, nuC_func, nuW_func, mCW, mWC ) 

    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))

    fsM = dadi.Numerics.reverse_array(fsO)


    fs = O*fsO+(1-O)*fsM

    return fs


