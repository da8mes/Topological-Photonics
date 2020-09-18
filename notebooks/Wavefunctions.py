# see also integrands.py
# The choices of function and variable names will make much more sense if you
# access the latest ET_Redux_V*.lyx 

import numpy as np
from numba import jit

LOOKUP_TABLE = np.array([
    1, 1, 2, 6, 24, 120, 720, 5040, 40320,
    362880, 3628800, 39916800, 479001600,
    6227020800, 87178291200, 1307674368000,
    20922789888000, 355687428096000, 6402373705728000,
    121645100408832000, 2432902008176640000], dtype='int64')


@jit(nopython=True)
def fast_factorial(n):
    if n > 20:
        raise ValueError
    return LOOKUP_TABLE[n]


def wf_lll_radial(indices_tuple, r):
    # the radial part of a wavefunction in the lowest Landau level
    # includes normalization factors so that a 2D integral of np.abs(wf_lll_radial())**2 yields 1
    # inputs:
    #    indices_tuple:  a tuple with a single integer entry, the angular momentum 
    #    r:              the unitless radial coordinate
    # outputs:
    #    the complex value of the wavefunction
    l = indices_tuple[0]
    return wf_lll_radial_jit(l, r)

@jit(nopython=True)
def wf_lll_radial_jit(l, r):
    # the radial part of a wavefunction in the lowest Landau level
    # includes normalization factors so that a 2D integral of np.abs(wf_lll_radial())**2 yields 1
    # inputs:
    #    indices_tuple:  a tuple with a single integer entry, the angular momentum 
    #    r:              the unitless radial coordinate
    # outputs:
    #    the complex value of the wavefunction
    return r**l * np.exp(-r**2/4) / np.sqrt(2*np.pi * 2**l * fast_factorial(l))
    
    
    