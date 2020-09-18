# see also wavefunctions.py
# The choices of function and variable names will make much more sense if you
# access the latest ET_Redux_V*.lyx 

import numpy as np
from numba import jit

# integrand for the two-mode normalization factor 
# takes wavefunction as an additional input
@jit(nopython=True)
def a_integrand_jit(wf, l0, l1, H, r1, r2, dtheta):
    eps = 10**-40
    dist_factor = ((r1 * np.sin(dtheta))**2 + (r1 * np.cos(dtheta) - r2)**2 + eps)**-3 + eps
    wf_factor = wf(l0, r1)*wf(l1, r2)*np.exp(1j*(l0 - l1)*dtheta) + wf(l0, r2)*wf(l1, r1) 
    
    integrand = 2*np.pi* (1 + np.abs(H))**(2) * r1*r2 / np.abs(H + dist_factor)**2 * np.abs(wf_factor)**2 # 
    
    return integrand


def a_integrand(wf, mode_indices, H, r1, r2, dtheta):
    # inputs:
    #    wf:                  callable;  the radial wavefunction, should take a tuple as first input and radius as second input
    #    mode_indices:        a list of *two* tuples; each tuple can be passed to radial_wavefunction along with a radius to get a wavefunction value. tuple[0] should be the angular momentum
    #    H:                   the dimensionless parameter characterizing the interaction strength relative to the Rydberg lifetime
    #    r1:                  radial position of particle 1
    #    r2:                  radial position of particle 2
    #    dtheta:              difference in angular positions between the particles
    #
    # outputs:
    #    integrand:           the value of the integrand at the given location
    
    mi = mode_indices
    
    dist_factor = ((r1 * np.sin(dtheta))**2 + (r1 * np.cos(dtheta) - r2)**2)**-3
    wf_factor = wf(mi[0], r1)*wf(mi[1], r2)*np.exp(1j*(mi[0][0] - mi[1][0])*dtheta) + wf(mi[0], r2)*wf(mi[1], r1) 
    
    integrand = 2*np.pi* (1 + np.abs(H))**(2) * r1*r2 / np.abs(H + dist_factor)**2 * np.abs(wf_factor)**2 # 
    
    return integrand
    


# integrand for the omega parameter
@jit(nopython=True)
def b_integrand_jit(wf, l0, l1, l2, l3, H, r1, r2, dtheta):
    eps = 10**-40
    dist_factor = ((r1 * np.sin(dtheta))**2 + (r1 * np.cos(dtheta) - r2)**2 + eps)**-3 + eps
    wf_factor_1 = wf(l0, r1) * wf(l1, r2) * np.exp(1j * l0 * dtheta)  +  wf(l0, r2) * wf(l1, r1) * np.exp(1j * l1 * dtheta)
    wf_factor_2 = wf(l2, r1) * wf(l3, r2) * np.exp(1j * l2 * dtheta)
    
    integrand = 2*np.pi* (1 + np.abs(H)) * r1*r2 * np.conj(wf_factor_1 / (H + dist_factor)) * wf_factor_2 
    
    return integrand


def b_integrand(wf, mode_indices, H, r1, r2, dtheta):
    # inputs:
    #    wf:                  callable;  the radial wavefunction, should take a tuple as first input and radius as second input
    #    mode_indices:        a list of *four* tuples; each tuple can be passed to radial_wavefunction along with a radius to get a wavefunction value. tuple[0] should be the angular momentum
    #    H:                   the dimensionless parameter characterizing the interaction strength relative to the Rydberg lifetime
    #    r1:                  radial position of particle 1
    #    r2:                  radial position of particle 2
    #    dtheta:              difference in angular positions between the particles
    #
    # outputs:
    #    integrand:           the value of the integrand at the given location
    
    mi = mode_indices
    
    dist_factor = ((r1 * np.sin(dtheta))**2 + (r1 * np.cos(dtheta) - r2)**2)**-3
    wf_factor_1 = wf(mi[0], r1) * wf(mi[1], r2) * np.exp(1j * mi[0][0] * dtheta)  +  wf(mi[0], r2) * wf(mi[1], r1) * np.exp(1j * mi[1][0] * dtheta)
    wf_factor_2 = wf(mi[2], r1) * wf(mi[3], r2) * np.exp(1j * mi[2][0] * dtheta)
    
    integrand = 2*np.pi* (1 + np.abs(H)) * r1*r2 * np.conj(wf_factor_1 / (H + dist_factor)) * wf_factor_2 
    
    return integrand



# integrand for the interaction parameter
@jit(nopython=True)
def c_integrand_jit(wf, l0, l1, l2, l3, H, r1, dr, dtheta):
    eps = 10**-40
    r2 = r1 + dr
    dist_factor =  (2*r1*(r2)*(1-np.cos(dtheta)) + dr**2 + eps)**-3 + eps
    wf_factor_1 = np.conj( wf(l0, r1) * wf(l1, r2) * np.exp(1j * l0 * dtheta)  +  wf(l0, r2) * wf(l1, r1) * np.exp(1j * l1 * dtheta) )
    wf_factor_2 =          wf(l2, r1) * wf(l3, r2) * np.exp(1j * l2 * dtheta)  +  wf(l2, r2) * wf(l3, r1) * np.exp(1j * l3 * dtheta) 
    
    integrand = 2*np.pi* (1 + np.abs(H))**(4/3) * r1*r2 * dist_factor / np.abs(H + dist_factor)**2 * wf_factor_1 * wf_factor_2
    
    return integrand


def c_integrand(wf, mode_indices, H, r1, dr, dtheta, eps=0):
    # inputs:
    #    wf:                  callable;  the radial wavefunction, should take a tuple as first input and radius as second input
    #    mode_indices:        a list of *four* tuples; each tuple can be passed to radial_wavefunction along with a radius to get a wavefunction value. tuple[0] should be the angular momentum
    #    H:                   the dimensionless parameter characterizing the interaction strength relative to the Rydberg lifetime
    #    r1:                  radial position of particle 1
    #    r2:                  radial position of particle 2
    #    dtheta:              difference in angular positions between the particles
    #    eps:                 a hack for trying to correct divide by zero errors
    #
    # outputs:
    #    integrand:           the value of the integrand at the given location
    
    mi = mode_indices
    
    r2 = r1 + dr
    dist_factor =  (2*r1*(r2)*(1-np.cos(dtheta)) + dr**2 + eps)**-3
    wf_factor_1 = np.conj( wf(mi[0], r1) * wf(mi[1], r2) * np.exp(1j * mi[0][0] * dtheta)  +  wf(mi[0], r2) * wf(mi[1], r1) * np.exp(1j * mi[1][0] * dtheta) )
    wf_factor_2 =          wf(mi[2], r1) * wf(mi[3], r2) * np.exp(1j * mi[2][0] * dtheta)  +  wf(mi[2], r2) * wf(mi[3], r1) * np.exp(1j * mi[3][0] * dtheta) 
    
    integrand = 2*np.pi* (1 + np.abs(H))**(4/3) * r1*r2 * dist_factor / np.abs(H + dist_factor)**2 * wf_factor_1 * wf_factor_2
    
    return integrand