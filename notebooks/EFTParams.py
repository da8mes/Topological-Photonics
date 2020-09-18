# The choices of function and variable names will make much more sense if you
# access the latest ET_Redux_V*.lyx 

import Integrands
import Wavefunctions
from scipy import integrate
import numpy as np


# IDEALLY DON'T USE THIS; USE THE PRECOMPUTE FUNCTIONALITY 
def dOmega_parameter(wf, mode_indices, Omega, H, epsabs=1.49e-08, epsrel=1.49e-08):
    print("WARNING: Not yet implemented!!")
    
    return None

    
    
# IDEALLY DON'T USE THIS; USE THE PRECOMPUTE FUNCTIONALITY 
def U_parameter(wf, mode_indices, C6, w, H, epsabs=1.49e-08, epsrel=1.49e-08):
    mi=mode_indices
    numerator   = C(wf, mi, H, epsabs=epsabs, epsrel=epsrel) 
    G1          = G_2mode(wf, mi[:2], H, epsabs=epsabs, epsrel=epsrel) 
    G2          = G_2mode(wf, mi[2:], H, epsabs=epsabs, epsrel=epsrel)
    denominator = 2 * G1 * G2 * np.sqrt((1+(mi[0][0]==mi[1][0]))*(1+(mi[2][0]==mi[3][0])))
    scale       = (1 + np.abs(H))**(2/3) * C6 * (w/2)**-6
    
    return scale*numerator/denominator

    
def G_2mode(wf, mode_indices, H, epsabs=1.49e-08, epsrel=1.49e-08):
    return np.sqrt(A(wf, mode_indices, H, epsabs=epsabs, epsrel=epsrel)/2)


def A(wf, mode_indices, H, epsabs=1.49e-08, epsrel=1.49e-08):
    # return complex value of the A integral defined in ET_Redux_V*.lyx 
    # ASSUMES that the first entry of each mode index tuple represents an angular momentum, and that the total of this angular momentum should be conserved, as in the Lowest Landau Level
    
    assert len(mode_indices) == 2
    
    result = integral(Integrands.a_integrand_jit, wf, mode_indices, H, epsabs=epsabs, epsrel=epsrel)
        
    return result




def B(wf, mode_indices, H, epsabs=1.49e-08, epsrel=1.49e-08):
    # return complex value of the B integral defined in ET_Redux_V*.lyx 
    # ASSUMES that the first entry of each mode index tuple represents an angular momentum, and that the total of this angular momentum should be conserved, as in the Lowest Landau Level
    
    assert len(mode_indices) == 4
    
    # only need to integrate if total angular momentum is conserved
    if mode_indices[0][0] + mode_indices[1][0] == mode_indices[2][0] + mode_indices[3][0]:
        result = integral(Integrands.b_integrand_jit, wf, mode_indices, H, epsabs=epsabs, epsrel=epsrel)
    else:
        result = 0
        
    return result




def C(wf, mode_indices, H, epsabs=1.49e-08, epsrel=1.49e-08):
    # return complex value of the C integral defined in ET_Redux_V*.lyx 
    # ASSUMES that the first entry of each mode index tuple represents an angular momentum, and that the total of this angular momentum should be conserved, as in the Lowest Landau Level
    
    assert len(mode_indices) == 4
    
    # only need to integrate if total angular momentum is conserved
    if mode_indices[0][0] + mode_indices[1][0] == mode_indices[2][0] + mode_indices[3][0]:
        result = integral_c(Integrands.c_integrand_jit, wf, mode_indices, H, epsabs=epsabs, epsrel=epsrel)
    else:
        result = 0
        
    return result
    
    
    

def integral(integrand_func, wf, mode_indices, H, epsabs=1.49e-08, epsrel=1.49e-08):
    theta_start = 0
    theta_end   = 2*np.pi

    r2_start    = lambda x: 0
    r2_end      = lambda x: np.inf

    r1_start    = lambda x, y: 0
    r1_end      = lambda x, y: np.inf

    if len(mode_indices) == 4:
        l0 = mode_indices[0][0]
        l1 = mode_indices[1][0]
        l2 = mode_indices[2][0]
        l3 = mode_indices[3][0]
        integrand_real = lambda r1, r2, dtheta: np.real(integrand_func(wf, l0, l1, l2, l3, H, r1, r2, dtheta))
        integrand_imag = lambda r1, r2, dtheta: np.imag(integrand_func(wf, l0, l1, l2, l3, H, r1, r2, dtheta))
    elif len(mode_indices) == 2:
        l0 = mode_indices[0][0]
        l1 = mode_indices[1][0]
        integrand_real = lambda r1, r2, dtheta: np.real(integrand_func(wf, l0, l1, H, r1, r2, dtheta))
        integrand_imag = lambda r1, r2, dtheta: np.imag(integrand_func(wf, l0, l1, H, r1, r2, dtheta))
    else:
        print("Error badly formed mode_indices")
        return None

    ## TEST ##
    integral_real, error_bound_real = integrate.tplquad(integrand_real, theta_start, theta_end, r2_start, r2_end, r1_start, r1_end, epsabs=epsabs, epsrel=epsrel)
    integral_imag, error_bound_imag = integrate.tplquad(integrand_imag, theta_start, theta_end, r2_start, r2_end, r1_start, r1_end, epsabs=epsabs, epsrel=epsrel)

    complete_complex_result = (integral_real + 1j * integral_imag) # include that factor of 2pi!
  
    return complete_complex_result




def integral_c(integrand_func, wf, mode_indices, H, epsabs=1.49e-08, epsrel=1.49e-08):
    # use theta as the innermost index
    theta_start = lambda r1, dr: -np.pi
    theta_end   = lambda r1, dr: np.pi

    # NOTE: here we use dr instead of r2, so the integral starts at -r1
    dr_start    = lambda r1: -r1
    dr_end      = lambda r1: np.inf

    r1_start    = 0
    r1_end      = np.inf
    
    l0 = mode_indices[0][0]
    l1 = mode_indices[1][0]
    l2 = mode_indices[2][0]
    l3 = mode_indices[3][0]

    integrand_real = lambda dtheta, dr, r1: np.real(integrand_func(wf, l0, l1, l2, l3, H, r1, dr, dtheta))
    integrand_imag = lambda dtheta, dr, r1: np.imag(integrand_func(wf, l0, l1, l2, l3, H, r1, dr, dtheta))

    ## TEST ##
    integral_real, error_bound_real = integrate.tplquad(integrand_real, r1_start, r1_end, dr_start, dr_end, theta_start, theta_end, epsabs=epsabs, epsrel=epsrel)
    integral_imag, error_bound_imag = integrate.tplquad(integrand_imag, r1_start, r1_end, dr_start, dr_end, theta_start, theta_end, epsabs=epsabs, epsrel=epsrel)

    complete_complex_result = (integral_real + 1j * integral_imag) # include that factor of 2pi!
  
    return complete_complex_result

