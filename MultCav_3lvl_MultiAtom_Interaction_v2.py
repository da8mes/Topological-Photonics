#!/usr/bin/env python
# coding: utf-8

# In[2]:


# 3 cavity modes, 3lvl, multiAtom, Interactions

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import chirp, find_peaks, peak_widths
from qutip import *
import time


# In[3]:


# Modeling multiple cavity modes with interaction term

cav_modes = 3 # num of cav modes
cav_max = 3 # this means can have 0, 1 or 2 photons
phot_init = 0 
#p_init = 0
#r_init = 0

# generic initial state
cav0 = basis(cav_max, n = phot_init)
p0 = basis(cav_max, n = phot_init)
r0 = basis(cav_max, n = phot_init)
#
cav = tensor([cav0 for i in range(cav_modes)])
exc = tensor([p0 for i in range(cav_modes)]) 
ryd = tensor([r0 for i in range(cav_modes)])
# |phot_n, p_n, r_ n > for n <= cav_modes
psi0 = tensor([cav, exc, ryd])

# collective bosonic a, p, r state destruction_ops
# list representation of operators to be tensored later
Id_cav = [qeye(cav_max) for i in range(cav_modes)]
Id_p = Id_cav
Id_r = Id_cav


a = ['' for i in range(cav_modes)] # to be populated by a0, a1, a2
p = ['' for i in range(cav_modes)] # to be populated by p0, p1, p2
r = ['' for i in range(cav_modes)] # to be populated by r0, r1, r2

# populate a
for i in range(cav_modes):
    temp = [qeye(cav_max) for j in range(cav_modes)]
    temp[i] = destroy(cav_max)
    temp = temp + Id_p + Id_r
    a[i] = tensor(temp)
    
# populate p
temp = []
for i in range(cav_modes):
    temp = [qeye(cav_max) for j in range(cav_modes)]
    temp[i] = destroy(cav_max)
    temp = Id_cav + temp + Id_r
    p[i] = tensor(temp)

# populate r  
temp = []
for i in range(cav_modes):
    temp = [qeye(cav_max) for j in range(cav_modes)]
    temp[i] = destroy(cav_max)
    r[i] = tensor(Id_cav + Id_p + temp)   

# cav-p exchange, p-r exchange, drive
pexchange = [(p[i].dag()*a[i] + a[i].dag()*p[i]) for i in range(cav_modes)]
rexchange = [(r[i].dag()*p[i] + p[i].dag()*r[i]) for i in range(cav_modes)]


# In[4]:

# parameters placeholders
c = [1, 1, 1]
delta_c = [1, 1, 1]
delta_e = 1
delta_2 = 1
G = [1, 1, 1] 
#U = np.zeros((cav_modes, cav_modes, cav_modes, cav_modes)) 
omega = 1
prb = 0.3 

# decay params
kappa = 1.4
gamma = 6
gamma_r = 0.1

# collapse operators
# NOTE: each mode gets its own collapse operator 
# but the decay param is identical for all modes
c_phot = [np.sqrt(kappa)*a_k for a_k in a]
c_p = [np.sqrt(gamma)*p_k for p_k in p]
c_r =  [np.sqrt(gamma_r)*r_k for r_k in r]

collapse = c_phot + c_p + c_r


# In[5]:

def g2(rho, H, dest, times, c=collapse, options=Options(tidy=False)):
    rho_m1 = dest * rho * dest.dag()
    norm = np.trace(rho_m1)
    rho_m1 = dest * (rho/norm) * dest.dag()
    rho_m1_t = mesolve(H, rho_m1, times, c, [dest.dag()*dest], options=options)
    
    num = rho_m1_t.expect[0] * norm
    den = np.power(np.trace(dest.dag()*dest*rho), 2)
    
    return num/den

# calculate H_eff given params
def H_eff(c, delta_c, delta_e, delta_2, G, omega, prb):
    # Effective Hamiltonian
    
    tot_drive = sum([c_k*a_k for c_k, a_k in zip(c, a)])
    tot_drive = prb*(tot_drive + tot_drive.dag())
    
    tot_phot = sum([coeff*op.dag()*op for coeff, op in zip(delta_c,a)])
    tot_p = delta_e * sum([op.dag()*op for op in p])
    tot_r = delta_2 * sum([op.dag()*op for op in r])
    
    tot_pexchange = sum([coeff*op for coeff, op in zip(G,pexchange)])
    tot_rexchange = omega * sum(rexchange)

    H_eff = tot_drive + tot_phot + tot_p + tot_r + tot_pexchange + tot_rexchange
    
    return H_eff

# How should I encode interaction strength?
# U = np.array(N_cav ^ 4) , (U_(a, b)(c, d)) where 
# N_cav - 1 >= a >= b >= 0 and N_cav - 1 >= c >= d >= 0 
# We start counting from 0 to match rydberg ops 
# U = np.zeros(cav_modes, cav_modes, cav_modes, cav_modes)
# 
# U(indx) for indx in indxs:

def H_ryd_int(U):
    # build indices
    temp = []
    for a in range(cav_modes):
        for b in range(a+1):
            temp.append((a,b))
            
    indxs = []
    for x in temp:
        for y in temp:
            indxs.append(tuple(list(x + y)))
            
    #print(indxs)
            
    H_ryd_int = sum([U[indx]*r[indx[0]].dag()*r[indx[1]].dag()*r[indx[2]]*r[indx[3]] for indx in indxs])
    
    return H_ryd_int

# del_omega == np.ones((Ncav, Ncav, Ncav, Ncav))
#
def H_exch_int(del_omega):
    # build indices 
    first = []
    for a in range(cav_modes):
        for b in range(a+1):
            first.append((a,b))
            
    second = []
    for c in range(cav_modes):
        for d in range(cav_modes):
            second.append((c,d))
    
    indxs = []
    for x in first:
        for y in second:
            indxs.append(tuple(list(x + y)))
            
    #print(indxs)
    
    H_exch_int = sum([del_omega[indx]*r[indx[0]].dag()*r[indx[1]].dag()*r[indx[2]]*p[indx[3]] for indx in indxs])
    H_exch_int += H_exch_int.dag()
    
    return H_exch_int 


# In[6]:


# test rydberg interaction term
# given cav_modes = 2

#U = np.ones((cav_modes, cav_modes, cav_modes, cav_modes)) 
#test = H_ryd_int(U)

#print(test)

####################

# test rydberg-p state interaction term 
# given cav_modes = 2

#del_omega = np.ones((cav_modes, cav_modes, cav_modes, cav_modes)) 
#test = H_exch_int(del_omega)

#test

#np.count_nonzero(test)

# This is missing 3 indices 


# In[ ]:


# Test [1] VRS central peak

# Testing for U_0000 (a[0]) mode

# params
c = [1, 0, 0]
delta = 0 # @ dark mode peak
delta_c = [delta, 0, 0] 
delta_e = delta
delta_2 = delta
G = [5, 0, 0] 
omega = 2
#U = np.zeros((cav_modes, cav_modes, cav_modes, cav_modes)) 
#U[(0,0,0,0)] = 0.0001 # U_0000 = 10 and everything else is 0
prb = 0.03

# Hamiltonian
H = H_eff(c, delta_c, delta_e, delta_2, G, omega, prb) #+ H_ryd_int(U)

print('Here!')

start = time.time()
# Somebody here is calling free() on a null pointer 
# which is throwing a segfault
rho_ss = steadystate(H, collapse, method='direct') # very slow, needs fixing
end = time.time()
print('SteadyState Computation Time: ', end - start)
times = np.linspace(0, 1, 3)

print('Here!')

# get g2(times) for a[0] 
g2_t = g2(rho_ss, H, a[0], times)

print('g2(0): ', g2_t[0])

# plot g2_t
fig=plt.figure()
plt.plot(times, g2_t)
plt.xlabel('t: Time')
plt.ylabel('$g^{(2)}(t)$')
plt.title('DARK' + '|mode='+ " ".join(str(x) for x in c) + '|U='+str(U[(0,0,0)]) + '| g='+str(G[0]) +  ' | $\Omega$ =' + str(omega) + ' | $\kappa$=' + str(kappa)  + ' | $\Gamma_R$=' + str(gamma_r) +' | $\Gamma_P$=' + str(gamma) + ' | PrbPwr=' + str(prb))


# In[ ]:

