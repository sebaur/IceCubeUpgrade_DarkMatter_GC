import numpy as np
from scipy import integrate 

convf = 3.08567758*10**21 # from kpc to cm

def NFW(r):
    rho_0 = 1.40 
    r_s = 16.1
    return rho_0 / ( r/r_s * (1+r/r_s)**2 ) * 1e7 * 1e-9 * 37.96
NFW = np.vectorize(NFW)

def Burkert(r):
    rho_0 = 4.13 
    r_s = 9.26
    return rho_0 / ( (1+r/r_s) * (1+(r/r_s)**2 ) ) * 1e7 * 1e-9 * 37.96
Burkert = np.vectorize(Burkert)

def Jfactor(profile, psi_in, exp):
    r_sol = 8.5
    s = np.linspace(0,50,10000)
    r = np.sqrt(r_sol**2+s**2-2*r_sol*s*np.cos(psi_in))
    if profile == 'NFW':
        rho = NFW(r)
    elif profile == 'Burkert':
        rho = Burkert(r)
    return integrate.trapz(np.power(rho,exp),s) * convf
Jfactor = np.vectorize(Jfactor)
