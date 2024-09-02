"""
some basic function needed for the non-linear power spectrum
"""

#import packages
import numpy as np
#from astropy import constants as const

from numba import jit, njit

@njit
def func_E_z(z, Omega_m_0, Omega_w_0):
    """
    return time dependence of the Hubble constants
    """
    #E2 = cosmo_dic['Omega_m_0']*(1+z)**3 + cosmo_dic['Omega_w_0']
    E2 = Omega_m_0*(1+z)**3 + Omega_w_0
    return np.sqrt(np.abs(E2))

@njit
def func_H_z(z, h, Omega_m_0, Omega_w_0):
    """
    H(z) is given in units of H0 eg km/(s * Mpc)
    """
    #z = cosmo_dic['z']
    #return 100*cosmo_dic['h'] * func_E_z(z, cosmo_dic)
    return 100 * h * func_E_z(z, Omega_m_0, Omega_w_0)

@njit
def func_Omega_comp_z(z, Omega_0, Omega_m_0, Omega_w_0):
    """
    returns matter componenet specified by Omega_0
    density parameter as a function of redshift
    """
    #z = cosmo_dic['z']
    return Omega_0 * (1+z)**3 / func_E_z(z, Omega_m_0, Omega_w_0)**2

@njit
def func_rho_comp_0(Omega_0):
    """
    present density of component Omega_0
    retruns demsity in units of solar_mass/Mpc^3*h^2
    """
    grav_const = 4.300917270069977e-09 #const.G.to('km**2*Mpc/(Msun*s**2)').value
    rho_crit_0 = 3. * 100**2 / (8. * np.pi * grav_const)
    
    return Omega_0  * rho_crit_0

@njit
def func_R_M(M, Omega_0):
    """
    M in solar_mass/h, where M is matter of the type Omega_0
    returns R in Mpc/h
    """
    R = (3 * M / (4 * np.pi * func_rho_comp_0(Omega_0)) )**(1./3.)
    return R  

