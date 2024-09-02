"""
functions for the critical overdensity in spherical/ellipsoidal collaps and halo overdensity
"""

#packages
import numpy as np
from scipy import integrate

#own packages
from .basic_cosmology import *

from numba import jit, njit
from solvers.integrators import *

@njit
def func_D_z_unnorm(z, Omega_m_0, Omega_w_0):
    """
    returns unnormalised grwoth function 
    """
    #def integrand(x):
    #    return (1+x)/func_E_z(x, Omega_m_0, Omega_w_0)**3
    
    z_array = np.linspace(z, 100, 2000)
    integrand = (1+z_array) / func_E_z(z_array, Omega_m_0, Omega_w_0)**3
    
    #factor = 5*cosmo_dic['Omega_m_0']/2
    factor = 5 * Omega_m_0 / 2
    #D = factor * func_E_z(z, Omega_m_0, Omega_w_0) * integrate.quad(integrand, z, np.inf)[0]
    #D = factor * func_E_z(z, Omega_m_0, Omega_w_0) * np.trapz(integrand, x=z_array) # tested to better than 0.5%
    D = factor * func_E_z(z, Omega_m_0, Omega_w_0) * trapz(integrand, z_array) # now with Numba trapz
    return D


@njit
def func_D_z_norm(z,  Omega_m_0, Omega_w_0):
    """
    returns normlised grwoth function, ie D(0)= 1
    this is used to scale the power spectrum for diffrent z's
    """
    normalisation = func_D_z_unnorm(0.,  Omega_m_0, Omega_w_0)
    growth = func_D_z_unnorm(z,  Omega_m_0, Omega_w_0)
    
    return growth/normalisation

@njit
def func_delta_c():
    """
    returns critical denity for spherical/ellepsoidal collapse for LCDM cosmos
    """                    
    return 1.686

@njit    
def func_Delta_vir(z, Omega_0, Omega_m_0, Omega_w_0):
    """
    halo overdensity for LCDM comos,
    make the change, that only matter of the type 
    Omega_0 is take into accound for the overdensity
    """
    #x = func_Omega_comp_z(cosmo_dic, Omega_0) -1
    #x = func_Omega_comp_z(cosmo_dic, cosmo_dic['Omega_m_0']) -1
    x = func_Omega_comp_z(z, Omega_m_0, Omega_m_0, Omega_w_0) - 1 # following the commented lines
    return (18*np.pi**2 + 82*x - 39*x**2) / (x+1)


@njit
def func_r_vir(z, M, Omega_0, Omega_m_0, Omega_w_0):
    """
    M in solar_mass/h where M is matter of the type Omega_0
    returns comoving virial radius in Mpc/h
    """
    Delta_vir = func_Delta_vir(z, Omega_0, Omega_m_0, Omega_w_0)
    return (3. * M / (4. * np.pi * func_rho_comp_0(Omega_0) * Delta_vir ))**(1./3.)

