"""
functions for the critical overdensity in spherical/ellipsoidal collaps and halo overdensity
"""

#packages
import numpy as np
from scipy import integrate

#own packages
from basic_cosmology import *


def func_D_z_unnorm(z, cosmo_dic):
    """
    returns unnormalised grwoth function 
    """
    def integrand(x):
        return (1+x)/func_E_z(x, cosmo_dic)**3
    
    factor = 5*cosmo_dic['Omega_m_0']/2
    D = factor * func_E_z(z, cosmo_dic) * integrate.quad(integrand, z, np.inf)[0]
    return D

def func_D_z_norm(z, cosmo_dic):
    """
    returns normlised grwoth function, ie D(0)= 1
    this is used to scale the power spectrum for diffrent z's
    """
    normalisation = func_D_z_unnorm(0., cosmo_dic)
    growth = func_D_z_unnorm(z, cosmo_dic)
    
    return growth/normalisation


def func_delta_c(cosmo_dic):
    """
    returns critical denity for spherical/ellepsoidal collapse for LCDM cosmos
    """                    
    return 1.686
    
def func_Delta_vir(cosmo_dic, Omega_0):
    """
    halo overdensity for LCDM comos,
    make the change, that only matter of the type 
    Omega_0 is take into accound for the overdensity
    """
    #x = func_Omega_comp_z(cosmo_dic, Omega_0) -1
    x = func_Omega_comp_z(cosmo_dic, cosmo_dic['Omega_m_0']) -1
    return (18*np.pi**2 + 82*x - 39*x**2) / (x+1)



def func_r_vir(M, cosmo_dic, Omega_0):
    """
    M in solar_mass/h where M is matter of the type Omega_0
    returns comoving virial radius in Mpc/h
    """
    Delta_vir = func_Delta_vir(cosmo_dic, Omega_0)
    return (3. * M / (4. * np.pi * func_rho_comp_0(Omega_0) * Delta_vir ))**(1./3.)

