"""
functions for the varaince of the power spectrum with difference parameters
"""

#packages
import numpy as np
from scipy import integrate, optimize

#own pachages
from basic_cosmology import *
from overdensities import *


def spherical_tophat_window_function(R, k, conditional_return=True):
    """
    k is in units of h/Mpc and R in Mpc/h
    returns window function for variance of the power spectrum
    """    
    x = np.outer(R, k)
    if conditional_return == True:
        # takes care of numerical issues for x -> 0 (did expand the term and checked)
        return np.where(x > 1.4e-6, (3.0 / x ** 3) * (np.sin(x) - x * np.cos(x)), 1.0)
    else:
        return 3.0 / x ** 3 * (np.sin(x) - x * np.cos(x))

    

def func_sigma_r(R, k, PS): 
    """
    R in Mpc/h, k is in units of h/Mpc and PS in (Mpc/h)^3.
    returns variance of the power spectrum filterd by a spherical top hat function
    """    
    integrand = PS * spherical_tophat_window_function(R, k) ** 2 * k ** 2
    sigma_squared = integrate.simps(y=integrand, x=k, axis=-1) / (2. * np.pi ** 2)
    
    return np.sqrt(sigma_squared)


def func_sigma_squared_damping_twohalo(k_sigma, PS_sigma): 
    """
    k is in units of h/Mpc and PS in (Mpc/h)^3.
    returns linear variance as in https://arxiv.org/pdf/2009.01858.pdf eq 10 
    in the limit R -> 0
    this is needed for the 
    """   
    integrand = PS_sigma
    sigma_squared = integrate.simps(y=integrand, x=k_sigma, axis=-1) / (2. * np.pi ** 2)
    
    return 1/3 * sigma_squared
  

def func_sigma_M(M, k, PS, cosmo_dic, Omega_0):
    """
    k is in units of h/Mpc, PS in (Mpc/h)^3 and M in solar_mass/h
    return variance of the chosen matter type, 
    NOTE: Omega_0 must match with chosen PS
    """    
    R = func_R_M(M, cosmo_dic, Omega_0)
    return func_sigma_r(R, k, PS)


def func_nu(M, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma):
    """
    k_sigma is in units of h/Mpc, PS_sigma in (Mpc/h)^3 and M in solar_mass/h,
    NOTE: Omega_0 must match with chosen PS_sigma
    """
    delta_c = func_delta_c(cosmo_dic)
    return delta_c/func_sigma_M(M, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma)


def func_R_nonlin(cosmo_dic, k, PS):
    """
    k is in units of h/Mpc and PS in (Mpc/h)^3
    returns the nonlinear scale in Mpc/h defined by
    1 = sigma(R_nonlin)
    """
    def find_root(R):
        return func_sigma_r(R, k, PS) - 1.
    R_nonlin = optimize.brentq(find_root, 1e-10, 10.)
    return R_nonlin