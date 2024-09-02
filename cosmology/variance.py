"""
functions for the varaince of the power spectrum with difference parameters
"""

#packages
import numpy as np
from scipy import integrate, optimize

#own pachages
from .basic_cosmology import *
from .overdensities import *

from solvers.root_finders import *
from numba import jit, njit

@njit
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

    
@njit
def func_sigma_r(R, k, PS): 
    """
    R in Mpc/h, k is in units of h/Mpc and PS in (Mpc/h)^3.
    returns variance of the power spectrum filterd by a spherical top hat function
    """    
    #integrand = PS * spherical_tophat_window_function(R, k) ** 2 * k ** 2
    #sigma_squared = integrate.simps(y=integrand, x=k, axis=-1) / (2. * np.pi ** 2)
    #sigma_squared = np.trapz(integrand, x=k, axis=-1) / (2. * np.pi ** 2) 
    if isinstance(R, (int, float)) == True:
        integrand = PS * spherical_tophat_window_function(R, k)[0] ** 2 * k ** 2
        sigma_squared = trapz(integrand, k) / (2. * np.pi ** 2)
    else:
        sigma_squared = np.zeros(len(R))
    
        for i in range(len(R)):
            integrand = PS * spherical_tophat_window_function(R[i], k)[0] ** 2 * k ** 2
            sigma_squared[i] = trapz(integrand, k) / (2. * np.pi ** 2)

    return np.sqrt(sigma_squared)

@jit(forceobj=True)
def func_sigma_squared_damping_twohalo(k, PS_cold): 
    """
    k is in units of h/Mpc and PS in (Mpc/h)^3.
    returns linear variance as in https://arxiv.org/pdf/2009.01858.pdf eq 10 
    in the limit R -> 0
    this is needed for the 
    """   
    integrand = PS_cold
    #sigma_squared = integrate.simps(y=integrand, x=k, axis=-1) / (2. * np.pi ** 2)
    sigma_squared = np.trapz(integrand, x=k_sigma, axis=-1) / (2. * np.pi ** 2) 
    
    return 1/3 * sigma_squared
  
@njit
def func_sigma_M(M, k, PS, Omega_0):
    """
    k is in units of h/Mpc, PS in (Mpc/h)^3 and M in solar_mass/h
    return variance of the chosen matter type, 
    NOTE: Omega_0 must match with chosen PS
    """    
    R = func_R_M(M, Omega_0)
    return func_sigma_r(R, k, PS)

@njit
def func_nu(M, k, PS_cold, Omega_0_sigma):
    """
    k is in units of h/Mpc, PS_cold in (Mpc/h)^3 and M in solar_mass/h,
    NOTE: Omega_0 must match with chosen PS_cold
    """
    delta_c = func_delta_c()
    return delta_c/func_sigma_M(M, k, PS_cold, Omega_0_sigma)

@njit
def func_R_nonlin(k, PS):
    """
    k is in units of h/Mpc and PS in (Mpc/h)^3
    returns the nonlinear scale in Mpc/h defined by
    1 = sigma(R_nonlin)
    """
    #def find_root(R):
    #    return func_sigma_r(R, k, PS) - 1.
    #R_nonlin = optimize.brentq(find_root, 1e-10, 10.)
    R_nonlin = newton(func_sigma_r, 1., y=1., args=[k, PS])

    return R_nonlin
