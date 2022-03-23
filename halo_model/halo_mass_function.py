"""
functions for the halo mass function
"""

#packages needed
import numpy as np
from scipy import integrate

#own pachages
from cosmology.basic_cosmology import *
from cosmology.variance import *

def func_sheth_tormen(v, A = 0.3222, p = 0.3, q = 0.707):
    """
    returns Sheth-Tormen fit function from 1999. A = 0.3222, p = 0.3 and q = 0.707 by default.
    v := delta_c/sigma(M), thus dimensionless by definition.
    """
    return A * (2. * q / np.pi) ** 0.5 * v * (1. + (q ** 0.5 * v) ** (-2. * p)) * np.exp((-q * v ** 2.) / 2.)


def func_term_derivative_sigma2_M(R, k, conditional_return=True):
    """
    returns integrand for the halo mass function as in my masterthesis eq. 4.24
    """
    x = np.outer(R, k)
    term1 = np.sin(x) * (1.0 - (3.0 / x ** 2)) + 3.0 / x * np.cos(x)
    term2 = np.sin(x) - x * np.cos(x)
    if conditional_return:
        # takes care of numerical issues for x -> 0 (did expand the term and checked)
        return np.where(x > 1e-3, term1 * term2, 0.0)
    else:
        return term1 * term2

    
def func_dlnsigma2_dlnM(M, k_sigma, PS_sigma, cosmo_dic, Omega_0, Omega_0_sigma):
    """
    k_sigma units of h/Mpc, M in solar_mass/h and PS_sigma in (Mpc/h)^3 
    returns integral for the halo mass function as in my masterthesis eq. 4.23
    """
    R = func_R_M(M, cosmo_dic, Omega_0) #translate M into R 
    sigma = func_sigma_M(M, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma)

    integrand = func_term_derivative_sigma2_M(R, k_sigma) * PS_sigma / k_sigma**2
    
    return integrate.simps(y=integrand, x = k_sigma, axis=-1) * 3./ ( np.pi**2 * sigma**2 * R**4)


def func_halo_mass_function(M, k_sigma, PS_sigma, cosmo_dic, Omega_0, Omega_0_sigma):
    """
    k_sigma in h/Mpc, M in solar_mass/h and PS_sigma in (Mpc/h)^3 
    returns halo mass functions in untits h^4/(M_sun Mpc^3)
    with teh definition as in my  masterthesis eq. 4.21
    """
    nu = func_nu(M, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma) 
    return 1./2. * func_rho_comp_0(Omega_0) / M**2 * func_sheth_tormen(nu) \
           * np.abs(func_dlnsigma2_dlnM(M, k_sigma, PS_sigma, cosmo_dic, Omega_0, Omega_0_sigma))