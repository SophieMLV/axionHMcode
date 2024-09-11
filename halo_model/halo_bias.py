"""
functions for the halo bias
"""

from cosmology.overdensities import func_delta_c
from cosmology.variance import func_nu

from numba import njit

#halo bias as in my notes
@njit
def func_halo_bias(M, k, PS, Omega_0, Omega_m_0, Omega_w_0, z, G_a, p = 0.3, q = 0.707):
    """
    k units of h/Mpc, M in solar_mass/h and PK in (Mpc/h)^3
    returns halo bias (dimensionless) as in eq 24 in https://arxiv.org/abs/2209.13445
    """
    delta_c = func_delta_c(z, Omega_m_0, Omega_w_0, G_a)
    nu = func_nu(M, k, PS, Omega_0, Omega_m_0, Omega_w_0, z, G_a)
    halo_bias = 1. + 1./delta_c * (q*nu**2 -1 +(2*p)/(1+(q*nu**2)**p ) )
    
    return halo_bias

