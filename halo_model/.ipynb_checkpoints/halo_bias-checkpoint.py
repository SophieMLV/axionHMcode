"""
functions for the halo bias
"""

#packages needed
import numpy as np
#own pachages
from cosmology.overdensities import *
from cosmology.variance import *

#halo bias as in my notes
def func_halo_bias(M, k_sigma, PS_sigma, cosmo_dic, Omega_0, p = 0.3, q = 0.707):
    """
    k_sigma units of h/Mpc, M in solar_mass/h and PK_sigma in (Mpc/h)^3
    returns halo bias (dimensionless) as in my masterthesis eq. 4.32
    """
    delta_c = func_delta_c(cosmo_dic)
    nu = func_nu(M, k_sigma, PS_sigma, cosmo_dic, Omega_0)
    halo_bias = 1. + 1./delta_c * (q*nu**2 -1 +(2*p)/(1+(q*nu**2)**p ) )
    return halo_bias

