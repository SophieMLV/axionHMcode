"""
functions for the critical overdensity in spherical/ellipsoidal collaps and halo overdensity
"""

import numpy as np
from scipy import integrate
from .basic_cosmology import func_E_z, func_Omega_comp_z, func_rho_comp_0

from numba import jit, njit
from solvers.integrators import *

@njit
def func_D_z_unnorm(z, Omega_m_0, Omega_w_0):
    """
    returns unnormalised grwoth function 
    """
    
    z_array = np.linspace(z, 100, 2000)
    integrand = (1+z_array) / func_E_z(z_array, Omega_m_0, Omega_w_0)**3
    
    factor = 5 * Omega_m_0 / 2
    D = factor * func_E_z(z, Omega_m_0, Omega_w_0) * trapz(integrand, z_array) # now with Numba trapz
    return D


@njit
def func_D_z_norm(z, Omega_m_0, Omega_w_0):
    """
    returns normlised grwoth function, ie D(0)= 1
    this is used to scale the power spectrum for diffrent z's
    """
    normalisation = func_D_z_unnorm(0.,  Omega_m_0, Omega_w_0)
    growth = func_D_z_unnorm(z,  Omega_m_0, Omega_w_0)
    
    return growth/normalisation

# @njit
def func_D_z_unnorm_int(z, Omega_m_0, Omega_w_0):
    """
    integral gwoth of the unmoralised growth function
    as in https://arxiv.org/pdf/2009.01858 eq. A5
    """

    f = lambda y, x: func_E_z(x, Omega_m_0, Omega_w_0)/(1+x)*(1+y)/func_E_z(y, Omega_m_0, Omega_w_0)**3
    G =  5 * Omega_m_0 / 2 * integrate.dblquad(f, z, 10000, lambda x:x, 10000)[0] # now with Numba trapz
    return G

@njit
def func_delta_c(z, Omega_ax_0, Omega_m_0, Omega_w_0, G_a, version):
    """
    returns critical denity for spherical/ellepsoidal collapse for LCDM cosmos
    """ 
    g_a = func_D_z_unnorm(z, Omega_m_0, Omega_w_0)*(1+z)
    p_10 = -0.0069
    p_11 = -0.0208
    p_12 = 0.0312
    p_13 = 0.0021
    p_20 = 0.0001
    p_21 = -0.0647
    p_22 = -0.0417
    p_23 = 0.0646
    f_1 = p_10 + p_11*(1-g_a) + p_12*(1-g_a)**2 + p_13*(1-G_a*(1+z))
    f_2 = p_20 + p_21*(1-g_a) + p_22*(1-g_a)**2 + p_23*(1-G_a*(1+z))

    Omega_m_z = func_Omega_comp_z(z, Omega_m_0, Omega_m_0, Omega_w_0)


    alpha_1 = 1
    alpha_2 = 0  
    f_frac = Omega_ax_0/Omega_m_0

    if version == 'dome':
        return 1.686 *(1-0.041*f_frac)* ( 1 + f_1*np.log10(Omega_m_z)**alpha_1 + f_2*np.log10(Omega_m_z)**alpha_2)
    else:
        return 1.686 * ( 1 + f_1*np.log10(Omega_m_z)**alpha_1 + f_2*np.log10(Omega_m_z)**alpha_2)


@njit    
def func_Delta_vir(z, Omega_ax_0, Omega_m_0, Omega_w_0, G_a, version):
    """
    halo overdensity for LCDM comos,
    make the change, that only matter of the type 
    Omega_0 is take into accound for the overdensity
    """
    g_a = func_D_z_unnorm(z, Omega_m_0, Omega_w_0)*(1+z)
    p_10 = -0.79
    p_11 = -10.17
    p_12 = 2.51
    p_13 = 6.51
    p_20 = -1.89
    p_21 = 0.38
    p_22 = 18.8
    p_23 = -15.87
    f_1 = p_10 + p_11*(1-g_a) + p_12*(1-g_a)**2 + p_13*(1-G_a*(1+z))
    f_2 = p_20 + p_21*(1-g_a) + p_22*(1-g_a)**2 + p_23*(1-G_a*(1+z))

    Omega_m_z = func_Omega_comp_z(z, Omega_m_0, Omega_m_0, Omega_w_0)

    f_frac = Omega_ax_0/Omega_m_0

    alpha_1 = 1
    alpha_2 = 2
    if version == 'dome':
        return 177.7 *(1+0.763*f_frac) * ( 1 + f_1*np.log10(Omega_m_z)**alpha_1 + f_2*np.log10(Omega_m_z)**alpha_2)
    else:
        return 177.7 * ( 1 + f_1*np.log10(Omega_m_z)**alpha_1 + f_2*np.log10(Omega_m_z)**alpha_2)


# @njit
def func_r_vir(z, M, Omega_ax_0, Omega_0, Omega_m_0, Omega_w_0, G_a, version):
    """
    M in solar_mass/h where M is matter of the type Omega_0
    returns comoving virial radius in Mpc/h
    """
    Delta_vir = func_Delta_vir(z, Omega_ax_0, Omega_m_0, Omega_w_0, G_a, version)
    return (3. * M / (4. * np.pi * func_rho_comp_0(Omega_0) * Delta_vir ))**(1./3.)