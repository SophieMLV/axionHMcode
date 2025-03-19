"""Crate a dictionary wuith all parameters in the HMCode2020 https://arxiv.org/abs/2009.01858 as in table 2"""

import numpy as np
from scipy import optimize
from scipy import interpolate, misc
from cosmology.variance import func_sigma_r
from cosmology.overdensities import func_delta_c
from halo_model.halo_mass_function import *

def func_k_d(k, PS_cold):
    """
    k is in units of h/Mpc and PS_cold in (Mpc/h)^3
    retruns small scale dumping wavenumber of the two halo term in h/Mpc
    as in HMCode2020: https://arxiv.org/abs/2009.01858 in table 2
    """
    sigma8 = func_sigma_r(8.0, k, PS_cold)
    return 0.05699 * sigma8**(-1.089)


def func_f_dewiggle(k, PS_cold):
    """
    k is in units of h/Mpc and PS_cold in (Mpc/h)^3
    retruns two halo damping parameter
    as in HMCode2020: https://arxiv.org/abs/2009.01858 in table 2
    """
    sigma8 = func_sigma_r(8.0, k, PS_cold)

    return 0.2696 * sigma8**(0.9403)


def func_k_star(k, PS_cold):
    """
    k is in units of h/Mpc and PS_cold in (Mpc/h)^3
    retruns large scale dumping wavedumber of the one halo term in h/Mpc
    as in HMCode2020: https://arxiv.org/abs/2009.01858 in table 2
    """
    sigma8 = func_sigma_r(8.0, k, PS_cold)
    return 0.056188 * sigma8**(-1.013)


def func_eta(k, PS_cold):
    """
    k is in units of h/Mpc and PS_cold in (Mpc/h)^3
    retruns halo bloating term 
    as in HMCode2020: https://arxiv.org/abs/2009.01858 in table 2
    """
    sigma8 = func_sigma_r(8.0, k, PS_cold)
    return 0.1281 * sigma8**(-0.3644) # for Bullock2001

def func_R_nonlin_2(cosmo_dic, k, PS):
    """
    k is in units of h/Mpc and PS in (Mpc/h)^3
    returns the nonlinear scale in Mpc/h defined by
    1 = sigma(R_nonlin)
    """
    delta_c = func_delta_c(cosmo_dic['z'], cosmo_dic['Omega_ax_0'], cosmo_dic['Omega_m_0'], cosmo_dic['Omega_w_0'], cosmo_dic['G_a'])
    def find_root(lnR):
        # return func_sigma_r(np.exp(lnR), k, PS) - delta_c
        return np.log(func_sigma_r(np.exp(lnR), k, PS) / delta_c)
    
    if find_root(-10) * find_root(1.) >= 0.:
        return -1 # if equation as no soln, return -1
    else:
        R_nonlin = optimize.brentq(find_root, -10, 1.)
        # print(func_sigma_r(np.exp(R_nonlin), k, PS) - delta_c)
        return np.exp(R_nonlin)

def func_alpha_param(cosmo_dic, k, PS_cold):
    """
    k is in units of h/Mpc and PS_cold in (Mpc/h)^3
    retruns smoothing parameter
    as in HMCode2020: https://arxiv.org/abs/2009.01858 in table 2
    """
    R = np.linspace(1e-6, 1e2, 5000)
    R_nonlin = func_R_nonlin_2(cosmo_dic, k, PS_cold)
    if R_nonlin <= 1e-5: # not exactly R[0] since we would want to ensure misc.derivative does not fail even when we are close to the edge of R
        neff = -3
    else:
        # ln_R = np.log(R)
            # ln_sigma_squared = np.log(func_sigma_r(R, k, PS_cold)**2)
            # func_lnsigma_lnR = interpolate.interp1d(ln_R, ln_sigma_squared, kind = 'cubic', bounds_error=False, fill_value=0.)
            # lnsigma_lnR = misc.derivative(func_lnsigma_lnR, np.log(R_nonlin))
    
        M_nonlin = (4 * np.pi * R_nonlin**3 * func_rho_comp_0(cosmo_dic['Omega_db_0']))/ 3
        dlnsigma2_dlnR = 3* func_dlnsigma2_dlnM(M_nonlin, k, PS_cold, cosmo_dic, cosmo_dic['Omega_db_0'])

        # neff = -3 - lnsigma_lnR
        neff = -3 - dlnsigma2_dlnR
        # print(neff, R_nonlin)
    
    f_ax = cosmo_dic['Omega_ax_0']/cosmo_dic['Omega_m_0']
    m_ax = cosmo_dic['m_ax']
    z = cosmo_dic['z']
    # CDM-CDM
    a1 = 0.124 # must be m_ax-dependent
    b1 = 0.0450
    c1 = 2.260e-01
    d1 = 1.13e+00
    s1 = 1 + a1*(10**(-24)/m_ax)**b1*f_ax**c1*(1+z)**d1
    
    # Cross
    a2 = 0.0487 # must be m_ax-dependent
    b2 = 0.0450
    c2 = 2.24e-01
    d2 = 2.21e+00
    s2 = 1 + a2*(10**(-24)/m_ax)**b2*f_ax**c2*(1+z)**d2
    
    return ((1.875 * (1.603)**neff) / s1, (1.875 * (1.603)**neff) / s2)

def HMCode_param_dic(cosmo_dic, k, PS_cold):
    """
    generate dictonary with parameetrs 
    as in HMCode2020: https://arxiv.org/abs/2009.01858 in table 2
    """
    param_model_dic = {}
    param_model_dic['k_d'] = func_k_d(k, PS_cold)
    param_model_dic['f'] = func_f_dewiggle(k, PS_cold)
    param_model_dic['n_d'] = 2.853
    param_model_dic['k_star'] = func_k_star(k, PS_cold)
    param_model_dic['eta'] = func_eta(k, PS_cold)
    param_model_dic['alpha'] = func_alpha_param(cosmo_dic, k, PS_cold)
    param_model_dic['c_min'] = 5.2
    
    return param_model_dic
