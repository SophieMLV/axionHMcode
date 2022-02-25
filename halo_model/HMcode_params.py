"""Crate a dictionary wuith all parameters in the HMCode2020 https://arxiv.org/abs/2009.01858 as in table 2"""

#packages needed
import numpy as np
from scipy import interpolate, misc

#own pachages
from cosmology.variance import *
from cosmology.overdensities import *

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
    return 0.1281 * sigma8**(-0.3644)

def func_R_nonlin_2(cosmo_dic, k, PS):
    """
    k is in units of h/Mpc and PS in (Mpc/h)^3
    returns the nonlinear scale in Mpc/h defined by
    1 = sigma(R_nonlin)
    """
    delta_c = func_delta_c(cosmo_dic)
    def find_root(R):
        return func_sigma_r(R, k, PS) - delta_c
    R_nonlin = optimize.brentq(find_root, 1e-10, 10.)
    return R_nonlin

    
def func_alpha_param(cosmo_dic, k, PS_cold, LCDM=True):
    """
    k is in units of h/Mpc and PS_cold in (Mpc/h)^3
    retruns smoothing parameter
    as in HMCode2020: https://arxiv.org/abs/2009.01858 in table 2
    """
    R = np.linspace(1e-3, 1e2, 1000)
    R_nonlin = func_R_nonlin_2(cosmo_dic, k, PS_cold)
    ln_R = np.log(R)
    ln_sigma_squared = np.log(func_sigma_r(R, k, PS_cold)**2)
    func_lnsigma_lnR = interpolate.interp1d(ln_R, ln_sigma_squared, kind = 'cubic')
    lnsigma_lnR = misc.derivative(func_lnsigma_lnR, np.log(R_nonlin))
    neff = -3 - lnsigma_lnR
    
    return 1.875 * 1.603**neff



def HMCode_param_dic(cosmo_dic, k_sigma, PS_cold_sigma):
    """
    generate dictonary with parameetrs 
    as in HMCode2020: https://arxiv.org/abs/2009.01858 in table 2
    """
    param_model_dic = {}
    param_model_dic['k_d'] = func_k_d(k_sigma, PS_cold_sigma)
    param_model_dic['f'] = func_f_dewiggle(k_sigma, PS_cold_sigma)
    param_model_dic['n_d'] = 2.853
    param_model_dic['k_star'] = func_k_star(k_sigma, PS_cold_sigma)
    param_model_dic['eta'] = func_eta(k_sigma, PS_cold_sigma)
    param_model_dic['alpha'] = func_alpha_param(cosmo_dic, k_sigma, PS_cold_sigma)
    
    return param_model_dic