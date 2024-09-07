"""
functions for the cold matter halo denity profile
"""

import numpy as np
import scipy 
from scipy import optimize, interpolate
from cosmology.basic_cosmology import func_rho_comp_0
from cosmology.overdensities import func_D_z_norm, func_r_vir, func_Delta_vir
from cosmology.variance import func_nu

# Global cache for the interpolator
_concentration_interpolator_cache = None
_zguess_interpolator_cache = None

#find z_formation given by Mead 2020 eq. 21
def func_z_formation(M, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma, f = 0.01):
    """
    k_sigma is in units of h/Mpc, PS_sigma in (Mpc/h)^3 and M in solar_mass/h
    returns the fromation redshift of a halo which is defined as
    in HMCode2020: https://arxiv.org/abs/2009.01858 in eq. 21
    """
    z = cosmo_dic['z']
    D_norm =  func_D_z_norm(z, cosmo_dic)
    def func_find_root(x, Mass):
        nu = func_nu(f*Mass, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma)
        return func_D_z_norm(x, cosmo_dic) - D_norm * nu
    if isinstance(M, (int, float)) == True:
        # Test if we find a root, if not by definition formation redshift is set to given z.
        if func_find_root(z, M)*func_find_root(100., M) > 0.:
            return z
        else:
            z_f = optimize.brentq(func_find_root, z, 100., args = (M))
            return z_f
    else:
        # Test also if we can find a root
        return np.array([optimize.brentq(func_find_root, z, 100., args=(m)) if func_find_root(z, m)*func_find_root(100., m) < 0 else z for m in M])

def func_conc_param(M, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma, recalc_c = False):
    """
    k_sigma is in units of h/Mpc, PS_sigma in (Mpc/h)^3 and M in solar_mass/h
    NOTE: Omega_0 must match with chosen PS_sigma
    returns the concentration parameter as defined in
    https://arxiv.org/abs/2009.01858 in eq. 20
    """
    global _concentration_interpolator_cache # CDM interpolator
    B = 5.196   
    concentrations = []    
    
    if _concentration_interpolator_cache is None or recalc_c:
        # First time the function is called, create the interpolator
        mass_range = np.logspace(7, 18, num=200)
        for M in mass_range:
            zf = func_z_formation(M, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma) 
            c = B*(1+zf)/(1.+cosmo_dic['z']) # Halo concentration; equation (20)
            concentrations.append(c)
        # Create an interpolation function
        _concentration_interpolator_cache = interpolate.interp1d(mass_range, concentrations, kind='linear', fill_value="extrapolate")
    
    # Use the cached interpolator
    c_eval = _concentration_interpolator_cache(M)
    return c_eval # Defined by r_vir/r_2

#function for the normaliation factor in NFW profile
def func_for_norm_factor(x):
    """
    normalisation function for the NFW profile
    """
    return (- x/(1+x)) + np.log(1+x)

#density profile in k space (fourietrafo)
def func_dens_profile_kspace(M, k, k_sigma, PS_sigma, cosmo_dic, hmcode_dic, Omega_0, Omega_0_sigma, eta_given = False, recalc_c = False):
    """
    k, k_sigma units of h/Mpc, M in solar_mass/h and PS, PS_sigma in (Mpc/h)^3 
    NOTE: be carefull, we have two k's: k is the k, where the function is evaluated 
    and k_sigma is needed for for sigma(M, z) (the same is true for the PSs)
    NOTE: Omega_0 must match with chosen PS_sigma
    returns Fourier trafo of NFW profile (dimensionless) at k as given in my masterthesis eq. 4.20
    """
    #eta is a halo shape parameter introduced my Mead in https://arxiv.org/abs/2009.01858 in Tab2
    if eta_given == True:
        eta = hmcode_dic['eta']
        nu = func_nu(M, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma)
    else:
        eta = np.array([0.]) 
        nu = 1.
    
    R_vir = np.atleast_1d(func_r_vir(M, cosmo_dic, Omega_0))
    concentration = np.atleast_1d(func_conc_param(M, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma, recalc_c = recalc_c))
    
    k_scaled = k * nu[:, None]**eta[:, None]
    k_R_vir = R_vir[:, None] * k_scaled
    a = (R_vir[:, None] / concentration[:, None]) * k_scaled
    
    def sin_integral(x):
        return scipy.special.sici(x)[0]
    def cos_integral(x):
        return scipy.special.sici(x)[1]
    
    summand1 = np.cos(a) * (cos_integral(a+k_R_vir) - cos_integral(a))
    summand2 = np.sin(a) * (sin_integral(a+k_R_vir) - sin_integral(a))
    summand3 = - np.sin(k_R_vir) / (a+k_R_vir)
    
    dens_profile_kspace = 1. / func_for_norm_factor(concentration)[:, None] * (summand1 + summand2 + summand3)
    return dens_profile_kspace 

#delta_char for the NFW profile
def func_delta_char(M, k_sigma, PS_sigma, cosmo_dic, hmcode_dic, Omega_0, Omega_0_sigma, eta_given = False): 
    """
    k_sigma units of h/Mpc, M in solar_mass/h and PS_sigma in (Mpc/h)^3 
    returns NFW profile in h^2 * M_sun/Mpc^3 at k in h/Mpc
    as given in my masterthesis eq. 4.15 (left)
    """
    #eta is a halo shape parameter introduced my Mead in https://arxiv.org/abs/2009.01858 in Tab2
    if eta_given == True:
        eta = hmcode_dic['eta']
        nu = func_nu(M, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma)
    else:
        eta = np.array([0.]) 
        nu = 1.   
    Delta_vir = func_Delta_vir(cosmo_dic, Omega_0)
    
    concentration = func_conc_param(M, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma) / (nu**eta)
    delta_char = func_rho_comp_0(Omega_0) * Delta_vir * concentration **3 / (3. * func_for_norm_factor(concentration))
    return delta_char


#density profile in real space
def NFW_profile(M, r, k_sigma, PS_sigma, cosmo_dic, hmcode_dic, Omega_0, Omega_0_sigma, eta_given = False):
    """
    r in units of Mpc/h, k_sigma in h/Mpc, M in solar_mass/h PS_sigma in (Mpc/h)^3
    returns NFW denisty profile in units of solar_mass/Mpc^3*h^2 at radius r
    """
    #eta is a halo shape parameter introduced my Mead in https://arxiv.org/abs/2009.01858 in Tab2
    if eta_given == True:
        eta = hmcode_dic['eta']
        nu = func_nu(M, k_sigma, PS_sigma, cosmo_dic, Omega_0)
    else:
        eta = np.array([0.]) 
        nu = 1.
        
    concentration = func_conc_param(M, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma) / (nu**eta)
    normalisation = func_delta_char(M, k_sigma, PS_sigma, cosmo_dic, hmcode_dic, Omega_0, Omega_0_sigma, eta_given = eta_given)
    r_s = func_r_vir(M, cosmo_dic, Omega_0) / concentration
    
    NFW_func = 1 /((r/r_s) * (1+r/r_s)**2)
    
    return normalisation * NFW_func