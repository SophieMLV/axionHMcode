"""
functions for the varaince of the power spectrum with difference parameters
"""

#packages
import numpy as np
from scipy import integrate, optimize

#own pachages
#at the moment in the same folder, change path when move the files
from basic_functions_cosmology import *
from overdensities import *

# variance power spectrum
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

    
#compute variance via simpson from Jurek
def func_sigma_r(R, k, PS): 
    """
    R in Mpc/h, k is in units of h/Mpc and PS in (Mpc/h)^3.
    returns variance of the power spectrum filterd by a spherical top hat function
    """    
    integrand = PS * spherical_tophat_window_function(R, k) ** 2 * k ** 2
    #print(np.max(spherical_tophat_window_function(R, k)))
    sigma_squared = integrate.simps(y=integrand, x=k, axis=-1) / (2. * np.pi ** 2)
    #dlnk = np.log(k[1]/k[0])
    #sigma_squared = integrate.trapz(y=integrand*k, dx=dlnk, axis=-1) / (2. * np.pi ** 2)
    
    return np.sqrt(sigma_squared)


#compute linear variance need in Mead2020
def func_sigma_squared_damping_twohalo(k_sigma, PS_sigma): 
    """
    k is in units of h/Mpc and PS in (Mpc/h)^3.
    returns variance of the linear filed filtered by a spherical top hat function in the limet R -> 0
    """   
    integrand = PS_sigma
    sigma_squared = integrate.simps(y=integrand, x=k_sigma, axis=-1) / (2. * np.pi ** 2)
    
    return 1/3 * sigma_squared


def func_R_M(M, cosmo_dic, Omega_0):
    """
    M in solar_mass/h, where M is matter of the type Omega_0
    returns R in Mpc/h
    """
    R = (3 * M / (4 * np.pi * func_rho_comp_0(Omega_0)) )**(1./3.)
    return R    


def func_sigma_M(M, k, PS, cosmo_dic, Omega_0):
    """
    k is in units of h/Mpc, PS in (Mpc/h)^3 and M in solar_mass/h
    return variance of the chosen matter type, 
    be carefull PS musst match the chosen matter type and redshift
    """    
    R = func_R_M(M, cosmo_dic, Omega_0)
    return func_sigma_r(R, k, PS)


def func_nu(M, k_sigma, PS_sigma, cosmo_dic, Omega_0, LCDM = True):
    """
    k_sigma is in units of h/Mpc, PS_sigma in (Mpc/h)^3 and M in solar_mass/h
    you can adjust which kind of critical density you want with the parameter LCDM
    returns variance of the chosen matter type, 
    be carefull PS musst match the chosen matter type and redshift
    """
    
    if LCDM == True:
        delta_c = func_delta_c_LCDM(cosmo_dic)
    else:
        delta_c = func_delta_c(cosmo_dic)
        
    return delta_c/func_sigma_M(M, k_sigma, PS_sigma, cosmo_dic, Omega_0)


def func_r_vir(M, cosmo_dic, Omega_0, LCDM = True):
    """
    M in solar_mass/h where M is matter of the type Omega_0
    you can adjust which kind of halo overdnsity you want with the parameter LCDM
    returns virial radius in Mpc/h
    """
    if LCDM == True:
        Delta_vir = func_Delta_vir_LCDM(cosmo_dic, Omega_0)
    else:
        Delta_vir = func_Delta_vir(cosmo_dic)
        
    return (3. * M / (4. * np.pi * func_rho_comp_0(Omega_0) * Delta_vir ))**(1./3.) # comoving virial radius

def func_R_nonlin(cosmo_dic, k, PS_cold, LCDM=True):
    """
    k is in units of h/Mpc and PS_cold in (Mpc/h)^3
    function needed for the smoothing parameter alpha, see function below
    """
    if LCDM == True:
        delta_c = func_delta_c_LCDM(cosmo_dic)#, cosmo_dic['Omega_db_0'])
    else:
        delta_c = func_delta_c(cosmo_dic)
    def find_root(R):
        return func_sigma_r(R, k, PS_cold) - 1.#1.686 #delta_c
    #print('UPPER', find_root(10))
    #print('LOWER', find_root(1e-40))
    R_nonlin = optimize.brentq(find_root, 1e-10, 10.)
    #print(R_nonlin)
    #k_nl = 0.14 * (1+cosmo_dic['z'])**(2/(2+cosmo_dic['ns'])) / cosmo_dic['h']
    #R_nonlin = 2*np.pi / k_nl
    return R_nonlin