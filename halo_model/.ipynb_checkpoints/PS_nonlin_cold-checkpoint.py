"""
functions for neutrino density profile
"""

#packages needed
import numpy as np
from scipy import integrate

#own pachages
from cosmology.basic_cosmology import *
from cold_density_profile import *
from halo_bias import *
from halo_mass_function import *



def func_non_lin_PS_matter(M, k, PS, k_sigma, PS_sigma, cosmo_dic, hmcode_dic, Omega_0, Omega_0_sigma, 
                           alpha = False, eta_given = False, nu_one_halo=False, ax_one_halo=False, one_halo_damping = False, two_halo_damping = False):    
    """ 
    The cold halo model se master thesis eq. 4.9 with (if set to True) the modifications of HMcode2020 https://arxiv.org/abs/2009.01858
    Since we work with axions, I indroduce the possibility to tread the axions as the HMcode2020 treates the neutrinos
    by substracting them from the one halo term.
    k, k_sigma units of h/Mpc, M in solar_mass/h and PS, PS_sigma in (Mpc/h)^3 
    NOTE:be carefull, we have two k's and PS's: k is the k, where the function is evaluated 
    and k_sigma is needed for for sigma(M, z), same for PS
    returns non-lin power spectrum of matter or cold matter in (Mpc/h)^3 at k
    as well as th one halo and two halo term
    """     
    dens_profile_arr = func_dens_profile_kspace(M, k, k_sigma, PS_sigma, cosmo_dic, hmcode_dic, Omega_0, Omega_0_sigma, eta_given = eta_given)
    halo_mass_func_arr = func_halo_mass_function(M, k_sigma, PS_sigma, cosmo_dic, Omega_0, Omega_0_sigma)
    
    integrand_arr_one = M[:, None]**2 * halo_mass_func_arr[:, None] * dens_profile_arr**2 
    #no neutrinos in halos
    if ax_one_halo == True:
        f_ax = cosmo_dic['Omega_ax_0']/cosmo_dic['Omega_m_0']
        one_halo = (1-f_ax)**2 * integrate.simps(integrand_arr_one, x = M, axis = 0)/ func_rho_comp_0(Omega_0)**2  
    else:
        one_halo = integrate.simps(integrand_arr_one, x = M, axis = 0)/ func_rho_comp_0(Omega_0)**2  
    #one halo damping
    if one_halo_damping == True:
        one_halo = one_halo * (k/hmcode_dic['k_star'])**4 / (1+(k/hmcode_dic['k_star'])**4)
    else:
        one_halo = one_halo
            
    #two halo damping and some extra factors in the two halo term to take care of nummerical issues.
    # see appendix A in https://arxiv.org/abs/2005.00009
    if two_halo_damping == True:
        halo_bias_arr = func_halo_bias(M, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma)
        integrand_arr_two = M[:, None] * halo_mass_func_arr[:, None] * halo_bias_arr[:, None] * dens_profile_arr
        
        #summand to take care of nummericals issues of the integral, see appendix A in https://arxiv.org/abs/2005.00009
        summand2 = func_dens_profile_kspace(np.min(M), k, k_sigma, PS_sigma, cosmo_dic, hmcode_dic, Omega_0, Omega_0_sigma, eta_given = eta_given) \
                               * ( 1 - integrate.simps(M[:, None] * halo_mass_func_arr[:, None] * halo_bias_arr[:, None], x = M, axis = 0) / func_rho_comp_0(Omega_0)) 
        factor2 = integrate.simps(integrand_arr_two, x = M, axis = 0) / func_rho_comp_0(Omega_0) + summand2 
        
        two_halo = PS * factor2**2 * (1-hmcode_dic['f'] * (k/hmcode_dic['k_d'])**hmcode_dic['n_d']/(1+(k/hmcode_dic['k_d'])**hmcode_dic['n_d']))
    else:
        halo_bias_arr = func_halo_bias(M, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma)
        integrand_arr_two = M[:, None] * halo_mass_func_arr[:, None] * halo_bias_arr[:, None] * dens_profile_arr
        #summand2 take care of nummericals issues of the integral, see appendix A in https://arxiv.org/abs/2005.00009
        summand2 = func_dens_profile_kspace(np.min(M), k, k_sigma, PS_sigma, cosmo_dic, hmcode_dic, Omega_0, Omega_0_sigma, eta_given = eta_given) \
                               * ( 1 - integrate.simps(M[:, None] * halo_mass_func_arr[:, None] * halo_bias_arr[:, None], x = M, axis = 0) / func_rho_comp_0(Omega_0)) 
        factor2 = integrate.simps(integrand_arr_two, x = M, axis = 0) / func_rho_comp_0(Omega_0) + summand2 
        
        two_halo = PS * factor2**2
    
    #smooth the transition
    if alpha == True:
        alpha_param = hmcode_dic['alpha']   
    else:
        alpha_param = 1
        
    
    return (one_halo**(alpha_param) + two_halo[0][:]**(alpha_param))**(1/alpha_param), one_halo, two_halo[0][:]





