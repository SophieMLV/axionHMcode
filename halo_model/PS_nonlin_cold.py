"""
functions for neutrino density profile
"""

from scipy import integrate

#own pachages
from cosmology.basic_cosmology import *
from .cold_density_profile import func_dens_profile_kspace
from .halo_bias import func_halo_bias
from .halo_mass_function import func_halo_mass_function



def func_non_lin_PS_matter(M, k, PS, cosmo_dic, hmcode_dic, Omega_0, 
                           alpha = False, eta_given = False, ax_one_halo=False, 
                           one_halo_damping = False, two_halo_damping = False, 
                           concentration_param = False, full_2h=False, 
                           axion_dic=None):
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
    # define the concentraion param for the cold matter profile
    if concentration_param == True:
        c_min = hmcode_dic['c_min']
    else:
        c_min = 4.

    # quantities for the one halo integral
    dens_profile_arr = func_dens_profile_kspace(M, k, PS, cosmo_dic, hmcode_dic, Omega_0, c_min, eta_given = eta_given, axion_dic=axion_dic)
    halo_mass_func_arr = func_halo_mass_function(M, k, PS, cosmo_dic, Omega_0)
    f_ax = cosmo_dic['Omega_ax_0']/cosmo_dic['Omega_m_0']
    
    integrand_arr_one = M[:, None]**2 * halo_mass_func_arr[:, None] * dens_profile_arr**2 
    # No neutrinos in halos
    if ax_one_halo == True:
        one_halo = (1-f_ax)**2 * integrate.simpson(integrand_arr_one, x = M, axis = 0)/ func_rho_comp_0(Omega_0)**2  
    else:
        one_halo = integrate.simpson(integrand_arr_one, x = M, axis = 0)/ func_rho_comp_0(Omega_0)**2  
    
    # One halo damping
    if one_halo_damping == True:
        one_halo = one_halo * (k/hmcode_dic['k_star'])**4 / (1+(k/hmcode_dic['k_star'])**4)
    else:
        one_halo = one_halo
            
    #two halo damping and some extra factors in the two halo term to take care of nummerical issues.
    # see appendix A in https://arxiv.org/abs/2005.00009

    if two_halo_damping == True:

        if full_2h == True:
            halo_bias_arr = func_halo_bias(M, k, PS, cosmo_dic['Omega_ax_0'], Omega_0, cosmo_dic['Omega_m_0'], cosmo_dic['Omega_w_0'], cosmo_dic['z'], cosmo_dic['G_a'])
            integrand_arr_two = M[:, None] * halo_mass_func_arr[:, None] * halo_bias_arr[:, None] * dens_profile_arr
            
            #summand to take care of nummericals issues of the integral, see appendix A in https://arxiv.org/abs/2005.00009
            summand2 = func_dens_profile_kspace(np.min(M), k, PS, cosmo_dic, hmcode_dic, Omega_0, 
                                                c_min, eta_given = eta_given, axion_dic=axion_dic) \
                                * ( 1 - integrate.simpson(M[:, None] * halo_mass_func_arr[:, None] * halo_bias_arr[:, None], x = M, axis = 0) / func_rho_comp_0(Omega_0)) 
            factor2 = integrate.simpson(integrand_arr_two, x = M, axis = 0) / func_rho_comp_0(Omega_0) + summand2 
            
            two_halo = PS * factor2[0]**2 * (1-hmcode_dic['f'] * (k/hmcode_dic['k_d'])**hmcode_dic['n_d']/(1+(k/hmcode_dic['k_d'])**hmcode_dic['n_d']))
        else:
            two_halo = PS * (1-hmcode_dic['f'] * (k/hmcode_dic['k_d'])**hmcode_dic['n_d']/(1+(k/hmcode_dic['k_d'])**hmcode_dic['n_d']))

    else:
        if full_2h == True:
            halo_bias_arr = func_halo_bias(M, k, PS, cosmo_dic['Omega_ax_0'], Omega_0, cosmo_dic['Omega_m_0'], cosmo_dic['Omega_w_0'], cosmo_dic['z'], cosmo_dic['G_a'])
            integrand_arr_two = M[:, None] * halo_mass_func_arr[:, None] * halo_bias_arr[:, None] * dens_profile_arr
            #summand2 take care of nummericals issues of the integral, see appendix A in https://arxiv.org/abs/2005.00009
            summand2 = func_dens_profile_kspace(np.min(M), k, PS, cosmo_dic, hmcode_dic, Omega_0, 
                                                c_min, eta_given = eta_given, axion_dic=axion_dic) \
                                * ( 1 - integrate.simpson(M[:, None] * halo_mass_func_arr[:, None] * halo_bias_arr[:, None], x = M, axis = 0) / func_rho_comp_0(Omega_0)) 
            factor2 = integrate.simpson(integrand_arr_two, x = M, axis = 0) / func_rho_comp_0(Omega_0) + summand2 
            
            two_halo = PS * factor2[0]**2
        else:
            two_halo = PS
    #smooth the transition
    if alpha == True:
        k_piv = 0.5*10**0 # h/cMpc
        Deltak = 0.1 # h/cMpc
        z = cosmo_dic['z']
        alphacdm = hmcode_dic['alpha'][0]
        alpha1cdm = max(alphacdm, 1.1/(1+z)) # make sure alpha is not too small on large scales
        if f_ax < 0.01:
            alpha1cdm = alphacdm
        alphacdm = alpha1cdm + (alphacdm - alpha1cdm)/(1 + np.exp((k_piv - k)/Deltak)) 
        # Logistic function, which provides a smooth transition between the desired αα values. Δk controls the width of the transition region around kpivot
    else:
        alphacdm = 1
        
    
    return (one_halo**(alphacdm) + two_halo**(alphacdm))**(1/alphacdm), one_halo, two_halo, halo_mass_func_arr, M, dens_profile_arr