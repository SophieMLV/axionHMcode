"""functions for my full axion halo model and simpified halo models"""

#packages needed
import numpy as np
from scipy import integrate
#own pachages
from cosmology.basic_cosmology import *
from .halo_bias import *
from .halo_mass_function import *
from .cold_density_profile import *
from .axion_density_profile import *
from .PS_nonlin_cold import *



def func_full_halo_model_ax(M, power_spec_dic, cosmo_dic, hmcode_dic, axion_dic, 
                            alpha = False, eta_given = False, one_halo_damping = True, 
                            two_halo_damping = False, concentration_param=False, full_2h=True):
    """ 
    Full Halo Model with Axions, see BLA for the full formula
    all modifications form HMcode2020 https://arxiv.org/abs/2009.01858 can be used (if wanted)
    by default use only the damping in the one halo term on large scales (for all parts) 
    to ensure the correct behaviour on large scales
    all other HMcode2020 modifications are only apllied to the cold part 
    because they are only calibrated to CDM.
    M in solar_mass/h
    in power_spec_dic and power_spec_dic_sigma the PS and k's are stored 
    and all units are either in (Mpc/h)^3 or h/Mpc
    NOTE: Two PS dicionaries are needed bacause we have to be carefull with the calulation of sigma
    returns total non-linear matter power spectrum in (Mpc/h)^3 at k
    """ 
    # define the concentraion param for the cold matter profile
    if concentration_param == True:
        c_min = hmcode_dic['c_min']
    else:
        c_min = 4.
    k = power_spec_dic['k']
    PS_cold = power_spec_dic['power_cold']
    PS_ax = power_spec_dic['power_axion']

    #########################################
    ############# Cold matter term ##########
    #########################################
    PS_cold_nonlin = func_non_lin_PS_matter(M, k, PS_cold, cosmo_dic, hmcode_dic, cosmo_dic['Omega_db_0'], 
                                            alpha = alpha, eta_given = eta_given, ax_one_halo=False, one_halo_damping = one_halo_damping, 
                                            two_halo_damping = two_halo_damping, concentration_param=concentration_param, full_2h=full_2h, axion_dic=axion_dic)[0]
    
    ########################################
    ##compute everything for one halo term##
    ########################################
    halo_mass_func_arr_2 = func_halo_mass_function(axion_dic['M_int'], k, PS_cold, cosmo_dic, 
                                                    cosmo_dic['Omega_db_0']) # for the integral over reduced array [M_cut, inf]
    dens_profile_cold_arr_2 = func_dens_profile_kspace(axion_dic['M_int'], k, PS_cold, cosmo_dic, hmcode_dic, cosmo_dic['Omega_db_0'], 
                                                       c_min, eta_given=eta_given, axion_dic=axion_dic) # for the integral over reduced array [M_cut, inf]
    dens_profile_ax_arr_2 = func_dens_profile_ax_kspace(k, axion_dic['M_int'], cosmo_dic, power_spec_dic, axion_dic['M_cut'], axion_dic['central_dens'], 
                                                        hmcode_dic, 
                                                        concentration_param=concentration_param, eta_given=False, axion_dic=axion_dic) # this integral in only in the reduced one
    
    ### cross one halo term ###
    #one halo damping, if wanted, you should, because otherwise you have porblems on large scales/small k's
    integrand_arr_one_halo_cross = axion_dic['M_int'][:, None] * axion_dic['M_ax'][:, None] * halo_mass_func_arr_2[:, None]\
                                   * dens_profile_cold_arr_2 * dens_profile_ax_arr_2 # integral over reduced array [M_cut, inf]   
    if one_halo_damping == True:
        one_halo_term_cross = (k/hmcode_dic['k_star'])**4 / (1+(k/hmcode_dic['k_star'])**4) \
                          * integrate.simpson(integrand_arr_one_halo_cross, x = axion_dic['M_int'], axis = 0) \
                          / (func_rho_comp_0(cosmo_dic['Omega_db_0']) 
                          * func_rho_comp_0(cosmo_dic['Omega_ax_0']) * axion_dic['frac_cluster'])
    else:
        one_halo_term_cross = integrate.simpson(integrand_arr_one_halo_cross, x = axion_dic['M_int'], axis = 0) \
                          / (func_rho_comp_0(cosmo_dic['Omega_db_0']) 
                          * func_rho_comp_0(cosmo_dic['Omega_ax_0']) * axion_dic['frac_cluster'])

    #### axion one halo term ###
    #one halo damping, if wanted, you should, because otherwise you have porblems on large scales/small k's
    integrand_arr_one_halo_ax =  axion_dic['M_ax'][:, None]**2 * halo_mass_func_arr_2[:, None] * np.array(dens_profile_ax_arr_2)**2 #for the integral over reduced array [M_cut, inf] 
    if one_halo_damping == True:
        one_halo_term_ax = (k/hmcode_dic['k_star'])**4 / (1+(k/hmcode_dic['k_star'])**4) \
                           * integrate.simpson(integrand_arr_one_halo_ax, x = axion_dic['M_int'], axis = 0) \
                           / (func_rho_comp_0(cosmo_dic['Omega_ax_0']) * axion_dic['frac_cluster'])**2 
    else:
        one_halo_term_ax = integrate.simpson(integrand_arr_one_halo_ax, x = axion_dic['M_int'], axis = 0) \
                           / (func_rho_comp_0(cosmo_dic['Omega_ax_0']) * axion_dic['frac_cluster'])**2 
    

    ########################################################################################
    ##compute all interdients for the two halo term if we want to use the full two halo term
    ########################################################################################
    if full_2h == True:
        ## individual ingredients ##
        halo_mass_func_arr = func_halo_mass_function(M, k, PS_cold, cosmo_dic, 
                                                    cosmo_dic['Omega_db_0']) # for the integral over M ~[0, inf]
        dens_profile_cold_arr = func_dens_profile_kspace(M, k, PS_cold, cosmo_dic, hmcode_dic, cosmo_dic['Omega_db_0'], 
                                                         c_min, eta_given = eta_given, axion_dic=axion_dic) # for the integral over M ~[0, inf]
        halo_bias_arr = func_halo_bias(M, k, PS_cold, cosmo_dic['Omega_db_0'], cosmo_dic['Omega_m_0'], cosmo_dic['Omega_w_0'], cosmo_dic['z'], cosmo_dic['G_a']) # for the integral over M ~[0, inf]
        halo_bias_arr_2 = func_halo_bias(axion_dic['M_int'], k, PS_cold, cosmo_dic['Omega_db_0'], cosmo_dic['Omega_m_0'], cosmo_dic['Omega_w_0'], cosmo_dic['z'], cosmo_dic['G_a']) # for the integral over reduced array [M_cut, inf]
    

        ## cross two halo term integrals ##
        integrand_arr_two_halo_cross = M[:, None] * halo_mass_func_arr[:, None] * halo_bias_arr[:, None] * dens_profile_cold_arr  # integral over M ~[0, inf]
        integrand_arr_two_halo_2_cross = axion_dic['M_ax'][:, None] * halo_mass_func_arr_2[:, None] * halo_bias_arr_2[:, None] * dens_profile_ax_arr_2 # integral over reduced array [M_cut, inf]
        #summand2_cross to take care of nummericals issues of the integral when going to M = 0, see appendix A in https://arxiv.org/abs/2005.00009
        summand2_cross = func_dens_profile_kspace(np.min(M), k, PS_cold, cosmo_dic, hmcode_dic, cosmo_dic['Omega_db_0'], 
                                                  c_min, eta_given = eta_given, axion_dic=axion_dic) \
                        * ( 1 - integrate.simpson(M[:, None] * halo_mass_func_arr[:, None] * halo_bias_arr[:, None], x = M, axis = 0) \
                        / func_rho_comp_0(cosmo_dic['Omega_db_0']) ) 
        factor2_cross = integrate.simpson(integrand_arr_two_halo_cross, x = M, axis = 0) / func_rho_comp_0(cosmo_dic['Omega_db_0']) + summand2_cross
            
        ### axion two halo term integral ##
        integrand_arr_two_halo_ax = axion_dic['M_ax'][:, None] * halo_mass_func_arr_2[:, None] * halo_bias_arr_2[:, None] * dens_profile_ax_arr_2 #for the integral over reduced array [M_cut, inf]
    
    ### cross two halo term ###
    #two halo damping, if wanted
    #if full two halo term is wanted some extra factors in the two halo term to take care of nummerical issues.
    #see appendix A in https://arxiv.org/abs/2005.00009
    if two_halo_damping == True:
        if full_2h == True:
            two_halo_term_cross =  PS_cold* (1-hmcode_dic['f'] * (k/hmcode_dic['k_d'])**hmcode_dic['n_d']/(1+(k/hmcode_dic['k_d'])**hmcode_dic['n_d']))\
                                * factor2_cross[0] \
                                * integrate.simpson(integrand_arr_two_halo_2_cross, x = axion_dic['M_int'], axis = 0)[:] \
                                / (func_rho_comp_0(cosmo_dic['Omega_ax_0']) * axion_dic['frac_cluster']) 
        else:
            two_halo_term_cross =  PS_cold* (1-hmcode_dic['f'] * (k/hmcode_dic['k_d'])**hmcode_dic['n_d']/(1+(k/hmcode_dic['k_d'])**hmcode_dic['n_d']))
        
    else:
        if full_2h == True:
            two_halo_term_cross =  PS_cold\
                                * factor2_cross[0] \
                                * integrate.simpson(integrand_arr_two_halo_2_cross, x = axion_dic['M_int'], axis = 0)[:] \
                                / (func_rho_comp_0(cosmo_dic['Omega_ax_0']) * axion_dic['frac_cluster']) 
        else:
            two_halo_term_cross =  PS_cold

    ### axion two halo term ###
    #two halo damping, if wanted
    #if full two halo term is wanted some extra factors in the two halo term to take care of nummerical issues.
    #see appendix A in https://arxiv.org/abs/2005.00009
    if two_halo_damping == True:
        if full_2h == True:
            two_halo_term_ax =  PS_cold * (1-hmcode_dic['f'] * (k/hmcode_dic['k_d'])**hmcode_dic['n_d']/(1+(k/hmcode_dic['k_d'])**hmcode_dic['n_d']))\
                            * integrate.simpson(integrand_arr_two_halo_ax, x = axion_dic['M_int'], axis = 0)[:]**2 \
                            / (func_rho_comp_0(cosmo_dic['Omega_ax_0']) * axion_dic['frac_cluster'])**2 
        else:
            two_halo_term_ax =  PS_cold * (1-hmcode_dic['f'] * (k/hmcode_dic['k_d'])**hmcode_dic['n_d']/(1+(k/hmcode_dic['k_d'])**hmcode_dic['n_d']))
    else:
        if full_2h == True:
            two_halo_term_ax =  PS_cold\
                            * integrate.simpson(integrand_arr_two_halo_ax, x = axion_dic['M_int'], axis = 0)[:]**2 \
                            / (func_rho_comp_0(cosmo_dic['Omega_ax_0']) * axion_dic['frac_cluster'])**2 
        else:
            two_halo_term_ax =  PS_cold
    
    ##########################################
    ############# Cross term #################
    ##########################################
    PS_total_cross = axion_dic['frac_cluster']*(one_halo_term_cross+two_halo_term_cross) + (1-axion_dic['frac_cluster'])*np.sqrt(PS_ax * PS_cold_nonlin)

    ###########################################
    ############ axion term ###################
    ###########################################
    PS_total_ax = axion_dic['frac_cluster']**2 * (one_halo_term_ax + two_halo_term_ax) \
                  + 2*(1- axion_dic['frac_cluster']) * axion_dic['frac_cluster'] * np.sqrt((one_halo_term_ax + two_halo_term_ax)*PS_ax) \
                  + (1- axion_dic['frac_cluster'])**2 * power_spec_dic['power_axion']  
    

    ##########################################################
    #####stick all together to the total matter non-lin PS####
    ##########################################################
    PS_total_matter = (cosmo_dic['Omega_db_0']/cosmo_dic['Omega_m_0'])**2 * PS_cold_nonlin \
                      + 2*cosmo_dic['Omega_db_0']*cosmo_dic['Omega_ax_0']/cosmo_dic['Omega_m_0']**2 * PS_total_cross \
                      + (cosmo_dic['Omega_ax_0']/cosmo_dic['Omega_m_0'])**2 * PS_total_ax
    
    return PS_total_matter, PS_cold_nonlin, PS_total_cross, PS_total_ax



