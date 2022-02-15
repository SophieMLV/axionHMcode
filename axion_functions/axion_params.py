#packages needed
import numpy as np
#own pachages
from cosmology.basic_cosmology import *
from axion_cutoff_mass import *
from halo_model.halo_bias import *
from halo_model.halo_mass_function import *
from halo_model.axion_density_profile import *


def func_axion_param_dic(M, cosmo_dic, power_spec_dic_sigma, eta_given=False):
    """
    generate dictionary with parameters for axions
    M_halo is the mass of the cold matter halo in solar_mass/h
    """
    #generate dictionary
    axion_param_dic = {}
    #cut of mass. Below this cold halo mass no axion halo exists
    axion_param_dic['M_cut'] = func_cut_mass_axion_halo(cosmo_dic, power_spec_dic_sigma)
    #cold halo masses for which an axion halo exists
    axion_param_dic['M_int'] = np.geomspace(axion_param_dic['M_cut'], np.max(M), num=len(M))
    #central density parameter, to ensure, that the axion halo has the correct mass
    axion_param_dic['central_dens'] = func_central_density_param(axion_param_dic['M_int'], cosmo_dic, power_spec_dic_sigma, axion_param_dic['M_cut'], eta_given=eta_given)
    
    #for some halo mass no central density parameter can be found, ie there is no axion halo for this halo mass and M_int must be reduced
    axion_param_dic['M_int'] = axion_param_dic['M_int'] * np.where(np.array(axion_param_dic['central_dens']) <= 0, 0, 1) #set mass to zero, if central densit param is zero
    axion_param_dic['M_int'] = axion_param_dic['M_int'][axion_param_dic['M_int'] != 0.0] #delete all zero halo masses
    axion_param_dic['central_dens'] = np.array(axion_param_dic['central_dens'])[np.array(axion_param_dic['central_dens']) != 0.0] #delete all zero central density params
    
    #axion halo mass, given my the cosmic abundance
    axion_param_dic['M_ax'] = cosmo_dic['Omega_ax_0']/cosmo_dic['Omega_db_0'] * axion_param_dic['M_int']
    
    #not all axions cluster. Compute clustered fraction by f = 1/\rho_ax * int_Mcut^\inf n(M)*b(M)*M_ax(M) dM, see masterthesis eq.
    k = power_spec_dic_sigma['k']
    PS_cold = power_spec_dic_sigma['cold']
    integrand_arr = func_halo_mass_function(axion_param_dic['M_int'], k, PS_cold , cosmo_dic, cosmo_dic['Omega_db_0'], cosmo_dic['Omega_db_0']) * \
                    axion_param_dic['M_ax'] * \
                    func_halo_bias(axion_param_dic['M_int'], k, PS_cold, cosmo_dic, cosmo_dic['Omega_db_0'])
    axion_param_dic['frac_cluster'] = integrate.simps(y= integrand_arr, x = axion_param_dic['M_int'])/ func_rho_comp_0(cosmo_dic['Omega_ax_0'])
    
    return axion_param_dic

