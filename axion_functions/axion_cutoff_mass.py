"""functions for axions and when they cluster, ie evolve non-linear"""

import numpy as np
from scipy import optimize

from cosmology.basic_cosmology import *
from cosmology.overdensities import *
from cosmology.variance import *
from halo_model.cold_density_profile import *


    
def func_halo_jeans_kscale(M, cosmo_dic, power_spec_dic_sigma):
    """
    halo jeans scale form https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.85.1158 eq4/5
    note, that we use a slightly different definition of r= \pi /k
    M in solar_mass/h
    return halo jeans scale in units of h/Mpc
    """
    prefactor = 1424.6287862379945 # computed directly, see matherthesis TBC for fomula
    
    z = cosmo_dic['z']
    m_ax = cosmo_dic['m_ax']
    omega_cold = cosmo_dic['omega_db_0']
    concentration = func_conc_param(M, power_spec_dic_sigma['k'], power_spec_dic_sigma['cold'], 
                                    cosmo_dic, cosmo_dic['Omega_db_0'])
    form_factor = func_for_norm_factor(concentration) 
    M_in_formula = M / cosmo_dic['h'] # note that we need solar_mass units for M in the formula
    return prefactor * (m_ax/1e-22)**(2/3) * (form_factor/concentration**2 *100/((np.log(11) - 10/11)))**(-1/3) \
           * (M_in_formula/1e10)**(1/9) * (omega_cold/0.12)**(2/9) / cosmo_dic['h'] * (1+z)**(1/3) 
    
    
def func_jeans_virial_ratio(M, cosmo_dic, power_spec_dic_sigma):
    """
    M is given in solar_mass/h
    returns the ratio between the halo jeans scale and the virial radius
    """
    halo_jeans_scale = np.pi / func_halo_jeans_kscale(M, cosmo_dic, power_spec_dic_sigma) # in Mpc/h
    virial_radius = func_r_vir(M, cosmo_dic, cosmo_dic['Omega_db_0']) # in Mpc/h
    
    return halo_jeans_scale/virial_radius


def func_cut_mass_axion_halo(cosmo_dic, power_spec_dic_sigma):
    """
    axions form a halo if the mass of the host cold matter halo
    satisfies r_hJ / r_vir = 1
    compute this mass here
    M is in solar_mass/h 
    returns minimal mass where axions are put intto halos in solar_mass/h
    """
    def func_find_root(m):
        return func_jeans_virial_ratio(m, cosmo_dic, power_spec_dic_sigma) -1
    
    Mass_min = optimize.brentq(func_find_root, 1e1, 1e18)
    return Mass_min

