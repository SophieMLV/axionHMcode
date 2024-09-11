"""functions for axions and when they cluster, ie evolve non-linear"""

import numpy as np
from cosmology.overdensities import func_r_vir
from halo_model.cold_density_profile import func_conc_param, func_for_norm_factor
from scipy.interpolate import RegularGridInterpolator
import warnings
    
def func_halo_jeans_kscale(M, cosmo_dic, power_spec_dic, axion_dic=None):
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
    concentration = func_conc_param(M, power_spec_dic['k'], power_spec_dic['power_cold'], 
                                    cosmo_dic, cosmo_dic['Omega_db_0'], axion_dic=axion_dic)
    form_factor = func_for_norm_factor(concentration) 
    M_in_formula = M / cosmo_dic['h'] # note that we need solar_mass units for M in the formula
    return prefactor * (m_ax/1e-22)**(2/3) * (form_factor/concentration**2 *100/((np.log(11) - 10/11)))**(-1/3) \
           * (M_in_formula/1e10)**(1/9) * (omega_cold/0.12)**(2/9) / cosmo_dic['h'] * (1+z)**(1/3) 
    
    
def func_jeans_virial_ratio(M, cosmo_dic, power_spec_dic, axion_dic=None):
    """
    M is given in solar_mass/h
    returns the ratio between the halo jeans scale and the virial radius
    """
    halo_jeans_scale = np.pi / func_halo_jeans_kscale(M, cosmo_dic, power_spec_dic, axion_dic=axion_dic) # in Mpc/h
    virial_radius = func_r_vir(cosmo_dic['z'], M, cosmo_dic['Omega_db_0'], cosmo_dic['Omega_m_0'], 
                               cosmo_dic['Omega_w_0'], cosmo_dic['G_a']) # in Mpc/h
    
    return halo_jeans_scale/virial_radius

def func_cut_mass_axion_halo(cosmo_dic, power_spec_dic, axion_dic=None):
    """
    axions form a halo if the mass of the host cold matter halo
    satisfies r_hJ / r_vir = 1
    compute this mass here
    M is in solar_mass/h 
    returns minimal mass where axions are put intto halos in solar_mass/h
    """
    def func_find_root(m):
        return func_jeans_virial_ratio(10**(m), cosmo_dic, power_spec_dic, axion_dic=axion_dic) -1
    
    Mass_min = optimize.brentq(func_find_root, 7, 17)

def func_beta2(cosmo_dic, power_spec_dic_sigma, axion_dic = None):
    """
    returns beta2 parameter in axion mass-cold mass relation Max(Mc) parametrization (steepness of suppression compared to cosmic average Ωa/Ωc Mc)
    """
    # Define the data
    z_values = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
    f_values = np.array([0.01, 0.1, 0.2, 0.3])
    
    # Data matrix for β2 values corresponding to each (z, f) pair
    beta2_matrix = np.array([
        [1.27, 0.98, 0.61, 0.30],  # z = 1.0
        [1.22, 0.87, 0.59, 0.32],  # z = 2.0
        [1.32, 0.93, 0.66, 0.40],  # z = 3.0
        [1.37, 0.96, 0.71, 0.45],  # z = 4.0
        [1.43, 1.04, 0.75, 0.51],  # z = 5.0
        [1.46, 1.02, 0.79, 0.58],  # z = 6.0
        [1.48, 1.04, 0.77, 0.63],  # z = 7.0
        [1.52, 1.04, 0.78, 0.62]   # z = 8.0
    ])
    
    # Create an interpolation function
    ip = RegularGridInterpolator((z_values, f_values), beta2_matrix, bounds_error=False, fill_value=None) # values outside the domain are extrapolated

    # Function to get β2 value for any (z, f)
    def get_beta2(z, f):
        if z < z_values.min() or z > z_values.max():
            warnings.warn("Warning. Redshift z = {:.2f} is not in the range [1.0, 8.0] where beta2 is defined. Extrapolate".format(z))
        if f < f_values.min() or f > f_values.max():
            warnings.warn("Warning. Axion fraction f = {:.2f} is not in the range [0.01, 0.3] where beta2 is defined. Extrapolate".format(f))
        p = [z, f]
        return ip([p])[0]
    z = cosmo_dic['z']
    f = cosmo_dic['omega_ax_0']/cosmo_dic['omega_m_0']
    beta2 = get_beta2(z, f)
    return beta2