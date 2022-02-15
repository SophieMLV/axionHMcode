'''
different scales of interests for the axions
'''

import numpy as np
from scipy import optimize
from scipy.interpolate import interp1d

import variance_functions
import overdensities
import PS_interpolate
import basic_functions_cosmology
import cold_density_profile

#### non-linear scales defined by \sigma=1 and k = \pi / r
def func_R_nonlin(cosmo_dic, k, PS_cold, LCDM=True):
    """
    k is in units of h/Mpc and PS_cold in (Mpc/h)^3
    returns nonlinear scale where sigma=1 in Mpc/h
    """
    if LCDM == True:
        delta_c = overdensities.func_delta_c_LCDM(cosmo_dic)
    else:
        delta_c = overdensities.func_delta_c(cosmo_dic)
    def find_root(R):
        return variance_functions.func_sigma_r(R, k, PS_cold) - 1.0#1.686 #delta_c

    R_nonlin = optimize.brentq(find_root, 1e-3, 10.)
    return R_nonlin

def func_k_nonlin(cosmo_dic, power_spec_dic, hubble=True):
    """
    k_sigma in h/Mpc and PK_sigma in (Mpc/h)^3
    returns the non-linear k either in 1/Mpc or h/Mpc depending on the parameter hubble
    """ 
    # in units h/Mpcs and (Mpc/h)^3
    k_interpolate, PS_interploate = PS_interpolate.func_interploate_PS(1e-5, power_spec_dic['k'], power_spec_dic['cold'] , cosmo_dic, cosmo_dic['Omega_db_0']) # in h/Mpc and (Mpc/h)^3
    R_nonlin = func_R_nonlin(cosmo_dic, k_interpolate, PS_interploate , LCDM=True) # in Mpc/h

    k_nonlin = np.pi / (R_nonlin) # in h/Mpc
    if hubble == True:
        return k_nonlin
    else:
        return k_nonlin * cosmos_dic['h']
    
    

    
#### half mode scales ####

def func_transfer_ratio(cosmo, transfer_dic, transfer_dic_LCDM):
    # be carefull transfer_dic['k'] is in units of h/Mpc and the tranfer_function in (Mpc)^2
    transfer_total_matter_LCDM_func = interp1d(transfer_dic_LCDM['k'], transfer_dic_LCDM['transfer_total'], kind='cubic')
 
    return transfer_dic['k'] , transfer_dic['transfer_total'] / transfer_total_matter_LCDM_func(transfer_dic['k']) #transfer_dic_LCDM['transfer_total']

def func_half_mode(cosmo, transfer_dic, transfer_dic_LCDM):
    # be carefull transfer_dic['k'] is in units of h/Mpc and the tranfer_function in (Mpc)^2
    
    k_in_func, transfer_ratio_in_func = func_transfer_ratio(cosmo, transfer_dic, transfer_dic_LCDM) # k in h/Mpc

    transfer_ratio_half_value =  transfer_ratio_in_func[0] - 0.5 * (transfer_ratio_in_func[0] - transfer_ratio_in_func[-5])
    
    func_ratio_of_k = interp1d(k_in_func, transfer_ratio_in_func, kind='cubic')
    
    def func_root(k):
        return func_ratio_of_k(k) - transfer_ratio_half_value
        
    half_mode_value = optimize.brentq(func_root, transfer_dic['k'][0], transfer_dic['k'][-1]) # in units ofh1/Mpc
    return half_mode_value #, transfer_ratio_half_value


#### jeans scale in a pure axion DM universe ####
def func_jeans_scale(cosmo_dic):
    """
    Jeans scale as in https://arxiv.org/abs/1510.07633 eq. 101
    m_ax in eV
    returns jeans scale in h/Mpc
    """
    cosmos_func = cosmo_dic.copy()
    omega_ax = cosmos_func['omega_ax_0']
    m_ax = cosmos_func['m_ax']
    z = cosmos_func['z']
    return 66.5 * (1+z)**(-1/4) * (omega_ax/0.12)**(1/4) * (m_ax/1e-22)**(1/2) / cosmo_dic['h']



##### halo jeans scale in a mixed DM universe ####
def func_halo_jeans_kscale(M, cosmo_dic, power_spec_dic):
    """
    halo jeans scale form https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.85.1158 eq4/5
    note, that we use a slightly different definition of r= \pi /k
    M in solar_mass/h
    return halo jeans scale in units of h/Mpc
    """
    
    #prefactor = np.pi / ( (np.pi/66.5 )**(4/3) * 200**(-2/9)  * 0.12**(1/9) * 3**(2/9) * 4**(1/9) * (2.8e11)**(1/9) * ((np.log(11) - 10/11)/100)**(1/3) * (1e10)**(-1/9) )
    prefactor = 1424.6287862379945
    
    z = cosmo_dic['z']
    m_ax = cosmo_dic['m_ax']
    omega_cold = cosmo_dic['omega_db_0']
    concentration = cold_density_profile.func_conc_param_2(M, power_spec_dic['k'], power_spec_dic['cold'], cosmo_dic, cosmo_dic['Omega_db_0'], LCDM = True)
    
    form_factor = cold_density_profile.func_for_norm_factor(concentration) 
    M_in_formula = M / cosmo_dic['h'] # note that we need solar_mass units for M in the formula
    return prefactor * (m_ax/1e-22)**(2/3) * (form_factor/concentration**2 *100/((np.log(11) - 10/11)))**(-1/3) * (M_in_formula/1e10)**(1/9) * (omega_cold/0.12)**(2/9) / cosmo_dic['h'] * (1+z)**(-1/4) 
    
    
    
