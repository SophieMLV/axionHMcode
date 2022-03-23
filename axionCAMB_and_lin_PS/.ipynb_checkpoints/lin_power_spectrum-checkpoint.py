"""
functions to compute the linear power spectrum 
form the tranfer functions given from axionCAMB
"""
import numpy as np
import sys
from PS_interpolate import *

def primordial_PS(k, initial_amplitude, k_pivot, scalar_index, hubble):
    """
    Input: Initial amplitude is dimensionless, k in h/Mpc and k_pivot in 1/Mpc.
    Returns the dimensionless primordial power spectrum.
    """
    return initial_amplitude * (k * (hubble) / k_pivot) ** (scalar_index - 1.0)


def transfer_to_PS(k, transfer_function, cosmo_dic):
    """
    k in in h/Mpc
    The transfer function from axionCAMB is given by the 'normal' transfer function 
    divided by k^2 with k in 1/Mpc, so has units 1/Mpc^2
    returns the linar power spectrum associated with the given tranfer function 
    in units Mpc^3/h^3 at k in h/Mpc.
    """
    return transfer_function ** 2 * (k * cosmo_dic['h']) ** 4 \
           * primordial_PS(k, cosmo_dic['As'], cosmo_dic['k_piv'], cosmo_dic['ns'], cosmo_dic['h']) \
           * 2 * np.pi ** 2 / k ** 3


def load_transfer_from_file(file_path):
    """
    Loads transfer function output file in file_path:
    Load all different transfer functions in the order as metioned in components_dic
    k will be in units h/Mpc and the transfer function by CAMB is divided by k^2 so in units (Mpc)^2.
    """
    components_dic = {'k': (0,), 'CDM': (1,), 'baryon': (2,), 'photon': (3,), 'massless neutrino': (4,),
                      'massive neutrino': (5,),
                      'axion': (6,), 'growth rate': (7,), 'total': (8,)}
    transfer_dic = {
        'k': np.loadtxt(file_path, unpack=True, usecols=components_dic['k']),
        'transfer_total': np.loadtxt(file_path, unpack=True, usecols=components_dic['total']),
        'growth_rate': np.loadtxt(file_path, unpack=True, usecols=components_dic['growth rate']),
        'transfer_CDM': np.loadtxt(file_path, unpack=True, usecols=components_dic['CDM']),
        'transfer_baryon': np.loadtxt(file_path, unpack=True, usecols=components_dic['baryon']),
        'transfer_axion': np.loadtxt(file_path, unpack=True, usecols=components_dic['axion']),
        'transfer_massless-neutrino': np.loadtxt(file_path, unpack=True,
                                                 usecols=components_dic['massless neutrino']),
        'transfer_massive-neutrino': np.loadtxt(file_path, unpack=True,
                                                usecols=components_dic['massive neutrino']),
    }
    return transfer_dic

def func_power_spec_dic(file_path, cosmo_dic):
    """
    crate a dictionary with the linear power spectra computed from the different
    tranfer functions from axionCAMB (stored in file in file_path) as a function of k. 
    The units for k is h/Mpc and for the PSs (Mpc/h)^3.
    The folowing linear power spectra will be computed:
    total, CDM, baryon, cold matter, axion
    """ 
    #call the tranfer_dic function to get the different tranfer functions
    transfer_dic = load_transfer_from_file(file_path)
    
    #create PS dictionary
    power_spec_dic = {}
    power_spec_dic['k'] = transfer_dic['k'] 
    power_spec_dic['power_total'] = transfer_to_PS(transfer_dic['k'], transfer_dic['transfer_total'], cosmo_dic) 
    power_spec_dic['power_CDM'] = transfer_to_PS(transfer_dic['k'], transfer_dic['transfer_CDM'], cosmo_dic) 
    power_spec_dic['power_baryon'] = transfer_to_PS(transfer_dic['k'], transfer_dic['transfer_baryon'], cosmo_dic) 
    #for the cold matter PS the cold matter transfer function is need: T_cold = (T_cdm*Omega_cdm + T_b*Omega_b)/(Omega_cdm+Omega_b)
    power_spec_dic['cold'] = transfer_to_PS(transfer_dic['k'],
                                            (cosmo_dic['Omega_b_0']*transfer_dic['transfer_baryon'] + \
                                             cosmo_dic['Omega_d_0']*transfer_dic['transfer_CDM'])/cosmo_dic['Omega_db_0'], 
                                             cosmo_dic)
    power_spec_dic['power_axion'] = transfer_to_PS(transfer_dic['k'], transfer_dic['transfer_axion'], cosmo_dic)
    return power_spec_dic

def func_power_spec_interp_dic(power_spec_dic, cosmo_dic):
    """
    crate a dictionary with the linear power spectra from func_power_spec_dic
    on a larger k range. These PSs are needed to compute the variance correctly.
    The k range depends on the smallest halo mass given in cosmo_dic
    The units for k is h/Mpc and for the PSs (Mpc/h)^3.
    The following linear power spectra will be computed:
    total, CDM, baryon, cold matter, axion
    """ 
    #create PS dictionary
    power_spec_interp_dic = {}
    power_spec_interp_dic['k'], _ = func_PS_interpolate_M(cosmo_dic['M_min'], power_spec_dic['k'], power_spec_dic['cold'], cosmo_dic, cosmo_dic['Omega_db_0'])
    power_spec_interp_dic['power_total'] = func_PS_interpolate(power_spec_interp_dic['k'], power_spec_dic['k'], power_spec_dic['power_total'])
    power_spec_interp_dic['power_CDM'] = func_PS_interpolate(power_spec_interp_dic['k'], power_spec_dic['k'], power_spec_dic['power_CDM'])
    power_spec_interp_dic['power_baryon'] = func_PS_interpolate(power_spec_interp_dic['k'], power_spec_dic['k'], power_spec_dic['power_baryon'])  
    power_spec_interp_dic['cold'] = func_PS_interpolate(power_spec_interp_dic['k'], power_spec_dic['k'], power_spec_dic['cold'])    
    power_spec_interp_dic['power_axion'] = func_PS_interpolate(power_spec_interp_dic['k'], power_spec_dic['k'], power_spec_dic['power_axion'])
    return power_spec_interp_dic