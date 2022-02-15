"""
some basic function needed for the non-linear power spectrum
"""

#import packages
import numpy as np
from astropy import constants as const



def func_E_z(z, cosmo_dic):
    """
    return time dependence of the Hubble constants
    """
    E2 = cosmo_dic['Omega_m_0']*(1+z)**3 + cosmo_dic['Omega_w_0']
    return np.sqrt(np.abs(E2))
 
def func_H_z(cosmo_dic):
    """
    H(z) is given in units of H0 eg km/(s * Mpc)
    """
    z = cosmo_dic['z']
    return 100*cosmo_dic['h'] * func_E_z(z, cosmo_dic)

def func_Omega_comp_z(cosmo_dic, Omega_0):
    """
    returns matter componenet specified by Omega_0
    density parameter as a function of redshift
    """
    z = cosmo_dic['z']
    return Omega_0 * (1+z)**3 / func_E_z(z, cosmo_dic)**2

def func_rho_comp_0(Omega_0):
    """
    present density of component Omega_0
    retruns demsity in units of solar_mass/Mpc^3*h^2
    """
    grav_const = const.G.to('km**2*Mpc/(Msun*s**2)').value
    rho_crit_0 = 3. * 100**2 / (8. * np.pi * grav_const)
    
    return Omega_0  * rho_crit_0


def func_R_M(M, cosmo_dic, Omega_0):
    """
    M in solar_mass/h, where M is matter of the type Omega_0
    returns R in Mpc/h
    """
    R = (3 * M / (4 * np.pi * func_rho_comp_0(Omega_0)) )**(1./3.)
    return R  

