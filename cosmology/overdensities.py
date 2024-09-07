"""
functions for the critical overdensity in spherical/ellipsoidal collaps and halo overdensity
"""

import numpy as np
from scipy import integrate
from basic_cosmology import func_E_z, func_Omega_comp_z, func_rho_comp_0


def func_D_z_unnorm(z, cosmo_dic):
    """
    returns unnormalised grwoth function 
    """
    def integrand(x):
        return (1+x)/func_E_z(x, cosmo_dic)**3
    
    factor = 5*cosmo_dic['Omega_m_0']/2
    D = factor * func_E_z(z, cosmo_dic) * integrate.quad(integrand, z, np.inf)[0]
    return D

def func_D_z_unnorm_int(z, cosmo_dic):
    f = lambda y, x: func_E_z(x, cosmo_dic)/(1+x)*(1+y)/func_E_z(y, cosmo_dic)**3
    G =  5 * cosmo_dic['Omega_m_0'] / 2 * integrate.dblquad(f, z, 10000, lambda x:x, 10000)[0] # now with Numba trapz
    return G

def func_D_z_norm(z, cosmo_dic):
    """
    returns normlised grwoth function, ie D(0)= 1
    this is used to scale the power spectrum for diffrent z's
    """
    normalisation = func_D_z_unnorm(0., cosmo_dic)
    growth = func_D_z_unnorm(z, cosmo_dic)
    
    return growth/normalisation

def func_delta_c_old(cosmo_dic):
    """
    returns critical denity for spherical/ellepsoidal collapse for LCDM cosmos
    """          
    return 1.686
    
def func_Delta_vir_old(cosmo_dic, Omega_0):
    """
    halo overdensity for LCDM comos,
    make the change, that only matter of the type 
    Omega_0 is take into accound for the overdensity
    """
    #x = func_Omega_comp_z(cosmo_dic, Omega_0) -1
    x = func_Omega_comp_z(cosmo_dic, cosmo_dic['Omega_m_0']) -1
    return (18*np.pi**2 + 82*x - 39*x**2) / (x+1)

def f1_delta_c(x, y):
    return -0.0069 - 0.0208*(1-x) + 0.0312*(1-x)**2 + 0.0021*(1-y)

def f2_delta_c(x, y):
    return 0.0001 - 0.0647*(1-x) - 0.0417*(1-x)**2 + 0.0646*(1-y)

def f3_delta_vir(x, y):
    return -0.79 - 10.17*(1-x) + 2.51*(1-x)**2 + 6.51*(1-y)

def f4_delta_vir(x, y):
    return -1.89 + 0.38*(1-x) + 18.8*(1-x)**2 - 15.87*(1-y)

def func_delta_c(cosmo_dic):
    """
    returns critical denity for spherical/ellepsoidal collapse for LCDM cosmos
    """          
    Omega_a = func_Omega_comp_z(cosmo_dic, cosmo_dic['Omega_m_0'])
    g = func_D_z_unnorm(cosmo_dic['z'], cosmo_dic)
    G = cosmo_dic['G']
    a = 1/(1+cosmo_dic['z'])
    f1 = f1_delta_c(g/a, G/a)
    f2 = f2_delta_c(g/a, G/a)
    f_frac = cosmo_dic['omega_ax_0']/cosmo_dic['omega_m_0']
    func_delta = 1.686*(1-0.041*f_frac)*(1+f1*np.log10(Omega_a)**1+f2*np.log10(Omega_a)**0)
    return func_delta
    
def func_Delta_vir(cosmo_dic, Omega_0):
    """
    halo overdensity for LCDM comos,
    make the change, that only matter of the type 
    Omega_0 is take into accound for the overdensity
    """
    Omega_a = func_Omega_comp_z(cosmo_dic, cosmo_dic['Omega_m_0'])
    g = func_D_z_unnorm(cosmo_dic['z'], cosmo_dic)
    G = cosmo_dic['G']
    a = 1/(1+cosmo_dic['z'])
    f3 = f3_delta_vir(g/a, G/a)
    f4 = f4_delta_vir(g/a, G/a)
    f_frac = cosmo_dic['omega_ax_0']/cosmo_dic['omega_m_0']
    func_Delta_vir = 177.7*(1+0.763*f_frac)*(1+f3*np.log10(Omega_a)**1+f4*np.log10(Omega_a)**2)
    return func_Delta_vir

def func_r_vir(M, cosmo_dic, Omega_0):
    """
    M in solar_mass/h where M is matter of the type Omega_0
    returns comoving virial radius in Mpc/h
    """
    Delta_vir = func_Delta_vir(cosmo_dic, Omega_0)
    return (3. * M / (4. * np.pi * func_rho_comp_0(Omega_0) * Delta_vir ))**(1./3.)