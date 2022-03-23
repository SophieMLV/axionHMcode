"""functions for the density profile of axions"""

import numpy as np
from scipy import optimize, integrate
from scipy.interpolate import interp1d
from astropy import constants as const

from cosmology.variance import *
from cosmology.overdensities import *
from cosmology.basic_cosmology import *
from HMcode_params import *
from cold_density_profile import *




def func_core_radius(M, cosmo_dic):
    """
    M in solar_mass/h
    computes the core radius of the soliton given 
    as in https://arxiv.org/abs/2007.08256 eq. 8
    returns r in Mpc/h
    """
    grav_const = const.G.to('m**3/(Msun*s**2)').value
    M_tot = (1 + cosmo_dic['omega_ax_0']/cosmo_dic['omega_db_0']) * M/cosmo_dic['h'] # in solar_mass
    r_vir = func_r_vir((1 + cosmo_dic['omega_ax_0']/cosmo_dic['omega_db_0']) * M, cosmo_dic, cosmo_dic['Omega_m_0']) / cosmo_dic['h'] * 3.086e+22 # in m
    v_vir = np.sqrt(grav_const*M_tot/r_vir) # in m/s
    
    h_bar = const.hbar.value
    m_ax = cosmo_dic['m_ax'] * 1.78266269594644e-36 # in kg
    r_core = 2 * np.pi * h_bar / (7.5 * m_ax * v_vir) # in m
    
    return r_core / (3.086e+22) * cosmo_dic['h'] * (1+cosmo_dic['z'])**(-1./2.) # in Mpc/h



def func_rho_soliton(r, M, cosmo_dic, rho_central_param):
    """
    soliton profile for axions as in https://arxiv.org/abs/1407.7762 eq.3
    but with core radius as in func_core_radius
    r in Mpc/h, M_vir in solar_mass/h and m_ax in eV
    the rho_central_param scales the central density, this is needed
    for the complete axion density profile, see my masterthesis eq TBC
    returns the soliton denity profile in solar_mass/pc^3 * h^2
    """
    m_ax = cosmo_dic['m_ax']
    z = cosmo_dic['z']
    A = (1+z) * 0.019 * (m_ax/1e-22)**(-2)
    x_c = func_core_radius(M, cosmo_dic) * 1e3 / cosmo_dic['h'] #in the formula we need units kpc
    r_in_formula = r * 1e3 / cosmo_dic['h']  #in the formula we need units kpc
    
    if isinstance(M, (int, float)) == True:
        return A * rho_central_param / ( x_c**4 * (1 + 0.091 * (r_in_formula/x_c)**2)**8 )\
               * cosmo_dic['h']**2 * 1e18 #transform from solar_mass/pc^3 to solar_mass/pc^3 * h^2
    else:
        return A * rho_central_param / ( np.outer(x_c, np.ones(len(r))) **4 * (1 + 0.091 * np.outer(1/x_c, r_in_formula)**2)**8 ) \
               * cosmo_dic['h']**2 * 1e18 #transform from solar_mass/pc^3 to solar_mass/pc^3 * h^2hubble units

    


def func_dens_profile_ax(r_arr, M, cosmo_dic, power_spec_dic_sigma, M_cut, rho_central_param, eta_given=False):
    """
    r_arr in Mpc/h, M and M_cut in solar_mass/h
    returns the axion density profile
    with a solition core and a NFW profile in the 
    outer region and with the free patameter 
    rho_central_param. This free parameter is set such that
    we get the correct mass of the soliton halo, 
    see func_central_density_param
    the density profile has units solar_mass/Mpc^3 * h^2
    see masterthesis sec. 5.2.3.
    """
    #distinguish whether M is an array or a scalar
    if isinstance(M, (int, float)) == True:
        #there is no axion halo, if the cold halo is below a cut-off
        if rho_central_param == 0 or M_cut > M:
            if isinstance(r_arr, (int, float)) == True:
                return 0.0
            else:
                return np.zeros(len(r_arr))
        else:
            hmcode_params = HMCode_param_dic(cosmo_dic, power_spec_dic_sigma['k'], power_spec_dic_sigma['cold'])
            NFW = cosmo_dic['omega_ax_0']/cosmo_dic['omega_db_0'] * \
                          NFW_profile(M, r_arr, power_spec_dic_sigma['k'], power_spec_dic_sigma['cold'], cosmo_dic, 
                                      hmcode_params, cosmo_dic['Omega_db_0'], cosmo_dic['Omega_db_0'], eta_given = eta_given)
            soliton = func_rho_soliton(r_arr, M, cosmo_dic, rho_central_param)
            
            idx_arr = np.argwhere(np.diff(np.sign(NFW - soliton))).flatten() #found the intersection points
            if len(idx_arr)<=0:
                return soliton
            else:
                return np.where(r_arr > r_arr[idx_arr[-1]], NFW, soliton)
    
    else:
        return_arr = []
        hmcode_params = HMCode_param_dic(cosmo_dic, power_spec_dic_sigma['k'], power_spec_dic_sigma['cold'])
        for idx, m in enumerate(M):
            if rho_central_param[idx] == 0 or M_cut > m:
                if isinstance(r_arr, (int, float)) == True:
                    return_arr.append(0.0)
                else:
                    return_arr.append(np.zeros(len(r_arr)))
            else:
                NFW = cosmo_dic['omega_ax_0']/cosmo_dic['omega_db_0'] * \
                              NFW_profile(m, r_arr, power_spec_dic_sigma['k'], power_spec_dic_sigma['cold'], cosmo_dic, 
                                          hmcode_params, cosmo_dic['Omega_db_0'], cosmo_dic['Omega_db_0'], eta_given = eta_given)
                soliton = func_rho_soliton(r_arr, m, cosmo_dic, rho_central_param[idx])
                
                idx_arr = np.argwhere(np.diff(np.sign(NFW - soliton))).flatten() #found the intersection points
                if len(idx_arr)<=0:
                    return_arr.append(soliton)
                else:
                    return_arr.append(np.where(r_arr > r_arr[idx_arr[-1]], NFW, soliton))
        return return_arr
                
        
def func_ax_halo_mass(M, cosmo_dic, power_spec_dic_sigma, M_cut, rho_central_param, eta_given=False):
    """
    M and M_cut in solar_mass/h
    The free parameter rho_central_param is set such that
    we get the correct mass of the soliton halo, 
    see func_central_density_param
    returns the axion halo mass by integrating the halo 
    density profile in units of solar_mass/h
    """
    #distinguish whether M is an array or a scalar
    if isinstance(M, (int, float)) == True:
        r_vir = func_r_vir(M, cosmo_dic, cosmo_dic['Omega_db_0'])
        r_arr = np.geomspace(1e-15, r_vir, num=1000)
        integrand = func_dens_profile_ax(r_arr, M, cosmo_dic, power_spec_dic_sigma, M_cut, rho_central_param, eta_given=eta_given) * r_arr**2
        return 4 * np.pi * integrate.simps(y=integrand, x = r_arr)
    else:
        return [4 * np.pi * integrate.simps(y=func_dens_profile_ax(np.geomspace(1e-15, func_r_vir(M[i], cosmo_dic, cosmo_dic['Omega_db_0']), num=1000), 
                                                                              M[i], cosmo_dic, power_spec_dic_sigma, M_cut, rho_central_param[i], eta_given=eta_given) * \
                                            np.geomspace(1e-15, func_r_vir(M[i], cosmo_dic, cosmo_dic['Omega_db_0']), num=1000)**2, 
                                            x=np.geomspace(1e-15, func_r_vir(M[i], cosmo_dic, cosmo_dic['Omega_db_0']), num=1000)) for i in range(len(M))]
        

def func_central_density_param(M, cosmo_dic, power_spec_dic_sigma, M_cut, eta_given=False):
    """
    M and M_cut in solar_mass/h
    The central density of the soliton profile 
    has to be change in such a way that the total
    mass of the axion halo matches the abundance,
    ie M_ax_halo = Omega_ax/Omega_cold * M_cold_halo
    """
    #distinguish whether M is an array or a scalar
    if isinstance(M, (int, float)) == True:
        if M < M_cut:
            return np.array(0.0)
        
        else:
            r_c = func_core_radius(M, cosmo_dic) 
            hmcode = HMCode_param_dic(cosmo_dic, power_spec_dic_sigma['k'], power_spec_dic_sigma['cold'])
            
            #need a gues to find the correct central_dens_param:
            #guess is set via Omega_ax/Omega_cold * M = int_0_rvir \rho *r^2 dr
            #so we need the soliton and NFW part
            def integrand_ax(x):
                return func_dens_profile_ax(x, M, cosmo_dic, power_spec_dic_sigma, M_cut, 1., eta_given=eta_given)*x**2
            integral_soliton = integrate.quad(integrand_ax, 0, r_c)[0]
            
            r_arr = np.geomspace(1e-10 , 2*r_c, 1000)
            integrand_cold = NFW_profile(M, r_arr, power_spec_dic_sigma['k'], power_spec_dic_sigma['cold'], cosmo_dic, hmcode, cosmo_dic['Omega_db_0'], 
                                         cosmo_dic['Omega_db_0'], eta_given = eta_given) \
                             *r_arr**2 * cosmo_dic['Omega_ax_0']/cosmo_dic['Omega_db_0']
            integral_NFW = integrate.simps(y=integrand_cold, x = r_arr)
            
            guess = (M +  integral_NFW) / integral_soliton
            
            #find the central density parameter by solving the eq: 
            #M_ax_halo = Omega_ax/Omega_cold * M_cold_halo
            def func_find_root(dens):
                return func_ax_halo_mass(M, cosmo_dic, power_spec_dic_sigma, M_cut, dens, eta_given=eta_given) - cosmo_dic['Omega_ax_0']/cosmo_dic['Omega_db_0'] * M
            dens_param = optimize.root(func_find_root, x0 = guess).x
            #sometimes the solution is not really a solution,
            #so set than the central density paameter to zero, ie no solution can be found
            if np.abs(guess - dens_param) > 100.:
                return 0.
            else:
                return float(dens_param)
    else:
        dens_param_arr = []
        r_c = func_core_radius(M, cosmo_dic)
        hmcode = HMCode_param_dic(cosmo_dic, power_spec_dic_sigma['k'], power_spec_dic_sigma['cold'])
        for idx, m in enumerate(M):
            if m < M_cut:
                dens_param_arr.append(0.)
            else:
                #need a gues to find the correct central_dens_param:
                #guess is set via Omega_ax/Omega_cold * M = int_0_rvir \rho *r^2 dr
                #so we need the soliton and NFW part
                def integrand_ax(x):
                    return func_dens_profile_ax(x, m, cosmo_dic, power_spec_dic_sigma, M_cut, 1., eta_given=eta_given)*x**2
                integral_soliton = integrate.quad(integrand_ax, 0, r_c[idx])[0]
                
                r_arr = np.geomspace(1e-15 , r_c[idx], 1000)
                integrand_cold = NFW_profile(m, r_arr, power_spec_dic_sigma['k'], power_spec_dic_sigma['cold'], cosmo_dic, hmcode, cosmo_dic['Omega_db_0'], 
                                             cosmo_dic['Omega_db_0'], eta_given = eta_given)*r_arr**2 * cosmo_dic['Omega_ax_0']/cosmo_dic['Omega_db_0'] 
                integral_NFW = integrate.simps(y=integrand_cold, x = r_arr)
                
                
                guess = integral_NFW / integral_soliton
                
                #find the central density parameter by solving the eq: 
                #M_ax_halo = Omega_ax/Omega_cold * M_cold_halo
                def func_find_root(dens):
                    return func_ax_halo_mass(m, cosmo_dic, power_spec_dic_sigma, M_cut, dens, eta_given=eta_given) - cosmo_dic['Omega_ax_0']/cosmo_dic['Omega_db_0'] * m
                dens_param = optimize.root(func_find_root, x0 = guess).x
                
                #sometimes the solution is not really a solution,
                #so set than the central density paameter to zero, ie so solution can be found
                if np.abs(guess - dens_param) > 100:
                    dens_param_arr.append(0.)
                else:
                    dens_param_arr.append(float(dens_param))
        return dens_param_arr
            


def func_dens_profile_ax_kspace(k, M, cosmo_dic, power_spec_dic_sigma, M_cut, central_dens_param, eta_given=False):
    """
    k in units of h/Mpc and M and M_cut in solar_mass/h
    The free parameter rho_central_param is set such that
    we get the correct mass of the soliton halo, 
    see func_central_density_param
    return kspace denisty profile for the axion halo
    the normalised density profile is demensionles
    """ 
    #the kspace density profile is defined via
    # \rho(k) = 4*\pi* int_0^r_vir \rho(r) * r^2 * sin(kr)/kr dr
    M_ax = func_ax_halo_mass(M, cosmo_dic, power_spec_dic_sigma, M_cut, central_dens_param)
    r_vir = func_r_vir(M, cosmo_dic, cosmo_dic['Omega_db_0'])
    
    #distinguish whether M is an array or a scalar
    if isinstance(M, (int, float)) == True:
        r_arr = np.geomspace(1e-15, r_vir, num=1000)
        dens_profile_arr = func_dens_profile_ax(r_arr, M, cosmo_dic, power_spec_dic_sigma, M_cut, central_dens_param, eta_given=eta_given) \
                           * r_arr**2 * np.sin(np.outer(k, r_arr)) / np.outer(k, r_arr)
        return list(4 * np.pi * integrate.simps(y=dens_profile_arr, x = r_arr, axis=-1) / M_ax)
        
    else:
        dens_profile_kspace_arr = []
        for idx, m in enumerate(M):
            if  M_ax[idx] == 0:
                dens_profile_kspace_arr.append(list(np.zeros(len(k))))
            else:
                r_arr = np.geomspace(1e-15, r_vir[idx], num=1000)
                dens_profile_arr = func_dens_profile_ax(r_arr, m, cosmo_dic, power_spec_dic_sigma, M_cut, central_dens_param[idx], eta_given=eta_given) \
                                   * r_arr**2 * np.sin(np.outer(k, r_arr)) / np.outer(k, r_arr)
                dens_kspace = list(4 * np.pi * integrate.simps(y=dens_profile_arr, x = r_arr, axis=-1) / M_ax[idx] )
                dens_profile_kspace_arr.append(dens_kspace)
        return dens_profile_kspace_arr