"""
interpolation functions for the power spectrum
"""

#packages needed
import numpy as np
from scipy import interpolate, misc, optimize

#own packages 
from cosmology.basic_cosmology import *

def log_interp1d(x,y,\
                 axis=-1,\
                 kind='linear',\
                 copy=True,\
                 bounds_error=None,\
                 fill_value=np.nan,\
                 assume_sorted=False):
    """
    log interpolater, code from Meads former python HMcode
    """
    logx = np.log(x)
    logy = np.log(y)
    lin_interp = interpolate.interp1d(logx,logy,\
                                    kind=kind,\
                                    axis=axis,\
                                    copy=copy,\
                                    bounds_error=bounds_error,\
                                    fill_value=fill_value,\
                                    assume_sorted=assume_sorted)
    log_interp = lambda z: np.exp(lin_interp(np.log(z)))
    return log_interp


def func_PS_interpolate(k_arr_value, k_arr, PS_arr):
    """
    k_arr, k_arr,value units of h/Mpc, M in solar_mass/h and PK in (Mpc/h)^3
    NOTE: we inerpolate the power spectrum PS_arr(k_arr)
    at the point in k_arr_value
    code from Meads former python HMcode
    """
    # Create P(k) interpolation function
    PS_func = log_interp1d(k_arr, PS_arr, kind='cubic')
    def func_PS_find(k_value):
        if(k_value<k_arr[0]):
            a=np.log(PS_arr[1]/PS_arr[0])/np.log(k_arr[1]/k_arr[0])
            b=np.log(PS_arr[0])-a*np.log(k_arr[0])
            return np.exp(a*np.log(k_value)+b)
        elif(k_value>k_arr[-1]):
            a=np.log(PS_arr[-2]/PS_arr[-1])/np.log(k_arr[-2]/k_arr[-1])
            b=np.log(PS_arr[-2])-a*np.log(k_arr[-2])
            return np.exp(a*np.log(k_value)+b)
        else:
            return PS_func(k_value)
        
    PS_interpolate = [func_PS_find(k_value) for k_value in k_arr_value]
    #for k_value in k_arr_value:
    #    PS_interpolate.append(func_PS_find(k_value))
    return np.array(PS_interpolate)

def func_PS_interpolate_M(M, k, PS, cosmo_dic, Omega_0):
    """
    k is in units of h/Mpc, PS in (Mpc/h)^3 and M in solar_mass/h
    returns interpolated PS and k. The interpolation range
    depends on M, because the interpolated PS is needed to compute
    the variance over a wide range of scales/masses correct
    NOTE_ PS and Omega_0 musst match
    """
    #k_max is chosen such that it is 100 times larger than the k that 
    #coresponds to the smallest mass: k = 2*\pi/R(M) 
    k_max = 100*2 * np.pi / func_R_M(np.min(M), cosmo_dic, Omega_0)

    if k_max < np.max(k):
        return k, PS
    else:
        k_interp = np.logspace(np.log10(np.min(k)), np.log10(k_max), 500)
        PS_interp = func_PS_interpolate(k_interp, k, PS)
        return k_interp, PS_interp



