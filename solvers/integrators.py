# Integration routines compatible with NUMBA

from numba import njit
import numpy as np

@njit
def trapz(y, x):
    '''
    Pure python version of trapezoid rule.
    Taken from https://berkeley-stat159-f17.github.io/stat159-f17/lectures/09-intro-numpy/trapezoid..html
    '''
    s = 0
    for i in range(1, len(x)):
        s += (x[i]-x[i-1])*(y[i]+y[i-1])
    return s/2

@njit
def cumtrapz(y, x):
    '''    
    Pure python version of trapezoid rule.
    Taken from https://berkeley-stat159-f17.github.io/stat159-f17/lectures/09-intro-numpy/trapezoid..html
    '''
    s = np.zeros(len(x))
    for i in range(1, len(x)):
        s[i] = s[i-1] + (x[i]-x[i-1])*(y[i]+y[i-1])
    return s/2
