"""
Load the input file into a dictionary used for the halo model and axionCAMB specifications
"""

import numpy as np
from scipy import integrate

def load_cosmology_input(params_path):
    """
    create a dictionary with all important cosmology parameter
    the input file is given in params_path and each parameter
    musst be given in a new line
    """
    cosmo_dic = {}
    file = open(params_path)
    for line in file:
        line = line.split('=')
        if line[0] == '\n':
            continue
        elif line[0].strip() == 'omega_b_0': 
            cosmo_dic['omega_b_0'] = eval(line[1].strip())
        elif line[0].strip() == 'omega_d_0': 
            total_dark = eval(line[1].strip())
        elif line[0].strip() == 'ax_fraction':
            cosmo_dic['omega_ax_0'] = total_dark * eval(line[1].strip())
            cosmo_dic['omega_d_0'] = total_dark - cosmo_dic['omega_ax_0']
            cosmo_dic['omega_db_0'] = cosmo_dic['omega_d_0'] + cosmo_dic['omega_b_0']
            cosmo_dic['omega_m_0'] = cosmo_dic['omega_db_0'] + cosmo_dic['omega_ax_0']
        elif line[0].strip() == 'm_ax':
            cosmo_dic['m_ax'] = eval(line[1].strip())
        elif line[0].strip() == 'h':
            cosmo_dic['h'] = eval(line[1].strip())
            cosmo_dic['Omega_b_0'] = cosmo_dic['omega_b_0']/cosmo_dic['h']**2
            cosmo_dic['Omega_d_0'] = cosmo_dic['omega_d_0']/cosmo_dic['h']**2
            cosmo_dic['Omega_ax_0'] = cosmo_dic['omega_ax_0']/cosmo_dic['h']**2
            cosmo_dic['Omega_db_0'] = cosmo_dic['omega_db_0']/cosmo_dic['h']**2
            cosmo_dic['Omega_m_0'] = cosmo_dic['Omega_db_0'] + cosmo_dic['Omega_ax_0']
            cosmo_dic['Omega_w_0'] = 1 - cosmo_dic['Omega_m_0']
        elif line[0].strip() == 'z':
            cosmo_dic['z'] = eval(line[1].strip())   
        elif line[0].strip() == 'M_min':
            cosmo_dic['M_min'] = eval(line[1].strip())
        elif line[0].strip() == 'M_max':
            cosmo_dic['M_max'] = eval(line[1].strip())
        elif line[0].strip() == 'ns':
            cosmo_dic['ns'] = eval(line[1].strip())
        elif line[0].strip() == 'As':
            cosmo_dic['As'] = eval(line[1].strip())  
        elif line[0].strip() == 'k_piv':
            cosmo_dic['k_piv'] = eval(line[1].strip())
        elif line[0].strip() == 'transfer_kmax':
            cosmo_dic['transfer_kmax'] = eval(line[1].strip())  
        elif line[0].split()[0] == '#':
            continue   
    file.close()
    return cosmo_dic

def load_LCDM_cosmology_input(params_path):
    """
    load the coresponding LCDM cosmology as in load_cosmology_inputs
    from the same input file in params_path
    by ignoring the axion fraction and hardcode omega_ax_0 = 1e-20.
    This is because axionCAMB needs omega_ax_0 > 0
    """
    cosmo_dic = {}
    file = open(params_path)
    for line in file:
        line = line.split('=')
        if line[0] == '\n':
            continue
        elif line[0].strip() == 'omega_b_0': 
            cosmo_dic['omega_b_0'] = eval(line[1].strip())
        elif line[0].strip() == 'omega_d_0':
            cosmo_dic['omega_d_0'] = eval(line[1].strip())
            cosmo_dic['omega_ax_0'] = 1e-20 #not equal to zero because axionCAMB needs omega_ax_0>0
            cosmo_dic['omega_db_0'] = cosmo_dic['omega_d_0'] + cosmo_dic['omega_b_0']
            cosmo_dic['omega_m_0'] = cosmo_dic['omega_db_0']
        elif line[0].strip() == 'm_ax':
            cosmo_dic['m_ax'] = eval(line[1].strip())
        elif line[0].strip() == 'h':
            cosmo_dic['h'] = eval(line[1].strip())
            cosmo_dic['Omega_b_0'] = cosmo_dic['omega_b_0']/cosmo_dic['h']**2
            cosmo_dic['Omega_d_0'] = cosmo_dic['omega_d_0']/cosmo_dic['h']**2
            cosmo_dic['Omega_db_0'] = cosmo_dic['omega_db_0']/cosmo_dic['h']**2
            cosmo_dic['Omega_m_0'] = cosmo_dic['Omega_db_0'] 
            cosmo_dic['Omega_w_0'] = 1 - cosmo_dic['Omega_m_0']
        elif line[0].strip() == 'z':
            cosmo_dic['z'] = eval(line[1].strip())   
        elif line[0].strip() == 'M_min':
            cosmo_dic['M_min'] = eval(line[1].strip())
        elif line[0].strip() == 'M_max':
            cosmo_dic['M_max'] = eval(line[1].strip())
        elif line[0].strip() == 'ns':
            cosmo_dic['ns'] = eval(line[1].strip())
        elif line[0].strip() == 'As':
            cosmo_dic['As'] = eval(line[1].strip())  
        elif line[0].strip() == 'k_piv':
            cosmo_dic['k_piv'] = eval(line[1].strip())
        elif line[0].strip() == 'transfer_kmax':
            cosmo_dic['transfer_kmax'] = eval(line[1].strip())  
        elif line[0].split()[0] == '#':
            continue   
    file.close()
    return cosmo_dic