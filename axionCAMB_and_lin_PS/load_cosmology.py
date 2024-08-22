"""
Load the input file into a dictionary used for the halo model and axionCAMB specifications
"""

import numpy as np
from scipy import integrate

# def load_cosmology_input(params_path):
#     """
#     create a dictionary with all important cosmology parameter
#     the input file is given in params_path and each parameter
#     musst be given in a new line
#     """
#     cosmo_dic = {}
#     # default values
#     cosmo_dic['M_min'] = 1e8
#     cosmo_dic['M_max'] = 1e17
#     cosmo_dic['transfer_kmax'] = 1e3
#     file = open(params_path)
#     for line in file:
#         line = line.split('=')
#         if line[0] == '\n':
#             continue
#         elif line[0].strip() == 'omega_b_0' or line[0].strip() == 'Omega_b_0': 
#             if line[0].strip() == 'omega_b_0':
#                 cosmo_dic['omega_b_0'] = eval(line[1].strip())
#             else:
#                 cosmo_dic['Omega_b_0'] = eval(line[1].strip())
#         elif line[0].strip() == 'omega_d_0' or line[0].strip() == 'Omega_d_0': 
#             if line[0].strip() == 'omega_d_0':
#                 total_dark = eval(line[1].strip()) # when reduced density param is given: Omega_i*h^2
#             else:
#                 Total_Dark = eval(line[1].strip()) # when denisty parame Omega is given
            
#         elif line[0].strip() == 'ax_fraction':
#             cosmo_dic['omega_ax_0'] = total_dark * eval(line[1].strip())
#             cosmo_dic['omega_d_0'] = total_dark - cosmo_dic['omega_ax_0']
#             cosmo_dic['omega_db_0'] = cosmo_dic['omega_d_0'] + cosmo_dic['omega_b_0']
#             cosmo_dic['omega_m_0'] = cosmo_dic['omega_db_0'] + cosmo_dic['omega_ax_0']
#         elif line[0].strip() == 'm_ax':
#             cosmo_dic['m_ax'] = eval(line[1].strip())
#         elif line[0].strip() == 'h':
#             cosmo_dic['h'] = eval(line[1].strip())
#             cosmo_dic['Omega_b_0'] = cosmo_dic['omega_b_0']/cosmo_dic['h']**2
#             cosmo_dic['Omega_d_0'] = cosmo_dic['omega_d_0']/cosmo_dic['h']**2
#             cosmo_dic['Omega_ax_0'] = cosmo_dic['omega_ax_0']/cosmo_dic['h']**2
#             cosmo_dic['Omega_db_0'] = cosmo_dic['omega_db_0']/cosmo_dic['h']**2
#             cosmo_dic['Omega_m_0'] = cosmo_dic['Omega_db_0'] + cosmo_dic['Omega_ax_0']
#             cosmo_dic['Omega_w_0'] = 1 - cosmo_dic['Omega_m_0']
#         elif line[0].strip() == 'z':
#             cosmo_dic['z'] = eval(line[1].strip())   
#         elif line[0].strip() == 'ns':
#             cosmo_dic['ns'] = eval(line[1].strip())
#         elif line[0].strip() == 'As':
#             cosmo_dic['As'] = eval(line[1].strip())  
#         elif line[0].strip() == 'k_piv':
#             cosmo_dic['k_piv'] = eval(line[1].strip())

#         elif line[0].strip() == 'transfer_kmax':
#             cosmo_dic['transfer_kmax'] = eval(line[1].strip()) 
#         elif line[0].strip() == 'M_min':
#             cosmo_dic['M_min'] = eval(line[1].strip())
#         elif line[0].strip() == 'M_max':
#             cosmo_dic['M_max'] = eval(line[1].strip()) 
#         elif line[0].split()[0] == '#':
#             continue   
#     file.close()
#     return cosmo_dic


def load_cosmology_input(params_path):
    """
    create a dictionary with all important cosmology parameter
    the input file is given in params_path and each parameter
    must be given in a new line
    """
    cosmo_dic = {}
    # default values
    cosmo_dic['M_min'] = 1e8
    cosmo_dic['M_max'] = 1e17
    cosmo_dic['transfer_kmax'] = 1e3

    total_dark = None
    Total_dark = None
    num_om = 0

    file = open(params_path)
    for line in file:
        line = line.split('=')
        if line[0] == '\n' or line[0].strip().startswith('#'):
            continue
        
        key = line[0].strip()
        value = eval(line[1].strip())

        if key == 'omega_b_0':
            cosmo_dic['omega_b_0'] = value
            num_om = num_om + 1
        elif key == 'Omega_b_0':
            cosmo_dic['Omega_b_0'] = value
            num_om = num_om + 1
        elif key == 'omega_d_0':
            num_om = num_om + 1
            total_dark = value
        elif key == 'Omega_d_0':
            num_om = num_om + 1
            Total_dark = value
        elif key == 'omega_m_0':
            num_om = num_om + 1
            cosmo_dic['omega_m_0'] = value
        elif key == 'Omega_m_0':
            num_om = num_om + 1
            cosmo_dic['Omega_m_0'] = value
        elif key == 'ax_fraction':
            ax_fraction = value
        elif key == 'm_ax':
            cosmo_dic['m_ax'] = value
        elif key == 'h':
            cosmo_dic['h'] = value
        elif key == 'z':
            cosmo_dic['z'] = value   
        elif key == 'ns':
            cosmo_dic['ns'] = value
        elif key == 'As':
            cosmo_dic['As'] = value  
        elif key == 'k_piv':
            cosmo_dic['k_piv'] = value
        elif key == 'transfer_kmax':
            cosmo_dic['transfer_kmax'] = value 
        elif key == 'M_min':
            cosmo_dic['M_min'] = value
        elif key == 'M_max':
            cosmo_dic['M_max'] = value 
        
    file.close()

    if num_om != 2:
        raise ValueError('Provide exactly two of the following parameters: omega_b_0, Omega_b_0, omega_d_0, Omega_d_0, omega_m_0, Omega_m_0')

    # when total matter Omega_m and baryons Omega_b are given, calculate the rest
    if 'Omega_m_0' in cosmo_dic and 'Omega_b_0' in cosmo_dic:
        cosmo_dic['omega_m_0'] = cosmo_dic['Omega_m_0'] * cosmo_dic['h']**2
        cosmo_dic['omega_b_0'] = cosmo_dic['Omega_b_0'] * cosmo_dic['h']**2
        Total_dark = cosmo_dic['Omega_m_0'] - cosmo_dic['Omega_b_0']
        cosmo_dic['omega_ax_0'] = Total_dark * ax_fraction
        cosmo_dic['Omega_ax_0'] = cosmo_dic['omega_ax_0'] / cosmo_dic['h']**2
        cosmo_dic['omega_d_0'] = total_dark - cosmo_dic['omega_ax_0']
        cosmo_dic['Omega_d_0'] = cosmo_dic['omega_d_0'] / cosmo_dic['h']**2
        cosmo_dic['omega_db_0'] = cosmo_dic['omega_d_0'] + cosmo_dic['omega_b_0']
        cosmo_dic['Omega_db_0'] = cosmo_dic['omega_db_0'] / cosmo_dic['h']**2
    elif 'omega_m_0' in cosmo_dic and 'omega_b_0' in cosmo_dic:
        cosmo_dic['Omega_m_0'] = cosmo_dic['omega_m_0'] / cosmo_dic['h']**2
        cosmo_dic['Omega_b_0'] = cosmo_dic['omega_b_0'] / cosmo_dic['h']**2
        total_dark = cosmo_dic['Omega_m_0'] - cosmo_dic['Omega_b_0']
        cosmo_dic['omega_ax_0'] = total_dark * ax_fraction
        cosmo_dic['Omega_ax_0'] = cosmo_dic['omega_ax_0'] / cosmo_dic['h']**2
        cosmo_dic['omega_d_0'] = total_dark - cosmo_dic['omega_ax_0']
        cosmo_dic['Omega_d_0'] = cosmo_dic['omega_d_0'] / cosmo_dic['h']**2
        cosmo_dic['omega_db_0'] = cosmo_dic['omega_d_0'] + cosmo_dic['omega_b_0']
        cosmo_dic['Omega_db_0'] = cosmo_dic['omega_db_0'] / cosmo_dic['h']**2

     # when total matter Omega_m and total dark matter Omega_d (contains CDM AND axions!) are given, calculate the rest
    if 'Omega_m_0' in cosmo_dic  and Total_dark is not None:
        cosmo_dic['omega_m_0'] = cosmo_dic['Omega_m_0'] * cosmo_dic['h']**2
        cosmo_dic['omega_d_0'] = cosmo_dic['Omega_d_0'] * cosmo_dic['h']**2

        
        cosmo_dic['Omega_ax_0'] = Total_dark * ax_fraction
        cosmo_dic['omega_ax_0'] = cosmo_dic['Omega_ax_0'] * cosmo_dic['h']**2


        cosmo_dic['Omega_d_0'] = Total_dark - cosmo_dic['Omega_ax_0']
        cosmo_dic['omega_d_0'] = cosmo_dic['Omega_d_0'] * cosmo_dic['h']**2

        cosmo_dic['Omega_b_0'] = cosmo_dic['Omega_m_0'] - Total_dark 
        cosmo_dic['omega_b_0'] = cosmo_dic['Omega_b_0'] * cosmo_dic['h']**2
        cosmo_dic['omega_db_0'] = cosmo_dic['omega_d_0'] + cosmo_dic['omega_b_0']
        cosmo_dic['Omega_db_0'] = cosmo_dic['omega_db_0'] / cosmo_dic['h']**2
    elif 'omega_m_0' in cosmo_dic and total_dark is not None:
        cosmo_dic['Omega_m_0'] = cosmo_dic['omega_m_0'] / cosmo_dic['h']**2
        
        cosmo_dic['omega_ax_0'] = total_dark * ax_fraction
        cosmo_dic['Omega_ax_0'] = cosmo_dic['omega_ax_0'] / cosmo_dic['h']**2


        cosmo_dic['omega_d_0'] = total_dark - cosmo_dic['omega_ax_0']
        cosmo_dic['Omega_d_0'] = cosmo_dic['omega_d_0'] / cosmo_dic['h']**2

        cosmo_dic['omega_b_0'] = cosmo_dic['Omega_m_0'] - total_dark 
        cosmo_dic['Omega_b_0'] = cosmo_dic['omega_b_0'] / cosmo_dic['h']**2
        cosmo_dic['omega_db_0'] = cosmo_dic['omega_d_0'] + cosmo_dic['omega_b_0']
        cosmo_dic['Omega_db_0'] = cosmo_dic['omega_db_0'] / cosmo_dic['h']**2


    # when total dark Omega_d (contains CDM AND axions!) and baryon Omega_b are given, calculate the rest
    if 'Omega_b_0' in cosmo_dic and Total_dark is not None:
        cosmo_dic['omega_b_0'] = cosmo_dic['Omega_b_0'] * cosmo_dic['h']**2
        
        cosmo_dic['Omega_ax_0'] = Total_dark * ax_fraction
        cosmo_dic['omega_ax_0'] = cosmo_dic['Omega_ax_0'] * cosmo_dic['h']**2

        cosmo_dic['Omega_d_0'] = Total_dark - cosmo_dic['Omega_ax_0']
        cosmo_dic['omega_d_0'] = cosmo_dic['Omega_d_0'] * cosmo_dic['h']**2

        cosmo_dic['Omega_m_0'] = cosmo_dic['Omega_b_0'] + total_dark 
        cosmo_dic['omega_m_0'] = cosmo_dic['Omega_m_0'] * cosmo_dic['h']**2

        cosmo_dic['omega_db_0'] = cosmo_dic['omega_d_0'] + cosmo_dic['omega_b_0']
        cosmo_dic['Omega_db_0'] = cosmo_dic['omega_db_0'] / cosmo_dic['h']**2

    elif 'omega_b_0' in cosmo_dic and total_dark is not None:
        cosmo_dic['Omega_b_0'] = cosmo_dic['omega_b_0'] / cosmo_dic['h']**2
        
        cosmo_dic['omega_ax_0'] = total_dark * ax_fraction
        cosmo_dic['Omega_ax_0'] = cosmo_dic['omega_ax_0'] / cosmo_dic['h']**2

        cosmo_dic['omega_d_0'] = total_dark - cosmo_dic['omega_ax_0']
        cosmo_dic['Omega_d_0'] = cosmo_dic['omega_d_0'] / cosmo_dic['h']**2

        cosmo_dic['omega_m_0'] = cosmo_dic['omega_b_0'] + total_dark 
        cosmo_dic['Omega_m_0'] = cosmo_dic['omega_m_0'] / cosmo_dic['h']**2

        cosmo_dic['omega_db_0'] = cosmo_dic['omega_d_0'] + cosmo_dic['omega_b_0']
        cosmo_dic['Omega_db_0'] = cosmo_dic['omega_db_0'] / cosmo_dic['h']**2

    cosmo_dic['Omega_w_0'] = 1 - cosmo_dic['Omega_m_0']

    
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