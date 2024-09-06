"""
Wrapper for calling axionCAMB from inside Python programs.
"""

#import packages
import os
import re
import subprocess
import sys
import fileinput
import numpy as np

#create param infile for axionCAMB
def axioncamb_params(params_path, arg_dic={}, print_info = True,
                ombh2=0.0224, omch2=0.12, omnuh2=0, omk=0, hubble=67.4,  
                omaxh2=0.001, m_ax=1.e-26, 
                scalar_amp__1___=2.1e-9, pivot_scalar=0.05, 
                scalar_spectral_index__1___=0.965,
                transfer_kmax=1000,
                transfer_redshift__1___=0,
                output_root='test', 
                get_scalar_cls='F',
                get_vector_cls='F', get_tensor_cls='F', get_transfer='T',
                do_lensing='F', do_nonlinear=0,
                l_max_scalar=4000, k_eta_max_scalar=5000,
                l_max_tensor=1500, k_eta_max_tensor=3000, use_physical='T',
                axion_isocurvature='F', alpha_ax=0, Hinf=13.7,
                w=-1, use_axfrac='F',
                cs2_lam=1, temp_cmb=2.726, helium_fraction=0.24,
                massless_neutrinos=3.046, massive_neutrinos=0,
                nu_mass_eigenstates=1, nu_mass_fractions=1, share_delta_neff='T',
                initial_power_num=1,
                pivot_tensor=0.05, 
                scalar_nrun__1___=0,
                tensor_spectral_index__1___=0,
                initial_ratio__1___=0, tens_ratio=0,
                reionization='T', re_use_optical_depth='T',
                re_optical_depth=0.06, # to be consistent when combining with CMB surveys
                re_redshift=11, re_delta_redshift=1.5, re_ionization_frac=-1,
                RECFAST_fudge=1.14, RECFAST_fudge_He=0.86, RECFAST_Heswitch=6,
                RECFAST_Hswitch='T', initial_condition=1,
                initial_vector='-1 0 0 0 0', vector_mode=0, COBE_normalize='F',
                CMB_outputscale=7.4311e12, transfer_high_precision='T', transfer_k_per_logint=0,
                transfer_num_redshifts=1, transfer_interp_matterpower='T', transfer_filename__1___='transfer_out.dat',
                transfer_matterpower__1___='matterpower.dat',
                scalar_output_file='scalCls.dat', vector_output_file='vecCls.dat',
                tensor_output_file='tensCls.dat', total_output_file='totCls.dat',
                lensed_output_file='lensedCls.dat',
                lensed_total_output_file='lensedtotCls.dat',
                lens_potential_output_file='lenspotentialCls.dat',
                FITS_filename='scalCls.fits', do_lensing_bispectrum='F',
                do_primordial_bispectrum='F', bispectrum_nfields=1,
                bispectrum_slice_base_L=0, bispectrum_ndelta=3,
                bispectrum_delta__1___=0, bispectrum_delta__2___=2,
                bispectrum_delta__3___=4, bispectrum_do_fisher='F',
                bispectrum_fisher_noise=0, bispectrum_fisher_noise_pol=0,
                bispectrum_fisher_fwhm_arcmin=7,
                bispectrum_full_output_file='',
                bispectrum_full_output_sparse='F',
                feedback_level=1, lensing_method=1, accurate_BB='F', massive_nu_approx=1,
                accurate_polarization='T', accurate_reionization='T',
                do_tensor_neutrinos='T', do_late_rad_truncation='T',
                number_of_threads=0, high_accuracy_default='F',
                 accuracy_boost=1, l_accuracy_boost=1, l_sample_boost=1):
    """
    Define a file with all parameters for running axionCAMB.
    Change some of them (if given) by the values in a given dictionary arg_dic
    file is saved in given path: params_path
    """
    # Get dict. of arguments 
    args = locals()

    # Get all parameters into the CAMB param.ini format
    camb_params_text = ""
    for key in args:
        keyname = key
        if "__" in key:  # Rename array parameters
            keyname = key.replace("___", ")").replace("__", "(")
        line_str = "=".join((keyname, str(args[key])))
        camb_params_text += line_str + "\n"

    # Output params file
    if print_info == True:
        print("Writing axionCAMB parameters to", params_path)
    f = open(params_path, 'w')
    f.write(camb_params_text)
    f.close()
    
    #change parameters by the given arg_dic
    if bool(arg_dic) == True:
        for line in fileinput.FileInput(params_path, inplace=1):
            sline=line.strip().split("=")
            if sline[0].startswith("ombh2"):
                sline[1]=str(arg_dic['omega_b_0']) #baryons
                line='='.join(sline) 
            elif sline[0].startswith("omch2"):
                sline[1]=str(arg_dic['omega_d_0']) #cold dark
                line='='.join(sline)
            elif sline[0].startswith("omaxh2"):
                sline[1]=str(arg_dic['omega_ax_0']) #axions
                line='='.join(sline)
            elif sline[0].startswith("m_ax"):
                sline[1]=str( arg_dic['m_ax']) #mass axion
                line='='.join(sline)
            elif sline[0].startswith("hubble"):
                sline[1]=str( arg_dic['h'] * 100.) #Hubble
                line='='.join(sline)
            elif sline[0].startswith("transfer_redshift(1)"):
                sline[1]=str( arg_dic['z']) #redshift
                line='='.join(sline)
            elif sline[0].startswith("scalar_spectral_index(1)"):
                sline[1]=str( arg_dic['ns']) #pestral index
                line='='.join(sline)
            elif sline[0].startswith("scalar_amp(1)"):
                sline[1]=str( arg_dic['As']) #scalar ampitude
                line='='.join(sline)
            elif sline[0].startswith("pivot_scalar"):
                sline[1]=str( arg_dic['k_piv']) #pivon scale in 1/kpc
                line='='.join(sline)
            elif sline[0].startswith("transfer_kmax"):
                sline[1]=str(arg_dic['transfer_kmax']) #max k to evaluate
                line='='.join(sline)
            else:
                line='='.join(sline) 
            print(line)

    

def run_axioncamb(params_path, camb_exec_dir, arg_dic, print_info = True):
    """
    Run axionCAMB, using a given (pre-written) params file (see axioncamb_params)
    stored in the file name params_path.
    The camb_exec_dir is the absolute path to the directory of your axionCAMB executable.
    The arg_dic musst be the dictionary retured from the function axioncamb_params
    All 
    """
    
    cwd = os.getcwd()
    file_path = cwd + "/" + params_path #file path for axionCAMB
    if print_info == True:
        print("Running CAMB on", file_path)
    output = subprocess.check_output([camb_exec_dir + "/" + "./camb", file_path]).decode(sys.stdout.encoding) #run axionCAMB

    #go through the output to save sigma8 of the PS in the parameter dictionary arg_dic
    if bool(np.array(arg_dic['z']).shape): # check if array (skipped if float)
        arg_dic['sigma_8'] = []
        
        for line in output.split("\n"):
            if "sigma8" in line:
                s8line = line[line.find('sigma8'):]  # Only get sigma8 part of string
                sigma8 = float(re.findall(r'\b\d+.+\b', s8line)[0])
                arg_dic['sigma_8'].append(sigma8) # save sigma8 in arg_dic
                
    else:
        for line in output.split("\n"):
            if "sigma8" in line:
                s8line = line[line.find('sigma8'):]  # Only get sigma8 part of string
                sigma8 = float(re.findall(r'\b\d+.+\b', s8line)[0])
                arg_dic['sigma_8'] = sigma8 # save sigma8 in arg_dic

    
def redshift_array(params_path, arg_dic, print_info=False):
    """
    Changes input file to calculate multiple redshifts with one axionCAMB call
    
    Usage:
    (1) Run axioncamb params function (assuming single redshift)
    (2) Overwrite arg_dic['z'] with z array
    (3) Run this function
    (4) Run run_axioncamb with updated arg_dic and input file
    """
    
    redshift_array = np.sort(arg_dic['z'])[::-1] # start with earlier redshift first
    with open(params_path, 'r') as input_file:
        lines = input_file.readlines()
    
    with open(params_path, 'w') as input_file:
        # remove all instances of redshift-specific out files 
        for line in lines:
            if "transfer_redshift" in line or "transfer_file" in line or "transfer_matter" in line:
                pass
            elif "transfer_num_redshifts" in line:
                sline=line.strip().split("=")
                sline[1]=str( len(arg_dic['z']))
                line='='.join(sline)
                input_file.write(line + '\n' )
            else:
                input_file.write(line)
            
        # refill fields with redshift array data
        for i in range(1, len(redshift_array)+1):
            input_file.write('\n' + 'transfer_redshift('+ str(i) + ')=' + str(redshift_array[i-1]))
            input_file.write('\n' + 'transfer_filename__' + str(i) + '__=transfer_out_' + str(i)  + '.dat')
            input_file.write('\n' + 'transfer_matterpower__' + str(i) + '__=matterpower_' + str(i) +'.dat')

    if print_info:
        with open(params_path, 'r') as input_file:
            lines = input_file.readlines()
            for line in lines:
                print(line)

    return
