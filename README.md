# axionHMcode

`axionHMcode` is a code to compute the non-linear matter power spectrum in a mixed dark matter cosmology with an ultra-light axion (ULA) component of the dark matter as described in [Vogt et al.](https://arxiv.org/abs/2209.13445). A very accurate halo model for a $\Lambda$CDM or massive neutrino cosmology is given in [Mead et al.](https://arxiv.org/abs/2009.01858) by the `HMCode-2020`. Since this model uses some of the fitting parameters and is inspired by this code the presented code is named after `HMcode-2020`. The code was updated and calibrated to simulations in [Dome et al.](https://arxiv.org/abs/2409.11469)

## Theory

The model computes the non-linear power spectrum by using the fully expanded power spectrum

![codesketch](eq_halo_model.png)


The cold part can be computed as usually with the standard halo model (see [Massara et al.](https://arxiv.org/abs/1410.6813) or [Mead et al.](https://arxiv.org/abs/2009.01858)). In contrast, the cross and axion parts have to take into account the non clustering of axions on small scales due to free-streaming. This is done by splitting the axion overdensity into a clustered and linear component. For details see [Massara et al.](https://arxiv.org/abs/1410.6813) where the same full treatment was used for massive neutrinos, but can be translated to any other warm/hot, i.e. free streaming / (partially) non-clustering, matter component. 

## How Does the Code Work?

The code expects an input file as given in `input_file.txt` which contains the information about the cosmology and save it in a pyhton dictionary. Alternatively, one can also give the corresponding dictionary directly. Have a look into `load_cosmology.py` to see what informations has to be stored in this dictionary. Furthermore, the linear matter power spectrum for the given cosmology for the different components, axions, CDM and baryons saved in a dictionary is needed (see example notebook for the structure). This code includes a wrapper for `axionCAMB`. So if you have a working `axionCAMB` executable, it can directly compute the linear power spectrum from the cosmology crated from the `input_file.txt`.

An example python file is given in `example_file.py`. To run the file you have to change the ‘input_file_path’ and the ‘axionCAMB_exe_path’ (complete path). If the paths are not correct the python code will produce an error message. Besides the non-linear total matter power spectrum, the example file also computes the non-linear power spectrum in a $\Lambda$CDM cosmology where the axion density is transformed into CDM density. Both power spectra are saved in a file named by the variable `datafile_path`. The units of the wavenumber and the power spectra are h/Mpc and (Mpc/h)^3 respectively. The code also produces a plot of the ratio between the MDM and $\Lambda$CDM linear and non-linear power spectra.

Moreover, for consistency the example notebook makes a comparison to `HMcode-2020` in the $\Lambda$CDM case. A comparison is made for the implementation in `CAMB` and the python package `hmcode`.


## HMCode-2020 Parameters

The `axionHMcode` can also use the parameters from the `HMCode-2020` in [Mead et al.](https://arxiv.org/abs/2009.01858) which improves the predictions in the case of a LCDM cosmology with massive neutrinos. The parameters can be switched on by setting the corresponding parameters to `True` or `False` in the function for the non-linear power spectrum (make sure the parameters are consistent in the different functions). The parameters are the smoothing parameters, `alpha`, the halo bloating term, `eta_given`, the one halo damping on large scales, `one_halo_damping`, and the two halo damping on large scales, `two_halo_damping`.

Because these parameters are not calibrated to mixed axion - cold dark matter cosmologies we are suggest to not use the parameters in the base version.


## Code Updates including calibration from [Dome et al.](https://arxiv.org/abs/2409.11469)

[Dome et al.](https://arxiv.org/abs/2409.11469) presented an improved version of `axionHMcode`, which introduced new parameters which are calibrated from simulations presented in this paper.

Moreover, several updated to speed up the code were made

The changes and bug fixed in the code are the following:
1. the axion-mass cold mass relation, $M_a(M_c)$, is now modeled as a broken power law (see Eq. 51 in [Dome et al.](https://arxiv.org/abs/2409.11469)) which ensures the old relation of $M_a(M_c) = \Omega_a/\Omega_c M_c$ above a defined cut-off mass. This new relation is inspired by simulations and is calibrated by them in the redshift range $1 < z < 8$ and for axion fractions of $0.01 < f_{\mathrm{ax}} = \Omega_{\mathrm{ax}} / \Omega_{\mathrm{m}} < 0.3$.
2. [Dome et al.](https://arxiv.org/abs/2409.11469) introduced new smoothing parameters $\alpha$ for the cold-cold power spectrum and the cross-power spectrum which depend on the axion mass and axion density. The exact form of these parameters was calibrated by simulations n the redshift range $1 < z < 8$ and for axion fractions of $0.01 < f_{\mathrm{ax}} = \Omega_{\mathrm{ax}} / \Omega_{\mathrm{m}} < 0.3$.
3. The two halo term is now calculated by the linear power spectrum only. The difference between the total two halo term is minor and the speed up is increased. There is still the option to use the full two halo term by setting the parameter `full_2h = True`. 
4. A bug in the cold density profile when using the halo bloating term was corrected.
5. The relations of the critical density threshold, $\delta_c$, and the virial overdensity, $\Delta_{\mathrm{vir}}$, are now calculated as in `HMCode-2020` (se  [Mead et al.](https://arxiv.org/abs/2009.01858) Eq. A1 and A2). This ensures that `axionHMcode` and `HMCode-2020` agree in the case of a $\Lambda$CDM cosmology.
6. The minimum concentration is now 5.196 as found in [Mead et al.](https://arxiv.org/abs/2009.01858).
7. Alex Laguë implemented Numba in this updated version to increase the speed of the code
8. Alex Laguë and Keir Rogers also included the optional parameters alpha_1, alpha_2, gamma_1, gamma_2 defined in [Dentler er al.](https://arxiv.org/abs/2111.01199) Eq. (36). To use them, just include them in your dictionary of cosmological parameters (the `cosmo_dic` file) before running the `params` and `power spectra` calculation using e.g. `cosmo_dic['alpha_1'] = X`. They are not yet included in the input file.

It is nor possible to use the base version as well as the calibrated version from [Dome et al.](https://arxiv.org/abs/2409.11469). For this one as to set the parameter `version` in the input file to `basic` or `dome`, respectively. 

When using the calibrated version from [Dome et al.](https://arxiv.org/abs/2409.11469) one has to take into account the redshift range as well as the axion fraction on which the model was calibrated. Therefore, we only recomment to use the `dome` version if $1 < z < 8$ and  $0.01 < f_{\mathrm{ax}} = \Omega_{\mathrm{ax}} / \Omega_{\mathrm{m}} < 0.3$. In addition in the `dome` version, we recoment to use `alpha = True` and `concentration_param = True`.

## Contact data

If you find any bugs or have any questions with respect to the code, please send me a message via github or open an issue.