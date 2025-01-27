# Before running, activate the virtual environment:
# source /mnt/users/crookmansourj/swiftenv/bin/activate
# This environment has swiftools 3.0.22

from get_xrt_from_pipeline import get_xrt_prods, group_spectra
from fit_xrt_spectra import run_spectral_fit
from get_results import plot_lightcurve_and_hr, get_results, plot_all_spectral_fit_results, plot_spectral_results, get_index_range
import subprocess
import glob
import os
import sys
sys.path.append(os.path.abspath("./uplims_analysis"))
from calculate_uplims import uplims_runner
import numpy as np
from input_parameters import analysis_type,target_coords,target_names,target_ids,segments_lc,segments_spec,MJD_range,models,nH,gamma,low_count_indexes,fix,uplims_IDs,uplims_count_rate_ar,uplims_MJDs,uplims_MJDs_er,uplims_fluxes,models_indexes,low_energy_fit_bound_keV

"""
Order in which to run:
(1) get_data
(2) plot_lc. This will output the upper limit count rate and MJDs (if any) that will need to be dealt with separately.
(3) fit_spec. If you have an nH value, you can set this to be fixed.
(4) Look at the spectral_fit_results and spectral_fit_residuals folders.
(5) If the fit is failing in some epochs because of too low counts, you may need to determine parameters to fix.
(6) You can run get_parameter_results with 
(7) If some of the fit parameters (e.g. nH) need to be changed, re-run fit_spec with the value of nH. 
(8) If upper limits were identified from the plot_lc routine, run get_uplim_fluxes (with the output nH from the previous step, if needs be).
(9) If needs be, run get_index_range, to get the range of indexes corresponding to a particular MJD range.
(10) Run get_spec_results (with the output of get_uplim_fluxes).
"""

# TODO: 
# scrape of observation IDs from the internet
# scrape target coordinates off the internet
# recommended way to make requests, to avoid overload
# Possibly use a yaml file for inputs?
# Get IDs for the upper limits automatically, using website?
# Make it such that the fit_spec function fetches the upper limit IDs, and so does not require this user input.
# Minimise the need for user input by calling from other functions to fetch data e.g. uplims
if __name__ in "__main__":

    print("Analysis type: ",  analysis_type)

    
    ## STEP 1
    ## GET LIGHTCURVES & SPECTRA
    if analysis_type=="step1": 

        # Iterate through Target name(s), id(s), and segments and get xrt_prods
        
        for target_name, target_id, segments_lc in zip(target_names, target_ids, segments_lc):
            print("TARGET: ", target_name)
            print("Getting light curve:")
            get_xrt_prods(target_id, target_name, target_coords, segments_lc, centroid=True, prod_type = 'lightcurve')
            print("-------------------------------------------------------------------------------------------------------")
        print("-------------------------------------------------------------------------------------------------------")
        
        for target_name, target_id, segments_spec in zip(target_names, target_ids, segments_spec):
            print("TARGET: ", target_name)
            try:
                print("Getting spectrum:")
                get_xrt_prods(target_id, target_name, target_coords, segments_spec, centroid=True, prod_type = 'spectrum')
            except:
                print("Error getting spectrum.")
            print("-------------------------------------------------------------------------------------------------------")
        
        # Group the spectra -- required before running spectral fits
        group_spectra()
        
    
    ## STEP 2 
    ## PLOT THE LIGHT CURVE AND HARDNESS RATIO
    elif analysis_type=="step2": 
        plot_lightcurve_and_hr() 

    ## STEP 3
    ## FIT THE SPECTRA AND GET INITIAL SPECTRAL FIT RESULTS
    elif analysis_type=="step3": 
        run_spectral_fit(models, uplims_IDs, low_energy_fit_bound_keV, fix={}) # fit the spectra
        get_results(models, fixing=False) # output the fit results
        plot_all_spectral_fit_results(models, fixing=False) # plot all the fit results

    ## STEP 4
    ## GET THE PARAMETER VALUES FROM THE SPECTRAL RESULTS, BASED ON CHOSEN MODELS
    elif analysis_type=="step4": 
        get_results(models, models_indexes, fixing=False)

    ## STEP 5
    ## FIT THE SPECTRA AND GET SPECTRAL FIT RESULTS -- WITH PARAMETERS TO FIX
    elif analysis_type=="step5": 
        run_spectral_fit(models, uplims_IDs, low_energy_fit_bound_keV, fix=fix) # fit the spectra
        get_results(models, fixing=True) # output the fit results

    ## STEP 6
    ## GET THE UPPER LIMIT FLUX VALUES
    elif analysis_type=="step6": 
        uplims_runner(uplims_IDs, uplims_count_rate_ar, target_coords, nH, gamma)

    ## STEP 7
    ## GET THE FINAL RESULTS, BASED ON CHOSEN MODELS
    elif analysis_type=="step7": 
        plot_spectral_results(models, models_indexes, uplims_IDs, uplims_MJDs, uplims_MJDs_er, uplims_fluxes, fixing=True)


    ## HELPER FUNCTIONS 
    ## GET THE INDEX RANGE
    # Find the index range in dates_MJD for a particular MJD_range, where MJD_range = [start_MJD, end_MJD)
    # The output index range is inclusive, i.e. [start_index, end_index]
    elif analysis_type=="get_index_range":
        get_index_range(MJD_range) 

    else:
        print("Analysis_type not valid.")
