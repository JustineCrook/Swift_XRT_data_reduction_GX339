# Before running, activate the virtual environment:
# source /mnt/users/crookmansourj/swiftenv/bin/activate
# This environment has swiftools 3.0.22

from get_xrt_from_pipeline import get_xrt_prods, group_spectra
from fit_xrt_spectra import run_spectral_fit
from get_results import plot_lightcurve_and_hr, get_results, plot_all_spectral_fit_results, plot_spectral_results, get_index_range, f_test
import subprocess
import glob
import os
import sys
sys.path.append(os.path.abspath("./uplims_analysis"))
from calculate_uplims import uplims_runner
import numpy as np
from input_parameters import analysis_type,target_coords,target_names,target_ids,segments_lc,segments_spec,MJD_range,models_unconstrained,models_constrained,nH,gamma,low_count_indexes,fix,uplims_IDs,uplims_count_rate_ar,uplims_MJDs,uplims_MJDs_er,uplims_fluxes,models_indexes,min_E_keV,min_counts_chi, simple_model, complex_model, models_indexes_av



if __name__ in "__main__":

    print("Analysis type: ",  analysis_type)


    ## STEP 1: GET LIGHTCURVES & SPECTRA
    if analysis_type=="get_data": 

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


    ## STEP 2: RE-BIN THE SPECTRA
    elif analysis_type=="group_spectra":
        # Group the spectra -- required before running spectral fits
        group_spectra(min_E_keV , min_counts_chi)
        
  
    ## STEP 3: PLOT THE LIGHT CURVE AND HARDNESS RATIO
    elif analysis_type=="plot_lc_and_hr": 
        plot_lightcurve_and_hr() 


    ## STEP 4: FIT THE SPECTRA WITH NO FIXED PARAMETERS, AND GET INITIAL SPECTRAL FIT RESULTS
    elif analysis_type=="unconstrained_fit": 
        run_spectral_fit(models_unconstrained, uplims_IDs, min_E_keV, fix={}) # fit the spectra
        get_results(models_unconstrained, fixing=False) # output the fit results
        plot_all_spectral_fit_results(models_unconstrained, fixing=False) # plot all the fit results


    ## STEP 5: GET THE PARAMETER VALUES FROM THE SPECTRAL RESULTS, BASED ON CHOSEN MODELS
    elif analysis_type=="get_param_averages": 
        get_results(models_unconstrained, models_indexes_av, fixing=False)


    ## STEP 6: FIT THE SPECTRA AND GET SPECTRAL FIT RESULTS -- WITH PARAMETERS TO FIX
    elif analysis_type=="constrained_fit": 
        run_spectral_fit(models_constrained, uplims_IDs, min_E_keV , fix=fix) # fit the spectra
        get_results(models_constrained, models_indexes, fixing=True) # output the fit results
        ##TODO: Add plotting of all results here too, for comparing them?

  
    ## STEP 7: GET THE UPPER LIMIT FLUX VALUES
    elif analysis_type=="get_uplims": 
        uplims_runner(uplims_IDs, uplims_count_rate_ar, target_coords, nH, gamma)

    ## STEP 8: GET THE FINAL RESULTS, BASED ON CHOSEN MODELS
    elif analysis_type=="get_final_results": 
        plot_spectral_results(models_constrained, models_indexes, uplims_IDs, uplims_MJDs, uplims_MJDs_er, uplims_fluxes, fixing=True)
  

    ## HELPER FUNCTIONS:

    ## F-TEST
    elif analysis_type=="f_test":
        f_test(simple_model, complex_model, fixing=True)


    ## GET THE INDEX RANGE
    # Find the index range in dates_MJD for a particular MJD_range, where MJD_range = [start_MJD, end_MJD)
    # The output index range is inclusive, i.e. [start_index, end_index]
    elif analysis_type=="get_index_range":
        get_index_range(MJD_range) 


    else:
        print("Analysis_type not valid.")
