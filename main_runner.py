# Before running, activate the virtual environment:
# source /mnt/users/crookmansourj/swiftenv/bin/activate

from get_xrt_from_pipeline import get_xrt_prods, group_spectra
from fit_xrt_spectra import run_spectral_fit
from get_results import plot_lightcurve_and_hr, json_to_dataframes, plot_all_spectral_fit_results, get_nh, plot_spectral_results, tabulate_final_results, get_index_range
import subprocess
import glob
import os
import sys
sys.path.append(os.path.abspath("./uplims_analysis"))
from calculate_uplims import uplims_runner
import numpy as np


"""
Order in which to run:
(1) get_data
(2) plot_lr. Note the upper limit count rate and MJDs (if any) that will need to be dealt with separately.
(3) fit_spec
(4) Look at the spectral_fit_results and spectral_fit_residuals folders.
(5) If the fit is failing in some epochs because of too low counts, I may need to determine parameters to fix.
(6) If nH was not specified at the start, run get_nh. 
(7) If some of the fit parameters (e.g. nH) need to be changed, re-run fit_spec with the value of nH. 
(8) If upper limits were identified from the plot_lr routine, run get_uplim_fluxes (with the output nH from the previous step, if needs be).
(9) If needs be, run get_index_range, to get the range of indexes corresponding to a particular MJD range.
(10) Run get_spec_results (with the output of get_uplim_fluxes).
"""
# TODO: 
# scraping of observation IDs from the internet
# sort out issue with centroiding with some of the spectra
# scrape target coordinates off the internet
if __name__ in "__main__":
    
    analysis_type = "fit_spec"

    ## get_data parameters
    # In our case, centroiding is failing for the spectra for observation 00016584014
    target_coords = [261.930583, -16.205322]
    target_names = ['SwiftJ1727d8m1613', 'Swift J1727.8-1613']
    target_ids = ['00089766', '00016584']
    segments = [[2, 3, 4, 5, 6, 7, 12], [2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 17, 18, 19]] # observation number
    centroiding = [np.full(7, True, dtype=bool), np.concatenate((np.full(10, True, dtype=bool), [False], np.full(4, True, dtype=bool)))]

    ## Determine the ID numbers
    #sorted_segments = [sorted(segment) for segment in segments] # Sort each segment in ascending order
    #IDs = []
    #for target_id, segment in zip(target_ids, sorted_segments):
    #    IDs.extend([f"{target_id}{seg:03}" for seg in segment])
    #IDs=np.array(IDs)
    
    ## fit_spec parameters
    # If we want to fix parameters: fix = {"nh": {"indices": "all", "value": 0.1989}, "gamma": {"indices": [-4, -3, -2, -1], "value": 1.7}, "Tin": {"indices": [-4, -3, -2, -1], "value": 0.5} }
    # If we don't want to fix any parameters: fix={}
    fix = {"nh": {"indices": "all", "value": 0.24}, "gamma": {"indices": [-4, -3, -2, -1], "value": 1.7}, "Tin": {"indices": [-4, -3, -2, -1], "value": 0.5} }

    ## Model choice parameters -- for get_nh and get_spec_results
    # The values are (inclusive) start and end indexes
    models_indexes={'powerlaw': [[7, 20]], 'diskbb': [[0, 5]], 'both':[[6, 6]]}

    ## get_uplims parameters
    ## returns uplims_fluxes
    uplims_obs_id_ar = ["00016584014", "00016584016"]
    uplims_count_rate_ar = [0.01155 , 0.046665]
    nH=0.2

    ## get_index_range
    MJD_range = []

    ## get_spec_results parameters
    uplims_MJDs=[60444.44583857,60448.53271186]
    #uplims_fluxes=[2.63826081e-12, 2.37226154e-12]
    uplims_fluxes = [4.43677120e-13 , 1.86076696e-12]


    ## GET LIGHTCURVES & SPECTRA
    if analysis_type=="get_data": 

        # Iterate through Target name(s), id(s), and segments and get xrt_prods
        for k, target_name in enumerate(target_names):
            print("TARGET: ", target_name)
            print("Getting light curve:")
            get_xrt_prods(target_ids[k], target_name, target_coords, segments[k], centroid=centroiding, prod_type = 'lightcurve')
            print()
            print("Getting spectrum:")
            get_xrt_prods(target_ids[k], target_name, target_coords, segments[k], centroid=centroiding, prod_type = 'spectrum')
            print("-------------------------------------------------------------------------------------------------------")
            pass
        group_spectra()
    
    ## PLOT THE LIGHT CURVE AND HARDNESS RATIO
    elif analysis_type=="plot_lr": 
        plot_lightcurve_and_hr()

    ## FIT THE SPECTRA
    elif analysis_type=="fit_spec": 
        run_spectral_fit(fix) # fit the spectra
        json_to_dataframes() # output the fit results
        plot_all_spectral_fit_results() # plot all the fit results

    ## GET NH VALUE, BASED ON CHOSEN MODELS
    elif analysis_type=="get_nh": 
        get_nh(models_indexes)

    ## GET THE UPPER LIMIT FLUX VALUES
    elif analysis_type=="get_uplim_fluxes": 
        uplims_runner(uplims_obs_id_ar, uplims_count_rate_ar, target_coords, nH)

    elif analysis_type=="get_index_range":
        get_index_range(MJD_range) 
   
    ## GET THE FINAL RESULTS, BASED ON CHOSEN MODELS
    elif analysis_type=="get_spec_results": 
        plot_spectral_results(uplims_MJDs, uplims_fluxes, models_indexes)
        tabulate_final_results(uplims_MJDs, uplims_fluxes, models_indexes)

    else:
        print("analysis_type not valid.")
