# Before running, activate the virtual environment:
# source /mnt/users/crookmansourj/swiftenv/bin/activate
# This environment has swiftools 3.0.22

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
(2) plot_lc. This will output the upper limit count rate and MJDs (if any) that will need to be dealt with separately.
(3) fit_spec. If you have an nH value, you can set this to be fixed.
(4) Look at the spectral_fit_results and spectral_fit_residuals folders.
(5) If the fit is failing in some epochs because of too low counts, you may need to determine parameters to fix.
(6) If nH was not specified at the start, run get_nh. This will calculate an average nH value from the fits. 
(7) If some of the fit parameters (e.g. nH) need to be changed, re-run fit_spec with the value of nH. 
(8) If upper limits were identified from the plot_lc routine, run get_uplim_fluxes (with the output nH from the previous step, if needs be).
(9) If needs be, run get_index_range, to get the range of indexes corresponding to a particular MJD range.
(10) Run get_spec_results (with the output of get_uplim_fluxes).
"""

# TODO: 
# scrape of observation IDs from the internet
# scrape target coordinates off the internet
if __name__ in "__main__":
    
    analysis_type = "get_data"

    ## get_data parameters
    # In our case, the spectrum is failing for observation 00016584014 (when I don't include doNotFit=True in the request)
    # Note that Andrew initially had observation 00016584003, but this does not exist
    # Also, for one of the observations, there are two data sets (pc and wt)
    target_coords = [261.930583, -16.205322]
    target_names = ['SwiftJ1727d8m1613', 'Swift J1727.8-1613']
    target_ids = ['00089766', '00016584']
    segments_lc   = [[2, 3, 4, 5, 6, 7, 12], [2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 17, 18, 19]] # observation number
    segments_spec = [[2, 3, 4, 5, 6, 7, 12], [2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 17, 18, 19]] # obs 0001658414 ?
   
    ## fit_spec parameters
    # If we want to fix parameters: fix = {"nh": {"indices": "all", "value": 0.1989}, "gamma": {"indices": [-4, -3, -2, -1], "value": 1.7}, "Tin": {"indices": [-4, -3, -2, -1], "value": 0.5} }
    # If we don't want to fix any parameters: fix={}
    fix = {"nh": {"indices": "all", "value": 0.243}, "gamma": {"indices": [-5, -4, -3, -2, -1], "value": 1.7}, "Tin": {"indices": [-5, -4, -3, -2, -1], "value": 0.5} }

    ## Model choice parameters -- for get_nh and get_spec_results
    # The values are (inclusive) start and end indexes
    models_indexes={'powerlaw': [[6, 17], [19, 21]], 'diskbb': [[0, 5]], 'both':[[]]}

    ## get_uplim_fluxes parameters
    ## This method returns uplims_fluxes
    uplims_obs_id_ar = ["00016584014", "00016584016"]
    uplims_count_rate_ar = [0.01155 , 0.046665] # sometimes [0.011345 0.0475  ]
    nH=0.243
    get_data=False # whether the raw data for the uplims needs to be fetched

    ## get_index_range
    MJD_range = []

    ## get_spec_results parameters
    uplims_MJDs=[60444.44583857,60448.53271186] # sometimes [60444.4458385841 60448.5327118861]
    uplims_fluxes = [4.5305222777569824e-13 , 1.9000751269281060e-12]


    ## GET LIGHTCURVES & SPECTRA
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
    
        
        # Group the spectra -- required before running spectral fits
        group_spectra()
    

    ## PLOT THE LIGHT CURVE AND HARDNESS RATIO
    elif analysis_type=="plot_lc": 
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
        uplims_runner(uplims_obs_id_ar, uplims_count_rate_ar, target_coords, nH, get_data)

    ## GET THE INDEX RANGE
    # Find the index range in dates_MJD for a particular MJD_range, where MJD_range = [start_MJD, end_MJD)
    # The output index range is inclusive, i.e. [start_index, end_index]
    elif analysis_type=="get_index_range":
        get_index_range(MJD_range) 
   
    ## GET THE FINAL RESULTS, BASED ON CHOSEN MODELS
    elif analysis_type=="get_spec_results": 
        plot_spectral_results(uplims_MJDs, uplims_fluxes, models_indexes)
        tabulate_final_results(uplims_MJDs, uplims_fluxes, models_indexes)

    else:
        print("analysis_type not valid.")
