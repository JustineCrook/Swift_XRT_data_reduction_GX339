# Before running, activate the virtual environment:
# source /mnt/users/crookmansourj/swiftenv/bin/activate


from get_xrt_from_pipeline import get_xrt_prods, group_spectra
from fit_xrt_spectra import run_spectral_fit
from get_results import plot_lightcurve_and_hr, json_to_dataframes, plot_all_spectral_fit_results, plot_spectral_results, tabulate_final_results
import subprocess
import glob
import os
import sys


# TODO: 
# scraping of observation IDs from the internet
# sort out issue with centroiding with some of the spectra
# srape target coordinates off the internet
def get_data():

    # Define the Swift XRT parameters for Swift J1727.8-2609
    target_coords = [261.930583, -16.205322]

    ## LIGHTCURVES
    target_names = ['SwiftJ1727d8m1613', 'Swift J1727.8-1613']
    target_ids = ['00089766', '00016584']
    segments = [[2, 3, 4, 5, 6, 7, 12], [2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 17, 18, 19]] # observation number

    # Iterate through Target name(s), id(s), and segments and get xrt_prods
    for k, target_name in enumerate(target_names):
        print("TARGET: ", target_name)
        get_xrt_prods(target_ids[k], target_name, target_coords, segments[k], centroid=True, prod_type = 'lightcurve')
        print("-------------------------------------------------------------------------------------------------------")
        pass
    
    ## SPECTRA
    # I split this in two groups for Swift J1727.8-1613, because the centroiding does not work for the spectrum algorithm for observations 14, 16, 18 (too faint)
    #target_names = ['SwiftJ1727d8m1613', 'Swift J1727.8-1613', 'Swift J1727.8-1613']
    target_names = ['Swift J1727.8-1613']
    #target_ids = ['00089766', '00016584','00016584']
    target_ids = ['00016584']
    #segments = [[2, 3, 4, 5, 6, 7, 12], [2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 17, 19], [18]] # observation number
    segments = [[14]] # 14,16
    #centroiding = [True, True, False]
    centroiding = [False]

    # Iterate through Target name(s), id(s), and segments and get xrt_prods
    for k, target_name in enumerate(target_names):
        print("TARGET: ", target_name)
        get_xrt_prods(target_ids[k], target_name, target_coords, segments[k], centroid=centroiding[k], prod_type = 'spectrum')
        print("-------------------------------------------------------------------------------------------------------")
        pass

    group_spectra()


def plot_lr():
    plot_lightcurve_and_hr()



# TODO:
# Add option to specify parameters to fix
def fit_spec():

    run_spectral_fit()
    json_to_dataframes()
    plot_all_spectral_fit_results()


## TODO
# Returns uplims
def get_uplim_fluxes():
    print("hello")



# TODO:
# Add option to specify date/index ranges -- otherwise the default calculated ones are used
# Add option to add uplims
def plot_spec():
    plot_spectral_results()
    tabulate_final_results()



       
if __name__ in "__main__":
    
    analysis_type = input("Enter type of analysis -- get_data, plot_lr, fit_spec, get_uplim_fluxes, plot_spec: ")
    
    ## TODO: options to add other inputs
    if analysis_type=="fit_spec": fit_spec()
    if analysis_type=="plot_lr": plot_lr() 
    else: print("Error")
