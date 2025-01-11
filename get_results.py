import subprocess
import glob
import json
import os
import numpy as np
np.set_printoptions(threshold=np.inf, precision=16)
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pickle as p
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition,mark_inset)
from astropy.time import Time
import scipy.stats 
import scipy.special 
#from scipy.special import iv
from astropy.constants import c
from astropy.table import Table
#from PyDynamic import interp1d_unc
#rom scipy.interpolate import interp1d
import pandas as pd
pd.set_option('display.width', 1000)  
pd.set_option('display.max_columns', None)  
import warnings
mpl.pyplot.close()



## TODO: Clean up plotting functions, and make helper functions -- e.g. maybe with a getdata
##TODO: Clean up json_to_dataframes and get_nh



##########################################################
## Matplotlib formatting 

mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['font.size'] = 14
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['ytick.minor.visible'] = True
mpl.rcParams['ytick.labelsize'] = 'medium'
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['ytick.minor.width'] = 1.5
mpl.rcParams['ytick.major.size'] = 6.0
mpl.rcParams['ytick.minor.size'] = 3.0
mpl.rcParams['xtick.top'] = True
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['xtick.labelsize'] = 'medium'
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['xtick.minor.width'] = 1.5
mpl.rcParams['xtick.major.size'] = 6.0
mpl.rcParams['xtick.minor.size'] = 3.0
mpl.rcParams['axes.linewidth'] = 1.5



##########################################################
## Helper functions

def FormatAxis(ax, mjd, dt = 10, interval = 60):
    '''Function for putting UTC on top of axis'''

    ax[0].set_xlabel('Observing Date (UTC)', fontfamily='serif')    
    ax[0].xaxis.set_major_locator(mdates.DayLocator(interval=interval))
    ax[0].set_xlim(Time(mjd[0] - dt, format='mjd').datetime, Time(mjd[-1] + dt, format='mjd').datetime)
    ax[0].xaxis.set_label_position('top') 
    xformatter = mdates.DateFormatter('%Y-%m-%d')
    plt.gcf().axes[0].xaxis.set_major_formatter(xformatter) #ax[0].xaxis.set_major_formatter(xformatter)
    ax[0].tick_params(axis='x', which='major',rotation=10,labeltop=True, labelbottom=False)

    # Format secondary x-axis
    mjd_ax = ax[-1].secondary_xaxis('bottom', functions=(plot2mjd, mjd2plot))
    mjd_ax.set_xlabel('Observing Date (MJD)', fontfamily='serif')  
    mjd_ax.tick_params(which='major', direction='in', length = 0.0, width = 0.0)
    plt.draw()

    # Extract the labels
    mjd_ticks = []
    labels = ax[0].get_xticklabels(which='major')
    for lab in labels:
        mjd_ticks.append(lab.get_text() + 'T00:00:00')

    # Line up MJD and Datetime labels 
    mjd_ticks = (Time(mjd_ticks, format='isot').mjd).astype(int)
    mjd_ax.set_xticks(mjd_ticks, labels = mjd_ticks)
    

def plot2mjd(t):
    '''Convert from matplotlib plot date to mjd'''
    return Time(t, format="plot_date", scale='utc').mjd


def mjd2plot(mjd):
    '''Convert from mjd to matplotlib plot'''
    return Time(mjd, format="mjd", scale='utc').plot_date

def mjd2utc(mjd):
    # Convert MJD to UTC using astropy
    t = Time(mjd, format='mjd', scale='utc')
    return t.iso  # Returns in ISO format (YYYY-MM-DD HH:MM:SS.sss)

def iso2mjd(iso_dates):
    # Convert ISO dates to MJD using astropy Time
    times = Time(iso_dates, format='isot', scale='utc')
    return times.mjd



# Input: 2D list -- a list of ranges in the form [start_index, end_index], where the indexes are inclusive
def ranges_okay(ranges):

    # Check validity of each range
    for range_pair in ranges:
        if len(range_pair) != 2:
            return False  # Each range must have exactly 2 values
        if range_pair[0] > range_pair[1]:
            return False  # The first value must be less or equal to than the second

    # Sort ranges based on the start value
    sorted_ranges = sorted(ranges, key=lambda x: x[0])
    #print(sorted_ranges)
    
    for i in range(len(sorted_ranges) - 1):
        # Check if the current range overlaps with the next range
        if sorted_ranges[i][1] > sorted_ranges[i + 1][0]:
            return False  # Overlap found

    return True  # No overlaps and all ranges are valid




##########################################################
## LIGHTCURVE PLOTTING FUNCTIONS


def get_lightcurve_data(verbose=True):

    # Get observation IDs
    obs_ids = glob.glob('./lightcurves/*')

    ###########################
    ## Get count data

    # Initialise arrays
    lc_mjd      = np.array([]) # date
    lc_cps      = np.array([]) # counts
    lc_cps_nerr = np.array([]) # counts lower error
    lc_cps_perr = np.array([]) # counts upper error

    # Iterate through the files
    for obs_id in obs_ids:

        lc_name = f'{obs_id}/curve.qdp'

        if os.path.exists(f'{obs_id}/curve_incbad.qdp'):
            lc_name = f'{obs_id}/curve_incbad.qdp'

        i = 0
        while True:
            try:
                data = Table.read(lc_name, format='ascii.qdp', table_id=i)
                lc_mjd = np.append(lc_mjd, data['col1'].data) # creates a copy of the array, so we need to re-assign it to the original variable
                lc_cps = np.append(lc_cps, data['col2'].data)
                lc_cps_nerr = np.append(lc_cps_nerr, abs(data['col2_nerr'].data))
                lc_cps_perr = np.append(lc_cps_perr, abs(data['col2_perr'].data))
                i+=1
            except IndexError:
                break
        
    if verbose: 
        print(mjd2utc(lc_mjd))

    # Order light curve arrays
    index = np.argsort(lc_mjd)
    lc_mjd = lc_mjd[index]
    lc_cps = lc_cps[index]
    lc_cps_nerr = lc_cps_nerr[index]
    lc_cps_perr = lc_cps_perr[index]


    ## Find upper limits
    # These are the indexes where the count error is zero.
    index_uplims = np.where(lc_cps_nerr == 0.0)
    time_uplims = lc_mjd[index_uplims]
    cps_uplims = lc_cps[index_uplims]
    time_uplims_utc = mjd2utc(time_uplims)

    if len(time_uplims)!=0 and verbose: 
        print()
        print("Only upper limits were obtained for some observations:")
        print("MJDs for upper limits: ", time_uplims)
        print("UTC for upper limits: ", time_uplims_utc)
        print("Count rates [counts/s] for the upper limits: ", cps_uplims)
        print("Get the obsIDs for the upper limits from the following website: https://www.swift.psu.edu/operations/obsSchedule.php")


    ###########################
    ## Get hardness ratios (HR)

    # Initialise arrays
    hr_mjd  = np.array([])
    hr      = np.array([])
    hr_err = np.array([])

    # Iterate through the files
    for obs_id in obs_ids:
        hr_name = f'{obs_id}/hardrat.qdp'
        
        # The hardness ratio is only in every third table -- the other two are the hard-band data and soft-band data 
        i = 2
        while True:
            try:
                data = Table.read(hr_name, format='ascii.qdp', table_id=i)
                hr_mjd = np.append(hr_mjd, data['col1'].data)
                hr = np.append(hr, data['col2'].data)
                hr_err= np.append(hr_err, abs(data['col2_err'].data))
                i+=3
            except IndexError:
                break

    # Order HR arrays
    index = np.argsort(hr_mjd)
    hr_mjd = hr_mjd[index]
    hr = hr[index]
    hr_err = hr_err[index]

    return lc_mjd, lc_cps, lc_cps_nerr, lc_cps_perr, index_uplims, hr_mjd, hr, hr_err



def plot_lightcurve_and_hr():
    ## Plot the light curves and hardness ratios

    mpl.rcParams['xtick.labelbottom'] = False

    # Get the data
    lc_mjd, lc_cps, lc_cps_nerr, lc_cps_perr, index_uplims, hr_mjd, hr, hr_err = get_lightcurve_data()

    # Set up figure
    fig, ax = plt.subplots(2, figsize=(8,8),  sharex='col', gridspec_kw={'hspace': 0.1}) 
    fig.set_facecolor('white')
    ax = np.atleast_1d(ax)

    # Light curve
    ax[0].errorbar(Time(lc_mjd, format='mjd').datetime, lc_cps, yerr = [lc_cps_nerr, lc_cps_perr], fmt='o--', color='k') # data points
    ax[0].errorbar(Time(lc_mjd[index_uplims], format='mjd').datetime, lc_cps[index_uplims ], yerr = 0, fmt='v', color='k', zorder=1000, mec='k', ms = 10, mew=2, mfc='white') # upper limits
    ax[0].set_ylabel('$Swift$-XRT Count Rate\n(count s$^{-1}$)')
    ax[0].set_yscale('log')

    # Hardness ratio
    mask= hr_err <= 0.05
    #hr_err_filtered = np.where(filter_mask, hr_err, np.nan)
    #hr_filtered = np.where(filter_mask, hr, np.nan)
    ax[1].errorbar(Time(hr_mjd[mask], format='mjd').datetime, hr[mask], yerr =hr_err[mask], fmt='o--', color='k')
    ax[1].set_ylabel('$Swift$-XRT HR \n[2-10 keV]/[0.5-2 keV]')
    ax[1].set_ylim(-0.2, 0.7)

    FormatAxis(ax, lc_mjd, interval = 30) # use helper function defined above

    plots_dir = "./plots/"
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)
        print(f"{plots_dir } has been created.")
    plt.savefig('./plots/lightcurve_and_hr.png', bbox_inches="tight")


##########################################################
## INITIAL SPECTRAL RESULTS FUNCTIONS

## Print the spectral fit results (json file) to a txt file in table format.
def json_to_dataframes():

    # Load JSON data from file
    json_path = "./spectral_fit_results/xrt_spectral_dict.json"
    with open(json_path, 'r') as file:
        data = json.load(file)
    
    # Initialise an empty list to store DataFrames
    dataframes = []

    f = open("./spectral_fit_results/fit_outputs.txt", "w")

    # Loop through each entry in the dictionary (i.e. for each model)
    models = ["powerlaw", "diskbb", "both"]
    for model_name in models: 
        df = pd.DataFrame(data[model_name])
        df = df[['obs_isot'] + ['obs_mjds']+ [col for col in df.columns if col != 'obs_isot' and col!='obs_mjds']]
        dataframes.append(df)
        f.write(model_name +"\n") # model name
        f.write(df.to_string(index=True))

        # Filter the rows where nH is not -1, cstat is False, and redchi2 is between 0.8 and 1.2
        mask = (df['nH'] != -1) & (df['cstat?'] == False) & (df['redchi2/redcstat'] >= 0.8) & (df['redchi2/redcstat'] <= 1.2)
        filtered_df = df[mask]
        # Calculate the average of the nH column for the filtered rows
        avg_nh = filtered_df['nH'].mean()
        # Write the average to the file
        f.write(f"\nAverage nH: {avg_nh}\n") 

        # Do weighted average
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try: # will only work if some of the values haven't been fixed
                avg_nh = np.average(filtered_df['nH'], weights = np.amax((filtered_df['nH_neg'], filtered_df['nH_pos']), axis=0).astype(float) **(-2)) # Computes the weighted average of nh, using the inverse square of the maximum uncretainty values as the weights.
                nh_neg_er = np.sum(filtered_df['nH_neg'].astype(float) ** (-2)) ** (-0.5)
                nh_pos_er = np.sum(filtered_df['nH_pos'].astype(float) ** (-2)) ** (-0.5)
                # (Note: I could also propagate into this the uncertainty due to the variance of the results)
                f.write(f"Weighted average nH: {avg_nh} + {nh_pos_er} - {nh_neg_er}\n") 
            except:
                print("Weighted average not calculated.")

        f.write("\n\n")  

    f.close()
    # return tuple(dataframes)  


## TODO: Possibly use a mask (same as above) for the HR
## TODO: Add different symbol for cstat points.
def plot_all_spectral_fit_results_helper(model_indexes=[0,1,2]):

    mpl.rcParams['xtick.labelbottom'] = False
    
    colours = {0:'blue', 1:'red', 2:'green'}
    all_models = {0: "powerlaw", 1: "diskbb", 2: "both"}
    model_names = [all_models[index] for index in model_indexes] # model_names = np.array(list(all_models.values()))

    ## Load in the files 
    with open('./spectral_fit_results/xrt_spectral_dict.json', 'r') as j:
        xrt_fit_dict = json.load(j)
    
    ## Convert each entry in the nester dictionary to a numpy array
    for key in xrt_fit_dict.keys():
        for keyi in xrt_fit_dict[key].keys():
            xrt_fit_dict[key][keyi] = np.array(xrt_fit_dict[key][keyi])

    dates_MJD = xrt_fit_dict['powerlaw']['obs_mjds'] # all the models have the same MJDs
    length = len(dates_MJD)

    # Make the plot
    fig, ax = plt.subplots(5, figsize=(10,11), sharex='col', gridspec_kw={'hspace': 0.1})#, 'height_ratios': [1.0, 0.4, 1.0, 1.0, 0.4,]})

    # Set plot constraints
    ax[0].set_yscale('log')
    ax[0].set_ylabel('Flux [1$-$10 keV]\n(erg s$^{-1}$ cm$^{-1}$)')
    ax[1].set_ylabel(r'HR $\left(\frac{[0.5-2\,\text{keV}]}{[2-10\,\text{keV}]}\right)$')
    ax[2].set_ylabel(r'$k_B T_\text{in}$ (keV)')
    ax[3].set_ylabel('$\Gamma$')
    ax[4].set_ylabel(r'$\chi^2_\text{red}$')
    ax[4].set_yscale('log')

    # Note that -1 indicates the fit was unsuccessful.
    masks =[]
    for model_name in all_models.values(): # make a mask for all the models, just for easier tracking of indexes
        mask = xrt_fit_dict[model_name]['nH']!=-1
        masks.append(mask)

    all_dates_MJD = xrt_fit_dict['powerlaw']['obs_mjds']
    
    # Flux
    # For the log values, we need to propagate uncertainties: if y = log_10(x), then x_unc = x * ln(10) * y_unc = 10**y * ln(10) * y_unc
    for i, model_name in zip(model_indexes, model_names):
        mask = masks[i]
        data = xrt_fit_dict[model_name]
        dates_MJD = all_dates_MJD[mask]
        if model_name=='powerlaw': flux, flux_neg, flux_pos = 1e-12*data['norm'][mask], 1e-12*data['norm_neg'][mask], 1e-12*data['norm_pos'][mask]
        elif model_name=='diskbb': flux, flux_neg, flux_pos = 10**data['lg10Flux'][mask], 10**data['lg10Flux'][mask]*np.log(10)*data['lg10Flux_neg'][mask], 10**data['lg10Flux'][mask]*np.log(10)*data['lg10Flux_neg'][mask]
        elif model_name == "both": flux, flux_neg, flux_pos = 1e-12*data['norm'][mask] + 10**data['lg10Flux'][mask], np.sqrt ( (1e-12*data['norm_neg'][mask])**2 + (10**data['lg10Flux'][mask]*np.log(10)*data['lg10Flux_neg'][mask])**2 ) , np.sqrt( (1e-12*data['norm_pos'][mask])**2 + (10**data['lg10Flux'][mask]*np.log(10)*data['lg10Flux_neg'][mask])**2 )
        ax[0].errorbar(Time(dates_MJD, format='mjd').datetime, flux, [flux_neg, flux_pos], fmt='o:',color='k', mfc=colours[i])

    # HR
    lc_mjd, lc_cps, lc_cps_nerr, lc_cps_perr, index_uplims, hr_mjd, hr, hr_err = get_lightcurve_data(verbose=False)
    ax[1].errorbar(Time(hr_mjd, format='mjd').datetime, hr, [hr_err, hr_err], fmt='o:', color='k',mfc='black')

    # Tin
    indexes= [i for i in model_indexes if i in [1,2]] # disbb and both
    for i in indexes: 
        model_name, mask = all_models[i], masks[i]
        data = xrt_fit_dict[model_name]
        dates_MJD, Tin, Tin_neg, Tin_pos = all_dates_MJD[mask], data['Tin'][mask], data['Tin_neg'][mask], data['Tin_pos'][mask]
        ax[2].errorbar(Time(dates_MJD, format='mjd').datetime, Tin, [Tin_neg, Tin_pos], fmt='o:',color='k', mfc=colours[i])

    # Gamma
    indexes= [i for i in model_indexes if i in [0,2]] # powerlaw and both
    for i in indexes: 
        model_name, mask = all_models[i], masks[i]
        data = xrt_fit_dict[model_name]
        dates_MJD, gamma, gamma_neg, gamma_pos = all_dates_MJD[mask], data['PhoIndex'][mask], data['PhoIndex_neg'][mask], data['PhoIndex_pos'][mask]
        ax[3].errorbar(Time(dates_MJD, format='mjd').datetime, gamma, [gamma_neg, gamma_pos], fmt='o:', color='k', mfc=colours[i])

    # Chi^2
    for i, model_name in zip(model_indexes, model_names):
        mask = xrt_fit_dict[model_name]['cstat?'] == False
        chi, dates_MJD = xrt_fit_dict[model_name]['redchi2/redcstat'][mask], all_dates_MJD[mask]
        ax[4].errorbar(Time(dates_MJD, format='mjd').datetime, chi, 0.0, fmt='o:', color='k', mfc=colours[i], label=model_name)

    FormatAxis(ax, all_dates_MJD, interval = 30) # use helper function defined above

    ax[4].legend(fontsize=11)
    
    str = "_".join(model_names)
    plt.savefig('./spectral_fit_results/all_fits_'+str+'.png')
    #plt.show()


def plot_all_spectral_fit_results():
    plot_all_spectral_fit_results_helper(model_indexes=[0,1,2])
    plot_all_spectral_fit_results_helper(model_indexes=[0])
    plot_all_spectral_fit_results_helper(model_indexes=[1])
    plot_all_spectral_fit_results_helper(model_indexes=[2])



# Find the index range in dates_MJD for a particular MJD_range, where MJD_range = [start_MJD, end_MJD)
# The output index range is inclusive, i.e. [start_index, end_index]
def get_index_range(MJD_range):

    ## Load in the files 
    with open('./spectral_fit_results/xrt_spectral_dict.json', 'r') as j:
        xrt_fit_dict = json.load(j)
    dates_MJD = xrt_fit_dict['powerlaw']['obs_mjds'] # all the models have the same MJDs
    
    start_val, end_val = MJD_range
    length= len(dates_MJD)

    if start_val>end_val:
        print("Start val must be less than end val")
        return
    
    # Find the start index
    if start_val == 0:
        start_index = 0
    elif start_val == np.inf:
        start_index = length
    else:
        start_index = np.searchsorted(dates_MJD, start_val, side='left')
    
    # Find the end index
    if end_val == np.inf:
        end_index = length
    else:
        end_index = np.searchsorted(dates_MJD, end_val, side='right')

    end_index-=1 # to be inclusive

    print("For the MJD range [{start_val}, {end_val}), use the following indexes (starting from 0) for the xrt_spectral_dict: start index (inclusive)= {start_index} and end index (inclusive) = {end_index}.")
    
    


# Once we have an idea which model we would like to use for each observation, we can run this with the models specified
def get_nh(models_indexes={'powerlaw': [[7, 20]], 'diskbb': [[0, 5]], 'both':[[6, 6]]}):

    # Load JSON data from file
    json_path = "./spectral_fit_results/xrt_spectral_dict.json"
    with open(json_path, 'r') as file:
        data = json.load(file)

    f = open("./spectral_fit_results/nh.txt", "w")

    models = ["powerlaw", "diskbb", "both"]
    for model_name in models: 
        df = pd.DataFrame(data[model_name])
        length = df.shape[0]
        mask = np.full(length, False)
        model_range = models_indexes[model_name]
        for period in model_range:
            start_index, end_index = period[0], period[1]  
            end_index = min(end_index, length - 1)  # Ensure end_index doesn't exceed the bounds of the array
            mask[start_index: end_index+1] = True # Replace the False with True in those positions in mask
        filtered_df = df[mask]
        f.write(model_name +"\n") 
        f.write(filtered_df.to_string(index=True))
        # Calculate the average of the nH column for the filtered rows
        avg_nh = filtered_df['nH'].mean()
        # Write the average to the file
        f.write(f"\nAverage nH: {avg_nh}\n") 

        # Do weighted average
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try: # will only work if some of the values haven't been fixed
                avg_nh = np.average(filtered_df['nH'], weights = np.amax((filtered_df['nH_neg'], filtered_df['nH_pos']), axis=0).astype(float) **(-2)) # Computes the weighted average of nh, using the inverse square of the maximum uncretainty values as the weights.
                nh_neg_er = np.sum(filtered_df['nH_neg'].astype(float) ** (-2)) ** (-0.5)
                nh_pos_er = np.sum(filtered_df['nH_pos'].astype(float) ** (-2)) ** (-0.5)
                # (Note: I could also propagate into this the uncertainty due to the variance of the results)
                f.write(f"Weighted average nH: {avg_nh} + {nh_pos_er} - {nh_neg_er}\n") 
            except:
                print("Weighted average not calculated.")
        f.write("\n\n")  

    f.close()


##########################################################
## FINAL SPECTRAL RESULTS FUNCTIONS


## TODO: Possibly use a mask (same as above) for the HR
## TODO: Remove the cstat plotting
## TODO: Add different symbol for cstat points.
def plot_spectral_results(uplims_MJDs=[60444.44583857,60448.53271186], uplims_fluxes=[2.63826081e-12, 2.37226154e-12], models_indexes={'powerlaw': [[7, 20]], 'diskbb': [[0, 5]], 'both':[[6, 6]]}): 
    
    mpl.rcParams['xtick.labelbottom'] = False

    colours = ['blue', 'red', 'green']
    model_names = ["powerlaw", "diskbb", "both"]

    if len(uplims_MJDs)!=len(uplims_fluxes):
        print("uplims_MJDs and uplims_fluxes should have the same length.")
        return

    # Check that none of the defined state ranges overlap
    ranges_list = [np.array(value) for value_list in models_indexes.values() for value in value_list]
    # Remove empty arrays
    ranges_list = [x for x in ranges_list if x.size > 0]
    #print(ranges_list)
    if not ranges_okay(ranges_list):
        print(f"The  ranges for splitting the observation have overlaps.")
        return

    ## Load in the files 
    with open('./spectral_fit_results/xrt_spectral_dict.json', 'r') as j:
        xrt_fit_dict = json.load(j)
    
    ## Convert each entry in the nester dictionary to a numpy array
    for key in xrt_fit_dict.keys():
        for keyi in xrt_fit_dict[key].keys():
            xrt_fit_dict[key][keyi] = np.array(xrt_fit_dict[key][keyi])

    dates_MJD = xrt_fit_dict['powerlaw']['obs_mjds']
    length = len(dates_MJD)

    # Make the plot
    fig, ax = plt.subplots(5, figsize=(12,17), sharex='col', gridspec_kw={'hspace': 0.1})#, 'height_ratios': [1.0, 0.4, 1.0, 1.0, 0.4,]})


    masks = []
    for model_name in model_names: 
        mask = np.full(length, False)
        model_range = models_indexes[model_name]
        for period in model_range:
            #print(period)
            if period==[]: continue
            start_index, end_index = period[0], period[1]   
            end_index = min(end_index, length - 1)  # Ensure end_index doesn't exceed the bounds of the array 
            mask[start_index: end_index+1] = True # Replace the False with True in those positions in mask
        masks.append(mask)

    all_dates_MJD = xrt_fit_dict['powerlaw']['obs_mjds']

    # Flux
    for i, model_name in enumerate(['powerlaw', 'diskbb', 'both']):
        mask = masks[i]
        data = xrt_fit_dict[model_name]
        dates_MJD = all_dates_MJD[mask]
        if model_name=='powerlaw': flux, flux_neg, flux_pos = 1e-12*data['norm'][mask], 1e-12*data['norm_neg'][mask], 1e-12*data['norm_pos'][mask]
        elif model_name=='diskbb': flux, flux_neg, flux_pos = 10**data['lg10Flux'][mask], 10**data['lg10Flux'][mask]*np.log(10)*data['lg10Flux_neg'][mask], 10**data['lg10Flux'][mask]*np.log(10)*data['lg10Flux_neg'][mask]
        elif model_name == "both": flux, flux_neg, flux_pos = 1e-12*data['norm'][mask] + 10**data['lg10Flux'][mask], np.sqrt ( (1e-12*data['norm_neg'][mask])**2 + (10**data['lg10Flux'][mask]*np.log(10)*data['lg10Flux_neg'][mask])**2 ) , np.sqrt( (1e-12*data['norm_pos'][mask])**2 + (10**data['lg10Flux'][mask]*np.log(10)*data['lg10Flux_neg'][mask])**2 )
        ax[0].errorbar(Time(dates_MJD, format='mjd').datetime, flux, [flux_neg, flux_pos], fmt='o',color='k', mfc=colours[i])


    # Add uplims
    ax[0].errorbar(Time(uplims_MJDs, format='mjd').datetime, uplims_fluxes, yerr = 0, fmt='v', color='k',  mec='k', ms = 6, mew=2, mfc='white') # zorder=1000,

    # HR
    lc_mjd, lc_cps, lc_cps_nerr, lc_cps_perr, index_uplims, hr_mjd, hr, hr_err = get_lightcurve_data(verbose=False)
    ax[1].errorbar(Time(hr_mjd, format='mjd').datetime, hr, [hr_err, hr_err], fmt='o', color='k',mfc='black')
    

    # Tin
    for i in [1,2]: # disbb and both
        model_name, mask = model_names[i], masks[i]
        data = xrt_fit_dict[model_name]
        dates_MJD, Tin, Tin_neg, Tin_pos = all_dates_MJD[mask], data['Tin'][mask], data['Tin_neg'][mask], data['Tin_pos'][mask]
        ax[2].errorbar(Time(dates_MJD, format='mjd').datetime, Tin, [Tin_neg, Tin_pos], fmt='o',color='k', mfc=colours[i])

    # Gamma
    for i in [0,2]: # powerlaw and both
        model_name, mask = model_names[i], masks[i]
        data = xrt_fit_dict[model_name]
        dates_MJD, gamma, gamma_neg, gamma_pos = all_dates_MJD[mask], data['PhoIndex'][mask], data['PhoIndex_neg'][mask], data['PhoIndex_pos'][mask]
        ax[3].errorbar(Time(dates_MJD, format='mjd').datetime, gamma, [gamma_neg, gamma_pos], fmt='o', color='k', mfc=colours[i])

    # Chi^2
    for i, model_name in enumerate(['powerlaw', 'diskbb', 'both']):
        mask = masks[i]
        mask_no_cstat = xrt_fit_dict[model_name]['cstat?'] == False
        mask = mask & mask_no_cstat
        chi = xrt_fit_dict[model_name]['redchi2/redcstat'][mask]
        dates_MJD = all_dates_MJD[mask]
        ax[4].errorbar(Time(dates_MJD, format='mjd').datetime, chi, 0.0, fmt='o', color='k', mfc=colours[i], label=model_name)
    
    # Set plot constraints
    ax[0].set_ylabel('Flux [1$-$10 keV]\n(erg s$^{-1}$ cm$^{-1}$)')
    ax[0].set_yscale('log')
    ax[1].set_ylabel(r'HR $\left(\frac{[0.5-2\,\text{keV}]}{[2-10\,\text{keV}]}\right)$')
    ax[2].set_ylabel(r'$k_B T_\text{in}$ (keV)')
    ax[3].set_ylabel('$\Gamma$')
    ax[4].set_ylabel(r'$\chi^2_\text{red}$/Cstat$_\text{red}$')
    ax[4].set_yscale('log')
    ax[4].legend(fontsize=11)
    
    FormatAxis(ax, all_dates_MJD, interval = 30) # use helper function defined above

    
    plots_dir = "./plots/"
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)
        print(f"{plots_dir } has been created.")
    plt.savefig('./plots/flux.png')
    #plt.show()



# For the ranges [start, end] is inclusive (i.e. includes the start and end index value)
## TODO: Only show chi^2 for chosen ranges, as above
## TODO: Add HR, as above
## TODO: Remove the cstat plotting
## TODO: Add different symbol for cstat points.
## TODO: I might also want to use different markers/colours for the various models used. In this case, since errrobar only accepts a single symbol, and I want the flux points to be connected, I will have to first plot the line connecting all the points, and then plot the errobars in a loop.
def plot_spectral_results_connected_lines(uplims_MJDs=[60444.44583857,60448.53271186], uplims_fluxes=[2.63826081e-12, 2.37226154e-12], models_indexes={'powerlaw': [[7, 20]], 'diskbb': [[0, 5]], 'both':[[6, 6]]}): 

    mpl.rcParams['xtick.labelbottom'] = False

    colours = ['white', 'black', 'grey']
    models = ["powerlaw", "diskbb", "both"]

    if len(uplims_MJDs)!=len(uplims_fluxes):
        print("uplims_MJDs and uplims_fluxes should have the same length.")
        return

    # Check that none of the defined state ranges overlap
    ranges_list = [np.array(value) for value_list in models_indexes.values() for value in value_list]
    if not ranges_okay(ranges_list):
        print(f"The  ranges for splitting the observation have overlaps.")
        return

    ## Load in the files 
    with open('./spectral_results/xrt_spectral_dict.json', 'r') as j:
        xrt_fit_dict = json.load(j)
    
    ## Convert each entry in the nester dictionary to a numpy array
    for key in xrt_fit_dict.keys():
        for keyi in xrt_fit_dict[key].keys():
            xrt_fit_dict[key][keyi] = np.array(xrt_fit_dict[key][keyi])

    dates_MJD = xrt_fit_dict['powerlaw']['obs_mjds']
    length = len(dates_MJD)

    # Make the plot
    fig, ax = plt.subplots(4, figsize=(10,11), sharex='col', gridspec_kw={'hspace': 0.1})#, 'height_ratios': [1.0, 0.4, 1.0, 1.0, 0.4,]})

    # Set plot constraints
    ax[0].set_yscale('log')
    ax[0].set_ylabel('Flux [1$-$10 keV]\n(erg s$^{-1}$ cm$^{-1}$)')
    ax[1].set_ylabel(r'$k_B T_\text{in}$ (keV)')
    ax[2].set_ylabel('$\Gamma$')
    ax[3].set_ylabel(r'$\chi^2_\text{red}$/Cstat$_\text{red}$')
    ax[3].set_yscale('log')
    
    flux_all, flux_all_pos, flux_all_neg, Tin_all, Tin_all_pos, Tin_all_neg, gamma_all, gamma_all_pos, gamma_all_neg = np.full((9,length), np.nan)

    for i, model_name in enumerate(models):
        print(model_name)
        mask = np.full(length, False)
        
        model_range = models_indexes[model_name]
        for period in model_range:
            start_index, end_index = period[0], period[1]    
            mask[start_index: end_index+1] = True # Replace the False with True in those positions in mask
        print(mask)

        # Flux
        if model_name=="powerlaw": flux, flux_neg, flux_pos = 1e-12 *xrt_fit_dict['powerlaw']['norm'], 1e-12 *xrt_fit_dict['powerlaw']['norm_neg'], 1e-12 *xrt_fit_dict['powerlaw']['norm_pos']
        if model_name=="diskbb": flux, flux_neg, flux_pos = 10 **xrt_fit_dict['diskbb']['lg10Flux'], 10 **xrt_fit_dict['diskbb']['lg10Flux'] * np.log(10) * xrt_fit_dict['diskbb']['lg10Flux_neg'], 10 **xrt_fit_dict['diskbb']['lg10Flux']* np.log(10) * xrt_fit_dict['diskbb']['lg10Flux_neg']
        if model_name=="both": flux, flux_neg, flux_pos = 1e-12 *xrt_fit_dict['both']['norm'] + 10 **xrt_fit_dict['both']['lg10Flux'], 1e-12 *xrt_fit_dict['both']['norm_neg'] + 10 **xrt_fit_dict['both']['lg10Flux'] * np.log(10) * xrt_fit_dict['both']['lg10Flux_neg'], 1e-12 *xrt_fit_dict['both']['norm_pos'] + 10 **xrt_fit_dict['both']['lg10Flux']* np.log(10) * xrt_fit_dict['both']['lg10Flux_neg']
        flux_all[mask] = flux[mask]
        flux_all_pos[mask] = flux_pos[mask]
        flux_all_neg[mask] = flux_neg[mask]
        
        # Tin
        if model_name=="diskbb" or model_name=="both": # diskbb
            Tin, Tin_neg, Tin_pos = xrt_fit_dict[model_name]['Tin'], xrt_fit_dict[model_name]['Tin_neg'], xrt_fit_dict[model_name]['Tin_pos']
        else: 
            Tin, Tin_neg, Tin_pos = np.full(length, np.nan), np.full(length, np.nan), np.full(length, np.nan)
        Tin_all[mask] = Tin[mask]
        Tin_all_neg[mask] = Tin_neg[mask]
        Tin_all_pos[mask] = Tin_pos[mask]

        # Gamma
        if model_name=="powerlaw" or model_name=="both": # powerlaw
            gamma, gamma_neg, gamma_pos = xrt_fit_dict[model_name]['PhoIndex'], xrt_fit_dict[model_name]['PhoIndex_neg'], xrt_fit_dict[model_name]['PhoIndex_pos']
        else: 
            gamma, gamma_neg, gamma_pos = np.full(length, np.nan), np.full(length, np.nan), np.full(length, np.nan)    
        gamma_all[mask] = gamma[mask]
        gamma_all_neg[mask] = gamma_neg[mask]
        gamma_all_pos[mask] = gamma_pos[mask]

        # Chi^2
        chi= xrt_fit_dict[model_name]['redchi2/redcstat']
        ax[3].errorbar(Time(dates_MJD, format='mjd').datetime, chi, 0.0, fmt='o:', color='k', mfc=colours[i], label=model_name)


        # Transition lines
        transitions = np.where(mask[:-1] != mask[1:])[0]
        for i in transitions:
            date1 = dates_MJD[i]
            date2 = dates_MJD[i + 1]
            avg_mjd = (date1 + date2) / 2
            avg_datetime = Time(avg_mjd, format='mjd').datetime
            for ax_i in ax: ax_i.axvline(avg_datetime, ls='--', color='k')
    
    
    def plot(j, y, yer_neg, yer_pos):
        indexes = np.isnan(y) # some of the values have remained NaN
        ax[j].errorbar(Time(dates_MJD[~indexes], format='mjd').datetime, y[~indexes], [yer_neg[~indexes], yer_pos[~indexes]], fmt='o:', color='k')
        indexes_fixed = (yer_neg == 0) & (yer_pos == 0) # parameters that were fixed
        ax[j].scatter(Time(dates_MJD[indexes_fixed], format='mjd').datetime, y[indexes_fixed], marker='s', color='k')
    plot(1, Tin_all, Tin_all_neg, Tin_all_pos)
    plot(2, gamma_all, gamma_all_neg, gamma_all_pos)

    
    ## Add uplims
    if uplims_fluxes!=None:
        flux_all= np.concatenate((flux_all, uplims_fluxes))
        dates_MJD= np.concatenate((dates_MJD, uplims_MJDs))
        flux_all_neg = np.concatenate((flux_all_neg, np.full(len(uplims_fluxes), np.nan)))
        flux_all_pos = np.concatenate((flux_all_pos, np.full(len(uplims_fluxes), np.nan)))
        indexes_ordered = np.argsort(dates_MJD)
        dates_MJD = dates_MJD[indexes_ordered]
        flux_all = flux_all[indexes_ordered]
    
    # Have to do it this way so that the upper limits connect with the rest of the points
    ax[0].errorbar(Time(dates_MJD, format='mjd').datetime, flux_all, [flux_all_neg, flux_all_pos], fmt='o:', color='k')
    # Plot uplims on top:
    ax[0].errorbar(Time(uplims_MJDs, format='mjd').datetime, uplims_fluxes, yerr = 0, fmt='v', color='k', zorder=1000, mec='k', ms = 10, mew=2, mfc='white') 

    FormatAxis(ax, dates_MJD, interval = 30) # use helper function defined above

    ax[3].legend()
    plt.savefig('./plots/flux_connected_lines.png')
    #plt.show()



## TODO
## Take in date ranges for the models and then tabulate final results
def tabulate_final_results(uplims_MJDs, uplims_fluxes, models_indexes):
    print("TODO: tabulate final results")


