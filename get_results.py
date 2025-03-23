import subprocess
import glob
import json
import os
import numpy as np
np.set_printoptions(threshold=np.inf, precision=16, linewidth=100)
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
from astropy.io import fits
from scipy import stats
#from PyDynamic import interp1d_unc
#rom scipy.interpolate import interp1d
import pandas as pd
pd.set_option('display.width', 1000)  
pd.set_option('display.max_columns', None)  
import warnings
mpl.pyplot.close()
import shutil


####################################################################################################################
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


transitions1 = [58800,58950,59301,59485]
transitions2 = [58800, 58920, 58958, 59295, 59310, 59480, 59492]


####################################################################################################################
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


####################################################################################################################
## LIGHTCURVE PLOTTING FUNCTIONS


def get_lightcurve_data(verbose=True):

    # Get observation IDs
    obs_ids = glob.glob('./lightcurves/*') # overall obsID (not the segments)


    ###########################
    ## Get count rate data

    ##TODO: Clean up numpy array vs list

    # Initialise arrays
    lc_mjd      = np.array([]) # date
    lc_mjd_exp  = np.array([]) # approximate exposure time
    lc_cps      = np.array([]) # count rate
    lc_cps_nerr = np.array([]) # count rate lower error
    lc_cps_perr = np.array([]) # count rate upper error
    obs_type = []

    # Iterate through the files
    for obs_id in obs_ids:

        lc_name = f'{obs_id}/curve.qdp'
        if os.path.exists(f'{obs_id}/curve_incbad.qdp'): # incbad includes data where centroiding could not be done
            lc_name = f'{obs_id}/curve_incbad.qdp'

        headers = []
        with open(lc_name, 'r') as f:
            for line in f: 
                if line.startswith('! '): headers.append(line[1:].strip()) 

        # Data is stored in different table_id values in the .qdp table, so iterate through this
        i = 0
        while True:
            try: 
                data = Table.read(lc_name, format='ascii.qdp', table_id=i)
                lc_mjd = np.append(lc_mjd, data['col1'].data) # creates a copy of the array, so we need to re-assign it to the original variable
                lc_mjd_exp = np.append(lc_mjd_exp, abs(data['col1_perr'].data) + abs(data['col1_nerr'].data))
                lc_cps = np.append(lc_cps, data['col2'].data)
                lc_cps_nerr = np.append(lc_cps_nerr, abs(data['col2_nerr'].data))
                lc_cps_perr = np.append(lc_cps_perr, abs(data['col2_perr'].data))
                obs_type.extend([headers[i]] * len(data))
                i+=1
            except IndexError:
                break
        
    #if verbose: 
    #    print("MJDs: ", mjd2utc(lc_mjd))

    # Order light curve arrays in time
    index = np.argsort(lc_mjd)
    lc_mjd = lc_mjd[index]
    lc_cps = lc_cps[index]
    lc_cps_nerr = lc_cps_nerr[index]
    lc_cps_perr = lc_cps_perr[index]
    lc_mjd_exp = lc_mjd_exp[index]
    obs_type = np.array(obs_type)[index]

    ## Find upper limits
    # These are located at the indexes where the count error is zero
    index_uplims = np.where(lc_cps_nerr == 0.0)
    time_uplims = lc_mjd[index_uplims]
    time_uplims_utc = mjd2utc(time_uplims)
    exp_uplims = lc_mjd_exp[index_uplims]
    cps_uplims = lc_cps[index_uplims]
    obs_type_uplims = obs_type[index_uplims]
    
    ##TODO: Nicer formatting
    if len(time_uplims)!=0 and verbose: 
        print()
        print("Only upper limits were obtained for some observations:")
        print("MJDs for upper limits: \n" +  np.array2string(time_uplims, separator=","))
        print("UTC for upper limits: \n" + np.array2string(time_uplims_utc, separator=","))
        print("Observation duration [MJD] for upper limits \n" +  np.array2string(exp_uplims, separator=","))
        print("Observation type for the upper limits: \n" + np.array2string(obs_type_uplims, separator=","))
        print("Count rates [counts/s] for the upper limits: \n" + np.array2string(cps_uplims, separator=","))
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
                i+=3 # every third table
            except IndexError:
                break

    # Order HR arrays in time
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
    fig, ax = plt.subplots(2, figsize=(20,8),  sharex='col', gridspec_kw={'hspace': 0.1}) 
    fig.set_facecolor('white')
    ax = np.atleast_1d(ax)

    # Light curve
    # Data points:
    ax[0].errorbar(Time(lc_mjd, format='mjd').datetime, lc_cps, yerr = [lc_cps_nerr, lc_cps_perr], fmt='o--', color='k', ms=3) 
    # Upper limits:
    ax[0].errorbar(Time(lc_mjd[index_uplims], format='mjd').datetime, lc_cps[index_uplims ], yerr = 0, fmt='v', color='k', zorder=1000, mec='k', ms = 3, mew=2, mfc='white') 
    ax[0].set_ylabel('$Swift$-XRT Count Rate\n(count s$^{-1}$)')
    
    for t in transitions1: # transition points
        ax[0].axvline(Time(t, format='mjd').datetime, color='yellow', linestyle='--', linewidth=1.5, label="Transitions1")
    for t in transitions2: # transition points
        ax[0].axvline(Time(t, format='mjd').datetime, color='purple', linestyle='--', linewidth=1.5, label="Transitions2")
    ax[0].set_yscale('log')

    # Hardness ratio
    # Filter the array to only include values where the error is low
    mask= hr_err <= 0.05
    ax[1].errorbar(Time(hr_mjd[mask], format='mjd').datetime, hr[mask], yerr =hr_err[mask], fmt='o--', color='k', ms=3)
    ax[1].set_ylabel('$Swift$-XRT HR \n[2-10 keV]/[0.5-2 keV]')
    for t in transitions1: # transition points
        ax[1].axvline(Time(t, format='mjd').datetime, color='yellow', linestyle='--', linewidth=1.5)
    for t in transitions2: # transition points
        ax[1].axvline(Time(t, format='mjd').datetime, color='purple', linestyle='--', linewidth=1.5)
    #ax[1].set_ylim(-0.2, 0.7)


    plt.legend()

    T = np.max(lc_mjd) - np.min(lc_mjd)
    dt = int(T/4)

    FormatAxis(ax, lc_mjd, interval = dt) # use helper function defined above

    # Save the figure
    plots_dir = "./final_results/"
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)
        print(f"{plots_dir } has been created.")
    plt.savefig(plots_dir+'lightcurve_and_hr.png', bbox_inches="tight")


####################################################################################################################
## SPECTRAL RESULTS FUNCTIONS


# Statistical test to see whether extra complexity is justified
# Note that if I am using the 2-step cflux approach and fixing the parameters in the second run, this gives the incorrect #dof
def f_test(simple_model, complex_model, fixing):

     # Load JSON data from file
    filename = './spectral_fit_results'
    if fixing: filename+='_fixing'
    filename+='/'
    json_path = filename+"xrt_spectral_dict.json"
    with open(json_path, 'r') as file:
        data = json.load(file)

    # Output text file
    f = open(filename+"f_test_results.txt", "w")

    # Simple model
    df1 = pd.DataFrame(data[simple_model])
    ID = df1["IDs"].to_numpy()
    chi1 = df1["chi2"].to_numpy()
    dof1 = df1["dof"].to_numpy()

    # Complex model
    # NOTE: with the 2-step model though, this ins't quite right
    df2 = pd.DataFrame(data[complex_model])
    chi2 = df2["chi2"].to_numpy()
    dof2 = df2["dof"].to_numpy()


    for i, (ID, chi1, dof1, chi2, dof2) in enumerate(zip(ID, chi1, dof1, chi2, dof2 )):
        
        f.write(f"{i} ID = {ID}: ")

        chidiff = chi1 - chi2
        dofdiff = dof1 - dof2

        if chi1==-1 and chi2==-1:
            f.write("Neither fit works.")
        elif chi1==-1:
            f.write(f"{complex_model} is preferred; chi^2/dof = {chi2/dof2:.5f} \n")
        elif chi2==-1:
            f.write(f"{simple_model} is preferred; chi^2/dof = {chi1/dof1:.5f}\n")
        
        elif dofdiff==0:
            f.write(f"same #dof ; {simple_model} chi^2/dof = {chi1/dof1:.5f} ; {complex_model} chi^2/dof = {chi2/dof2:.5f} \n")


        else: # both models succeeded in fitting, and have different numbers of dof
        # Test to see whether the added complexity of complex_model is statistically justified
    
            F = (chidiff ) / chi2 * (dof2 / dofdiff)
            p_value = 1 - stats.f.cdf(F, dofdiff, dof2)

            if p_value < 0.001: # the added complexity is statistically justified
                f.write(f"{complex_model} is preferred; chi^2/dof = {chi2/dof2:.5f}\n")
            else:
                f.write(f"{simple_model} is preferred; chi^2/dof = {chi1/dof1:.5f}\n")

    f.close()


def plot_parameter(t, param_values, er_neg, er_pos, param_name, avg, weighted_avg, mod_name, fixing):

    filename = './spectral_fit_results'
    if fixing: filename+='_fixing'
    filename+='/'

    # parameter as a function of time
    plt.figure(figsize=(20, 10))
    #plt.scatter(t, param_values, s = 3, color="blue")
    plt.errorbar(t, param_values, yerr = [er_neg, er_pos], fmt='.', color="blue", markersize=3)
    plt.ylabel(param_name)
    plt.xlabel("MJD")
    plt.axhline(avg, color='red', linestyle='--', linewidth=1.5, label=f'Mean = {avg:.4f}')
    for t in transitions1: # transition points
        plt.axvline(t, color='yellow', linestyle='--', linewidth=1.5)
    for t in transitions2: # transition points
        plt.axvline(Time(t, format='mjd').datetime, color='purple', linestyle='--', linewidth=1.5)
    if weighted_avg!=None: plt.axhline(weighted_avg, color='green', linestyle='--', linewidth=1.5, label=f'Weighted mean = {weighted_avg:.4f}')
    plt.legend()
    plt.savefig(f"{filename}{param_name}_time_series_{mod_name}.png")

    # parameter distribution
    plt.figure(figsize=(12, 10))
    plt.hist(param_values, bins=30, color='blue', alpha=0.7, edgecolor='black')
    plt.axvline(avg, color='red', linestyle='dotted', linewidth=2, label=f'Mean = {avg:.4f}')
    if weighted_avg!=None: plt.axvline(weighted_avg, color='green', linestyle='--', linewidth=1.5, label=f'Weighted mean = {weighted_avg:.4f}')
    plt.xlabel(param_name)
    plt.ylabel('Frequency')
    plt.legend()
    plt.savefig(f"{filename}{param_name}_distribution_{mod_name}.png")


## Print the spectral fit results (json file) to a txt file in table format.
# Once we have an idea which model we would like to use for each observation, we can run this with the models specified
# The model_indexes are an inclusive range and correspond the time-ordered spectra that are fit
def get_results(models, models_indexes=[], count_threshold = 50, fixing=False):

    if any(model not in ['powerlaw', 'pegged_powerlaw', 'diskbb', 'powerlaw+diskbb', 'pegged_powerlaw+diskbb'] for model in models):
        raise ValueError(f"Error: Invalid model(s) found in array. Allowed values are: {['powerlaw', 'pegged_powerlaw', 'diskbb', 'powerlaw+diskbb', 'pegged_powerlaw+diskbb']}")

    # Load JSON data from file
    filename = './spectral_fit_results'
    if fixing: filename+='_fixing'
    filename+='/'
    json_path = filename+"xrt_spectral_dict.json"
    with open(json_path, 'r') as file:
        data = json.load(file)

    # Output text file
    f = open(filename+"fit_outputs.txt", "w")

    # Initialise an empty list to store DataFrames
    dataframes = []

    # Loop through each entry in the dictionary (i.e. for each model)
    for i, model_name in enumerate(models): 

        df = pd.DataFrame(data[model_name])
        df = df[['isot_i'] + ['mjd_i']+ ['mjd_f'] + [col for col in df.columns if col != 'isot_i' and col!='mjd_i' and col!='mjd_f']]
        dataframes.append(df)

        if i==0: # i.e. first model, print the low-count spectra
            counts = df['counts'].to_numpy()
            mask = counts < count_threshold
            obs = df['IDs'].to_numpy()[mask]
            indexes = np.where(mask)[0] 
            f.write("For the following observations, the fitting parameters should be fixed, due to low number of counts: \n")
            f.write("IDs: " + np.array2string(obs, separator=",") +"\n")
            f.write("Indexes: " + np.array2string(indexes, separator=",") + "\n")
            f.write("\n\n\n")

        f.write("All results: \n")
        f.write(model_name +"\n") # model name
        f.write(df.to_string(index=True))
        f.write("\n")

        # Mask when index ranges are specified for this model
        length = len(df['isot_i'])
        try: 
            model_range = models_indexes[i]
            mask_filter = np.full(length, False)
            for period in model_range:
                if period==[]: continue
                start_index, end_index = period[0], period[1]   
                end_index = min(end_index, length - 1)  # Ensure end_index doesn't exceed the bounds of the array 
                mask_filter[start_index: end_index+1] = True # Replace the False with True in those positions in mask
        except: # there are no model_indexes specified
            mask_filter = np.full(length, True)
        
        ## Calculate the average fit parameters -- nH / gamma / Tin
        # To do this, filter the rows where nH is not -1 (i.e. fit was successful), cstat is False (i.e. >=300 counts), and redchi2 is between 0.8 and 1.2
        # In other words, we only use high-quality fits when determining the average of these parameters
        mask_valid = (df['nH'] != -1) & (df['cstat?'] == False) & (df['redchi2'] >= 0.8) & (df['redchi2'] <= 1.2) 
        
        ## Filter the dataframe
        mask = mask_filter & mask_valid
        filtered_df = df[mask]
        t = filtered_df['mjd_i'].to_numpy()
        parameters = ['nH'] + [col for col in ['PhoIndex', 'Tin'] if col in filtered_df.columns] # list of relevant parameter names

        true_indexes = np.where(mask)[0]
        f.write("Indexes used for calculating parameters: \n" + np.array2string(true_indexes, separator=",") + "\n")

        
        ##TODO: Only average when the parameter was not fixed?
        for parameter in parameters:
            
            # Calculate the average of the column for the filtered rows
            values = filtered_df[parameter].to_numpy()
            avg = values.mean()
            f.write(f"Average {parameter}: {avg}\n") 

            # Do weighted average
            weighted_avg = None
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try: # will only work if some of the parameter values haven't been fixed
                    # Compute the weighted average of nh, using the inverse square of the maximum uncertainty values as the weights
                    weighted_avg = np.average(filtered_df[parameter], weights = np.amax((filtered_df[f'{parameter}_neg'], filtered_df[f'{parameter}_pos']), axis=0).astype(float) **(-2)) 
                    neg_er = np.sum(filtered_df[f'{parameter}_neg'].astype(float) ** (-2)) ** (-0.5)
                    pos_er = np.sum(filtered_df[f'{parameter}_pos'].astype(float) ** (-2)) ** (-0.5)
                    # (Note: I could also propagate into this the uncertainty due to the variance of the results)
                    f.write(f"Weighted average {parameter}: {weighted_avg} + {pos_er} - {neg_er}\n") 
                except:
                    f.write(f"Weighted average of {parameter} not calculated.")

            # Plotting results
            # Use all values for range specified for this model, except those for fits that failed
            mask_valid = (df['nH'] != -1) # model failed
            mask = mask_filter & mask_valid
            filtered_df = df[mask]
            t = filtered_df['mjd_i'].to_numpy()
            values = filtered_df[parameter].to_numpy()
            er_neg = filtered_df[parameter+"_neg"].to_numpy()
            er_pos = filtered_df[parameter+"_pos"].to_numpy()
            plot_parameter(t, values, er_neg, er_pos, parameter, avg, weighted_avg, model_name, fixing)
        
        f.write("\n\n")  

    f.close()
    # return tuple(dataframes)  



## TODO: Add different symbol for cstat points?
def plot_spectral_results(models, models_indexes=[], uplims_IDs=[], uplims_MJDs=[], uplims_MJDs_er=[], uplims_fluxes=[], fixing=False):

    mpl.rcParams['xtick.labelbottom'] = False
    colours = ['blue', 'red', 'green', 'purple']

    if not (len(uplims_IDs) == len(uplims_MJDs) == len(uplims_MJDs_er) == len(uplims_fluxes)):
        raise ValueError("The lengths of uplims_IDs, uplims_MJDs, uplims_MJDs_er, and uplims_fluxes must be the same.")

    # Check that none of the defined state ranges overlap
    ranges_list = [np.array(value) for value_list in models_indexes for value in value_list]
    ranges_list = [x for x in ranges_list if x.size > 0] # Remove empty arrays
    if not ranges_okay(ranges_list):
        print(f"The ranges have overlaps.")
        return

    ## Load in the results file
    filename = './spectral_fit_results'
    if fixing: filename+='_fixing'
    filename+='/'
    with open(filename+'xrt_spectral_dict.json', 'r') as j:
        xrt_fit_dict = json.load(j)
    
    # Convert each entry in the nester dictionary to a numpy array
    for key in xrt_fit_dict.keys():
        for keyi in xrt_fit_dict[key].keys():
            xrt_fit_dict[key][keyi] = np.array(xrt_fit_dict[key][keyi])

    all_dates_MJD_start= xrt_fit_dict[models[0]]['mjd_i'] # all the models have the same MJDs; this is the start MJD
    all_dates_MJD_end = xrt_fit_dict[models[0]]['mjd_f'] # all the models have the same MJDs; this is the end MJD
    
    all_dt_MJD = all_dates_MJD_end - all_dates_MJD_start
    all_dates_MJD = 0.5*(all_dates_MJD_start + all_dates_MJD_end) # middle MJD

    all_dt_MJD_alt = xrt_fit_dict[models[0]]['exp [s]']/(24*60*60)
    all_dates_MJD_alt = all_dates_MJD_start + 0.5*all_dt_MJD_alt # middle MJD

    if models_indexes!=[]: # Initialise dataframe to store final results
        df = pd.DataFrame(columns=["obs_id", "middle_mjds", "mjd_range", "flux", "flux_er_neg", "flux_er_pos", "uplims", "model", "cstat_bool", "redchi2"])

    # Make the plot
    fig, ax = plt.subplots(6, figsize=(30,18), sharex='col', gridspec_kw={'hspace': 0.05})#, 'height_ratios': [1.0, 0.4, 1.0, 1.0, 0.4,]})
    fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

    # HR
    # Filter the array to only include values where the error is low
    lc_mjd, lc_cps, lc_cps_nerr, lc_cps_perr, index_uplims, hr_mjd, hr, hr_err = get_lightcurve_data(verbose=False)
    mask= hr_err <= 0.05
    ax[1].errorbar(Time(hr_mjd[mask], format='mjd').datetime, hr[mask], [hr_err[mask], hr_err[mask]], fmt='o', color='k',mfc='black')


    length = len(all_dates_MJD)    
    for i, model in enumerate(models): # make a mask for all the models, just for easier tracking of indexes

        print("Model: ", model)

        mask_valid = xrt_fit_dict[model]['nH']!=-1 # -1 indicates the fit was unsuccessful
        
        # Mask when index ranges are specified
        try: 
            model_range = models_indexes[i]
            mask_filter = np.full(length, False)
            for period in model_range:
                if period==[]: continue
                start_index, end_index = period[0], period[1]   
                end_index = min(end_index, length - 1)  # Ensure end_index doesn't exceed the bounds of the array 
                mask_filter[start_index: end_index+1] = True # Replace the False with True in those positions in mask
        except: # there are no model_indexes specified
            mask_filter = np.full(length, True)
        
        mask = mask_valid & mask_filter
        print("Mask: ", mask)
        print()

        dates_MJD = all_dates_MJD[mask] # middle MJD
        exposures = all_dt_MJD[mask]
        IDs = data['IDs'][mask]

        data = xrt_fit_dict[model]
        
        # Flux
        # For the log values, we need to propagate uncertainties: if y = log_10(x), then x_unc = x * ln(10) * y_unc = 10**y * ln(10) * y_unc
        if model=='pegged_powerlaw': flux, flux_neg, flux_pos = 1e-12*data['norm'][mask], 1e-12*data['norm_neg'][mask], 1e-12*data['norm_pos'][mask]
        elif model=='powerlaw' or model=='diskbb' or model=='powerlaw+diskbb': flux, flux_neg, flux_pos = 10**data['lg10Flux'][mask], 10**data['lg10Flux'][mask]*np.log(10)*data['lg10Flux_neg'][mask], 10**data['lg10Flux'][mask]*np.log(10)*data['lg10Flux_neg'][mask]
        elif model == "pegged_powerlaw+diskbb": flux, flux_neg, flux_pos = 1e-12*data['norm'][mask] + 10**data['lg10Flux'][mask], np.sqrt ( (1e-12*data['norm_neg'][mask])**2 + (10**data['lg10Flux'][mask]*np.log(10)*data['lg10Flux_neg'][mask])**2 ) , np.sqrt( (1e-12*data['norm_pos'][mask])**2 + (10**data['lg10Flux'][mask]*np.log(10)*data['lg10Flux_neg'][mask])**2 )
        ax[0].errorbar(Time(dates_MJD, format='mjd').datetime, flux, [flux_neg, flux_pos], fmt='o',color='k', mfc=colours[i])

        #print(len(data['IDs'][mask]))
        #print()


        if models_indexes!=[]:
            df = pd.concat([df, pd.DataFrame({
            "obs_id": IDs,
            "middle_mjds": dates_MJD,
            "mjd_range": exposures,
            "flux": flux,
            "flux_er_neg": flux_neg,
            "flux_er_pos": flux_pos,
            "uplims": np.nan,
            "model": model,
            "cstat_bool": data['cstat?'][mask],
            "redchi2": data['redchi2'][mask]
            })], ignore_index=True)


        # nH
        nH, nH_neg, nH_pos = data["nH"][mask], data["nH_neg"][mask], data["nH_pos"][mask]
        ax[2].errorbar(Time(dates_MJD, format='mjd').datetime, nH, [nH_neg, nH_pos], fmt='o',color='k', mfc=colours[i])



        # Tin; only for models 'diskbb' and 'pegged_powerlaw+diskbb' and 'powerlaw+diskbb'
        if model=="diskbb" or model=="pegged_powerlaw+diskbb" or model=="powerlaw+diskbb": 
            Tin, Tin_neg, Tin_pos = data['Tin'][mask], data['Tin_neg'][mask], data['Tin_pos'][mask]
            fixed_mask = (Tin_neg == 0) & (Tin_pos == 0)  # plot square if parameter was fixed during fitting
            ax[3].errorbar(Time(dates_MJD[fixed_mask], format='mjd').datetime, Tin[fixed_mask], [Tin_neg[fixed_mask], Tin_pos[fixed_mask]], fmt='s', color='k', mfc=colours[i])
            ax[3].errorbar(Time(dates_MJD[~fixed_mask], format='mjd').datetime, Tin[~fixed_mask], [Tin_neg[~fixed_mask], Tin_pos[~fixed_mask]], fmt='o', color='k', mfc=colours[i])

        # Gamma; only for models 'pegged_powerlaw', 'powerlaw', 'powerlaw+diskbb', 'pegged_powerlaw+diskbb'
        if model=="powerlaw" or model=="pegged_powerlaw" or model=="powerlaw+diskbb" or model=="pegged_powerlaw+diskbb": 
            gamma, gamma_neg, gamma_pos = data['PhoIndex'][mask], data['PhoIndex_neg'][mask], data['PhoIndex_pos'][mask]
            fixed_mask = (gamma_neg == 0) & (gamma_pos == 0) # plot square if parameter was fixed during fitting
            ax[4].errorbar(Time(dates_MJD[fixed_mask], format='mjd').datetime, gamma[fixed_mask], [gamma_neg[fixed_mask], gamma_pos[fixed_mask]], fmt='s', color='k', mfc=colours[i])
            ax[4].errorbar(Time(dates_MJD[~fixed_mask], format='mjd').datetime, gamma[~fixed_mask], [gamma_neg[~fixed_mask], gamma_pos[~fixed_mask]], fmt='s', color='k', mfc=colours[i])

        # Chi^2
        if models_indexes==[]: mask = np.full(length, True) # show all chi^2 results, even when it is above the threshold for the fit and error calculation to be considered successful
        mask_no_cstat = xrt_fit_dict[model]['cstat?'] == False # exclude cstat points
        mask = mask & mask_no_cstat
        chi, dates_MJD = xrt_fit_dict[model]['redchi2'][mask], all_dates_MJD[mask]
        ax[5].errorbar(Time(dates_MJD, format='mjd').datetime, chi, 0.0, fmt='o', color='k', mfc=colours[i], label=model)


        if models_indexes!=[]: # i.e. getting final results
            
            # Store final fits in a separate folder, for evaluation
            final_resid_dir = "./final_residuals/"

            if not os.path.exists(final_resid_dir):
                os.makedirs(final_resid_dir)
                print(f"{final_resid_dir} has been created.")

            plots_dir = "./spectral_fit_residuals"
            if fixing: plots_dir+="_fixing"
            plots_dir+="/"

            for i, ID in enumerate(IDs):

                search_pattern = os.path.join(plots_dir, f"*{ID}_{model}_log.png")
                if i==0: print(search_pattern )

                # Find matching files
                matching_files = glob.glob(search_pattern)

                if matching_files:
                    shutil.copy(matching_files[0], final_resid_dir)
                    print(f"Copied {matching_files[0]} to {final_resid_dir}")
                else:
                    print("No matching file found.")


    if uplims_MJDs: # Add uplims
        ax[0].errorbar(Time(uplims_MJDs, format='mjd').datetime, uplims_fluxes, yerr = 0, fmt='v', color='k',  mec='k', ms = 6, mew=2, mfc='white') # zorder=1000,
        

        if models_indexes!=[]:
            df = pd.concat([df, pd.DataFrame({
            "obs_id": uplims_IDs,
            "middle_mjds": uplims_MJDs,
            "mjd_range": np.array(uplims_MJDs_er) * 2,
            "flux": np.nan,
            "flux_er_neg": np.nan,
            "flux_er_pos": np.nan,
            "uplims": uplims_fluxes,
            "model": None,
            "cstat_bool": None,
            "redchi2": None
            })], ignore_index=True)


    # Set plot constraints
    ax[0].set_yscale('log')
    ax[0].set_ylabel('Flux [1$-$10 keV]\n(erg s$^{-1}$ cm$^{-1}$)')
    ax[1].set_ylabel(r'HR $\left(\frac{[0.5-2\,\text{keV}]}{[2-10\,\text{keV}]}\right)$')
    ax[2].set_ylabel(r'$n_\text{H}$ ($\times 10^{22}$)')
    ax[3].set_ylabel(r'$k_B T_\text{in}$ (keV)')
    ax[4].set_ylabel('$\Gamma$')
    ax[5].set_ylabel(r'$\chi^2_\text{red}$')
    ax[5].set_yscale('log')
    ax[5].legend(fontsize=11)

    for i in range(6):
        for t in transitions1: # transition points
            ax[i].axvline(Time(t, format='mjd').datetime, color='yellow', linestyle='--', linewidth=1.5)
        for t in transitions2: # transition points
            ax[i].axvline(Time(t, format='mjd').datetime, color='purple', linestyle='--', linewidth=1.5)

    T = all_dates_MJD[-1] - all_dates_MJD[0]
    dt = int(T/4)
    FormatAxis(ax, all_dates_MJD, interval = dt) # use helper function defined above

    
    if models_indexes==[]: # plots just for testing
        name = "_".join(models)
        plt.savefig(filename + 'all_fits_'+name+'.png')
    

    else: # final plots, once model index selection has been made
        filename = "./final_results/"

        if not os.path.exists(filename):
            os.makedirs(filename)
            print(f"{filename} has been created.")

        plt.savefig(filename+'final_fit_selection.png')

        # Tabulate the results
        df = df.sort_values(by="middle_mjds").reset_index(drop=True)
        with open(filename + "final_fit_selection.txt", "w") as f:
            f.write(df.to_string(index=False))  
            f.write("\n")
        # Save as a CSV file
        df.to_csv(filename + "final_fit_selection.csv", index=False)
    
    #plt.show()


def plot_all_spectral_fit_results(models, fixing):

    if any(model not in ['powerlaw', 'pegged_powerlaw', 'diskbb', 'powerlaw+diskbb', 'pegged_powerlaw+diskbb']  for model in models):
        raise ValueError(f"Error: Invalid model(s) found in array. Allowed values are: {['powerlaw', 'pegged_powerlaw', 'diskbb', 'powerlaw+diskbb', 'pegged_powerlaw+diskbb'] }")

    for model in models: plot_spectral_results([model], fixing=fixing) # for each model individually
    plot_spectral_results(models, fixing=fixing) # for all the models together





# Find the index range in dates_MJD for a particular MJD_range, where MJD_range = [start_MJD, end_MJD)
# The output index range is inclusive, i.e. [start_index, end_index]
def get_index_range(MJD_range):

    ## Load in the files 

    try: # this will only work if the fit has already been run
        # Note that the data here has been time-ordered
        with open('./spectral_fit_results/xrt_spectral_dict.json', 'r') as j:
            xrt_fit_dict = json.load(j)
        model = list(xrt_fit_dict.keys())[0]
        dates_MJD = xrt_fit_dict[model]['obs_mjds'] # all the models have the same MJDs

    except: 
        raise RuntimeError("Run spectral fitting first.")

    
    start_val, end_val = MJD_range
    length= len(dates_MJD)

    if start_val>end_val:
        print("Start val must be less than end val")
        return
    
    # Find the start index
    if start_val == 0: start_index = 0
    elif start_val == np.inf: start_index = length
    else: start_index = np.searchsorted(dates_MJD, start_val, side='left')
    
    # Find the end index
    if end_val == np.inf: end_index = length
    else: end_index = np.searchsorted(dates_MJD, end_val, side='right')

    end_index-=1 # to be an inclusive range

    print("For the MJD range [{start_val}, {end_val}), use the following indexes (starting from 0) for the xrt_spectral_dict: start index (inclusive)= {start_index} and end index (inclusive) = {end_index}.")
    
    



####################################################################################################################

