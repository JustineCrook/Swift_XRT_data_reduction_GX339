import json
import glob
import os
import sys 
import subprocess
import numpy as np
from astropy.io import fits
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd

from xspec import * 


def iso2mjd(iso_dates):
    # Convert ISO dates to MJD using astropy Time
    times = Time(iso_dates, format='isot', scale='utc')
    return times.mjd


## Initialise different spectral models, based on the input model_index.
# Note: The model-fitting for the more complex compound models is more complicated, and hands-on fitting would be preferred here. 
# Since we are only interested in the flux (not a detailed fit), the simpler models are likely appropriate, at the cost of chi-squared.
# In some cases, emission/absorption lines are visible, and we may also need a more manual fit. 
def initialise_model(model, flux_guess, fix_names, fix_values, parameters=None):
    
    ## Initialise the parameters
    # We can enter ranges for the parameter by passing in a string containing <param>,<delta>,<hard min>,<soft min>,<soft max>,<hard max> where:
    # param: trial param value used initially in fit
    # delta: step size used in the numerical determination of the derivatives during the fitting process. This value may be overriden for all parameters by the xset delta command option, which will apply a proportional rather than a fixed delta.
    # Or, if we specify a value and then -1, it is fixes to this value.

    # Emin and Emax are set to the energy range over which we want the flux to be calculated.
    Emin = 1.0 
    Emax = 10.0

    # Defaults for starting parameters and ranges
    # Note: These bounds are very wide; e.g., Gamma 0 to 4, whereas the BH LMXBs are almost always 1-3. 
    # So in our fitting, if any of the values get pegged at the bounds, it either means; (i) the spectrum is insensitive to the pegged model, and thus it might not be necessary; (ii) we are stuck in a local minimum and may require more hands-on spectral fitting.
    if parameters==None:
        parameters = {
        "nh": '0.2,,0.005,0.005,10,10',  # absorption, range 0.005-10 <<<<<<<<<<<<<<<<<<<<<<<<<
        "gamma": '1.7,,0.0,0.0,4.0,4.0',  # spectral index, range 0-4 
        "Tin": '0.5,,0.05,0.05,4.0,4.0',  # temperature, range 0.05-4
        "norm1": '1.0',
        "norm2": '1.0'
        }
    # Set parameters to the values specified in 'fix'
    for name, value in zip(fix_names, fix_values):
        if name in parameters: parameters[name] = f"{value} -1"

    nh = parameters["nh"]
    gamma = parameters["gamma"]
    Tin = parameters["Tin"]
    norm1 = parameters["norm1"]
    norm2 = parameters["norm2"]


    ## Initialise the models. We make use of:
    # tbabs (see XSPEC manual pp. 348-349): Absorption due to the ISM including molecules and grains. Allows the user just to vary the hydrogen column
    # pgpwrlw (see XSPEC manual pg. 297): Power law with pegged normalisation. A(E) = K E^{-alpha}
    # cflux (see XSPEC manual pp. 359-360): A convolution model to calculate the flux of other model components.
    # diskbb (see XSPEC manual pg. 249): Multiple blackbody accretion disk model.

    ## Absorbed power law with pegged normalisation
    # par1: nH; equivalent hydrogen column, in units of 10^{22} atoms/cm^2 
    # par2: alpha; photon index of pwrlw (dimensionless)
    # par3: lower peg energy range (keV)
    # par4: upper peg energy range (keV)
    # par5: norm; flux (in units of 10^{-12} ergs/cm^2/s over the energy par2-par3). If par3 = par4, it is the flux in micro-Jy at par3.
    if model == "pegged_powerlaw":
        mod = 'tbabs * pegpwrlw' 
        initial_pars = {1:nh, 
                        2: gamma, 3:Emin, 4:Emax, 5:f'{flux_guess}'}

    ## Absorbed power law 
    # par1: nH; equivalent hydrogen column, in units of 10^{22} atoms/cm^2 
    # par2: Emin; Minimum energy over which flux is calculated.
    # par3: Emax; Maximum energy over which flux is calculated.
    # par4: lg10Flux; log (base 10) flux in erg/cm^2/s
    # par5: alpha; photon index of pwrlw (dimensionless)
    # par6: norm; K, photons/keV/cm2/s at 1 keV.
    elif model == "powerlaw":
        pwrlw_flux = np.log10(flux_guess * 1e-12)
        mod = 'tbabs * cflux * powerlaw' 
        initial_pars = {1:nh, 
                        2: Emin, 3: Emax, 4: f'{pwrlw_flux}',
                        5:gamma, 6:'1.0 -1'}

    ## Absorbed disk blackbody
    # par1: nH; equivalent hydrogen column, in units of 10^{22} atoms/cm^2 
    # par2: Emin; Minimum energy over which flux is calculated.
    # par3: Emax; Maximum energy over which flux is calculated.
    # par4: lg10Flux; log (base 10) flux in erg/cm^2/s
    # par5: temperature of the inner disk radius (keV)
    # par6: norm; (Rin/ D10)^2 cos0, where where Rin is an apparent inner disk radius in km, D10 the distance to the source in units of 10 kpc, and the angle of the disk ( = 0 is face-on). On the correction factor between the apparent inner disk radius and the realistic radius, see e.g. Kubota et al. 1998.
    elif model == "diskbb":
        bb_flux = np.log10(flux_guess * 1e-12)
        mod = 'tbabs * cflux *  diskbb'
        initial_pars = {1:nh, 
                        2:Emin, 3:Emax, 4:f'{bb_flux}', 5: Tin, 6:'1.0 -1'}
 
    ## Absorbed powerlaw + disk blackbody
    # par1: nH; equivalent hydrogen column, in units of 10^{22} atoms/cm^2 
    # par2: alpha; photon index of pwrlw (dimensionless)
    # par3: lower peg energy range (keV) for pwrlw
    # par4: upper peg energy range (keV) for pwrlw
    # par5: norm; pwrlw flux (in units of 10^{-12} ergs/cm^2/s over the energy par2-par3). If par3 = par4, it is the flux in micro-Jy at par3.
    # par6: Emin; Minimum energy over which bb flux is calculated.
    # par7: Emax; Maximum energy over which bb flux is calculated.
    # par8: lg10Flux; log (base 10) bb flux in erg/cm^2/s
    # par9: temperature of the inner disk radius (keV)
    # par10: norm for bb; (Rin/ D10)^2 cos0, where where Rin is an apparent inner disk radius in km, D10 the distance to the source in units of 10 kpc, and the angle of the disk ( = 0 is face-on). On the correction factor between the apparent inner disk radius and the realistic radius, see e.g. Kubota et al. 1998.
    elif model == "pegged_powerlaw+diskbb":
        mod = 'tbabs * (pegpwrlw + cflux *  diskbb)'
        bb_flux = np.log10(flux_guess * 0.1 * 1e-12)
        pwrlw_flux = flux_guess*0.9
        initial_pars = {1:nh, 
                        2: gamma, 3:Emin, 4:Emax, 5:f'{pwrlw_flux}',
                        6:Emin, 7:Emax, 8:f'{bb_flux}', 9: Tin, 10:'1.0 -1'}


    ## Absorbed powerlaw + disk blackbody
    # par1: nH; equivalent hydrogen column, in units of 10^{22} atoms/cm^2 
    # par2: alpha; photon index of pwrlw (dimensionless)
    # par3: norm; K, photons/keV/cm2/s at 1 keV.
    # par4: temperature of the inner disk radius (keV)
    # par5: norm for bb; (Rin/ D10)^2 cos0, where where Rin is an apparent inner disk radius in km, D10 the distance to the source in units of 10 kpc, and the angle of the disk ( = 0 is face-on). On the correction factor between the apparent inner disk radius and the realistic radius, see e.g. Kubota et al. 1998.
    elif model=="powerlaw+diskbb":
        mod = 'tbabs*(powerlaw + diskbb)'
        initial_pars = {1:nh, 
                        2:gamma, 3:'1.0',
                        4: Tin, 5:'1.0'}


    ## Absorbed powerlaw + disk blackbody
    # par1: nH; equivalent hydrogen column, in units of 10^{22} atoms/cm^2 
    # par2: Emin; Minimum energy over which flux is calculated.
    # par3: Emax; Maximum energy over which flux is calculated.
    # par4: lg10Flux; log (base 10) flux in erg/cm^2/s
    # par5: alpha; photon index of pwrlw (dimensionless)
    # par6: norm; K, photons/keV/cm2/s at 1 keV.
    # par7: temperature of the inner disk radius (keV)
    # par8: norm for bb; (Rin/ D10)^2 cos0, where where Rin is an apparent inner disk radius in km, D10 the distance to the source in units of 10 kpc, and the angle of the disk ( = 0 is face-on). On the correction factor between the apparent inner disk radius and the realistic radius, see e.g. Kubota et al. 1998.
    elif model=="cflux_(powerlaw+diskbb)":
        mod = 'tbabs*cflux(powerlaw + diskbb)'
        flux = np.log10(flux_guess * 1e-12)
        initial_pars = {1:nh, 
                        2: Emin, 3: Emax, 4: f'{flux}',
                        5:gamma, 6:norm1,
                        7: Tin, 8:norm2}

    else:
        sys.exit('Incorrect index: Please input either 0, 1, or 2')

    print(initial_pars)

    return mod, initial_pars # return the model and parameters


## Helper method to plot the residuals for the fits
def plot_resid(spectrum_name, mod_name, fixing):

    mpl.rcParams['xtick.labelbottom'] = True

    plots_dir = "./spectral_fit_residuals"
    if fixing: plots_dir+="_fixing"
    plots_dir+="/"
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)

    Plot('data resid') # data from most recent Fit that was performed

    # Get coordinates from plot:
    x = Plot.x()
    rates = Plot.y()
    yer = Plot.yErr()
    xer = Plot.xErr()
    folded = Plot.model()
    #resids = np.array(rates) - np.array(folded)
    resids = Plot.y(1,2)

    dataLabels = Plot.labels(1)
    residLabels = Plot.labels(2)

    fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={'height_ratios': [3, 1]})

    ax[0].errorbar(x, rates,xerr=xer, yerr=yer, fmt='ro',ms=1, label="data", elinewidth=0.2)
    ax[0].plot(x, folded, 'b', label="model", linewidth=1, zorder=100)
    ax[0].set_ylabel(dataLabels[1]) #ax[0].set_ylabel(r'counts/cm$^2$/sec/keV')
    ax[0].legend(fontsize=11)

    ax[1].plot(x, resids, 'g', label="residuals (data-model)", linewidth=1)
    ax[1].set_ylabel(residLabels[1])
    ax[1].set_xlabel(residLabels[0]) #ax[1].set_xlabel('Energy [keV]')
    ax[1].legend(fontsize=11)

    ax[0].set_title("Spectrum: " + spectrum_name + " & model: " + mod_name)

    # Save the results
    plt.savefig(plots_dir+ spectrum_name+"_"+mod_name)
    
    # Also save results on log scales
    ax[0].set_yscale("log")
    ax[0].set_xscale("log")
    plt.savefig(plots_dir+ spectrum_name+"_"+mod_name+"_log")

    # Clear the plot
    plt.clf()



# The parameter fix specifies which parameters need to be fixed during the fit and for which spectral files.
# The parameter models specifies which models to fit to the data -- the code fits each of these to all the data.
# Note: We don't fit the entire energy range of data because the response of the instrument is worse at the edges of the band.
##TODO: Checking the fix parameters
##TODO: Fill all with -1 by default (of correct length), to clean up the code
def run_spectral_fit(models= ['pegged_powerlaw'], uplim_files=[], low_energy_fit_bound_keV = 0.5, count_threshold = 50, fix={}, add_systematic = False):

    if any(model not in ['powerlaw', 'pegged_powerlaw', 'diskbb', 'powerlaw+diskbb', 'pegged_powerlaw+diskbb'] for model in models):
        raise ValueError(f"Error: Invalid model(s) found in array. Allowed values are: {['powerlaw', 'pegged_powerlaw', 'diskbb', 'powerlaw+diskbb', 'pegged_powerlaw+diskbb']}")
    
    if fix =={}: fixing = False
    else: fixing = True

    # Folder for results
    folder = './spectral_fit_results'
    if fixing: folder+='_fixing'
    folder+='/'
    if not os.path.exists(folder):
        os.makedirs(folder)

    # Initialise XSPEC
    # Xset: storage class for XSPEC settings
    Xset.abund = "wilm"    # set the abundance table used in the plasma emission and photoelectric absorption models; 'wilm' is an XSPEC built-in abundance table
    Xset.xsect = "vern"    # set the photoelectric absorption cross-sections
    Xset.chatter = 0 # set the console chatter level... change this to 10 for more chatter
    Xset.logChatter = 20
    Xset.openLog("xspec.log") 
    #Xset.parallel.error = 1 # previously 10, i.e. use up to 10 parallel processes during the Fit.error() command runs... however, this caused issues with the parameter saving being done before the values were updated during the error calculation

    # Fit: manager class for setting properties and running a fit
    Fit.query = "yes"      # when nIterations is reached, continue the fit without stopping to query
    Fit.nIterations = 100  # fit iterations for each attempt
    Fit.statTest = "chi" # default

    # Plot: manager class for performing XSPEC plots
    Plot.xAxis = "KeV"     # x-axis for plotting set to energy instead of channel
    Plot.device = '/null'
    Plot.yLog = True

    # Load in the spectral files
    spectral_files = glob.glob('./spectra/*final.pi') # the spectra after the grouping script has been run
    n_spectra = len(spectral_files) 

    
    # Initialise and fill variables
    obs_isots  = [] # to hold the observation dates
    end_obs_isots = []
    bin_counts = [] # to hold number of counts in grouped bin -- in order to decide whether chi^2 or C-stats should be used
    tot_counts = []
    exp = []
    for i, spectrum_file in enumerate(spectral_files):
        header = fits.getheader(spectrum_file, ext = 1)
        obs_isots.append(header['DATE-OBS']) # start date
        end_obs_isots.append(header['DATE-END']) # start date
        bin_counts.append(header['COUNTGRP']) # min number of counts in the bins
        tot_counts.append(header['COUNTS']) #... should be same as s.rate[0] * s.exposure
        exp.append(header['EXPOSURE']) # exposure length... same as s.exposure
    # Convert dates to MJDs
    obs_mjds = iso2mjd(np.array(obs_isots)) 
    end_obs_mjds = iso2mjd(np.array(end_obs_isots)) 


    # Order lists by time
    sort_index = np.argsort(obs_mjds)
    spectral_files = np.array(spectral_files)[sort_index]
    obs_isots = np.array(obs_isots)[sort_index]
    obs_mjds = obs_mjds[sort_index]
    end_obs_isots = np.array(end_obs_isots)[sort_index]
    end_obs_mjds = end_obs_mjds[sort_index]
    bin_counts = np.array(bin_counts)[sort_index]
    tot_counts = np.array(tot_counts)[sort_index]
    exp = np.array(exp)[sort_index]
    
    # Also get other parameters of interest
    IDs = np.array([sf[14:-8] for sf in spectral_files])
    mode = np.array([file[-10:-8] for file in spectral_files]) # PC or WT
    cstat = np.array([True if tot_count==1 else False for tot_count in tot_counts])

    # Filtering -- ignore spectra with very low counts
    # We do not want to fit the observations identified as uplims from the lightcurve generation (as we deal with these separately)
    # We also decide to not fit WT observations with fewer than 15 counts (as these are probably spurious)
    # Also do not fit PC observations with fewer than 3 counts
    # Lastly, we do not fit observations where
    # The 'good' boolean is specified above
    mask_detections = np.array([False if ID in uplim_files else True for ID in IDs])
    mask_valid = ((mode == 'pc') & (tot_counts >= 3)) | ((mode == 'wt') & (tot_counts >= 15))
    mask = mask_detections & mask_valid 

    # Print out obs_isot, obs_mjds, IDs, counts, exp [s] for the ignored spectra
    xrt_dict_no_fit = {"obs_isot": obs_isots[~mask],"obs_mjds": obs_mjds[~mask],"IDs": IDs[~mask],"counts": tot_counts[~mask],"exp": exp[~mask]}
    df = pd.DataFrame(xrt_dict_no_fit)
    df.to_string(folder+"not_fit.txt", index=False, justify="left")

    # Get the parameters for the observations that we will fit
    spectral_files = spectral_files[mask]
    obs_isots = obs_isots[mask]
    obs_mjds = obs_mjds[mask]
    end_obs_isots = end_obs_isots[mask]
    end_obs_mjds = end_obs_mjds[mask]
    bin_counts = bin_counts[mask]
    tot_counts = tot_counts[mask]
    exp = exp[mask]
    IDs = IDs[mask]
    mode = mode[mask]
    cstat = cstat[mask]


    # Initialise dictionary to hold results
    xrt_dict = {model: {"IDs": IDs.tolist(), "isot_i": obs_isots.tolist(), "mjd_i": obs_mjds.tolist(), "mjd_f": end_obs_mjds.tolist(), "exp [s]": exp.tolist(), "counts": tot_counts.tolist(), "cstat?":cstat.tolist(), "chi2": [], "dof": [], "redchi2": []} for model in models}

    ## Iterate through each spectrum (i.e. each observation), fit it, get the fit parameters, and then save these.
    # AllData: container for all loaded data sets (objects of class Spectrum).
    # AllModels: container for all Model objects.
    # XSPEC models are defined by creating a Model object.
    # Note: getattr(object, 'x') is completely equivalent to object.x, but used in the case we don't know exactly what 'x' is.
    for k, spectrum_file in enumerate(spectral_files):

        print('\n\n', spectrum_file)

        # Process the 'fix' input dictionary
        fix_names = [] # names of parameters to fix for this spectrum
        fix_values = [] # values to fix these parameters to for this spectrum
        for param, details in fix.items():
            indices = details["indices"]
            value = details["value"]
            # Only include this parameter if the index corresponds to the current spectrum
            if indices == "all" or k in indices or k in np.array(indices)+n_spectra: # the last check is in case some of the values were negative (i.e. counted from the end of the array)
                fix_names.append(param)
                fix_values.append(value)

        
        # Load the spectrum into the AllData object, removing any previously loaded data sets... We also run AllData.clear() before, just to be sure.
        # If the rmf and arf files are located in the same folder (which they are), these will automatically be loaded.
        AllData.clear()
        AllData(spectrum_file)
        # AllData.show() # check current state of the AllData container

        # Define range of energy to use -- i.e. ignore certain channels. This is done for all the loaded spectra.
        # The floating point values are assumed to be in keV, as this is what was set for the Plot.xAxis value above.
        AllData.ignore('*:10.0-**')
        AllData.ignore(f'*:**-{low_energy_fit_bound_keV}')
        AllData.notice(f'*:{low_energy_fit_bound_keV}-10.0') # so that edge points are used
        AllData.ignore('bad') # ignore data that is bad -- using the quality column

        # If too much data is ignored, the spectrum won't fit 
        # Also print out a warning 
        if (len(AllData(1).ignored) == len(AllData(1).noticed)): print("ERROR: too much bad data")

        # Choose whether to use chi-squared or cash statistics, based on the number of bins in the data
        # statMethod: type of statistic to minimise 
        # Chi-squared stats assumes that all the spectral channels are Gaussian distributed and that the estimate for the variance is uncorrelated with the observed counts.
        # If the data are Poisson, then the C-statistic is better.
        if bin_counts[k] == 1:
            Fit.statMethod = 'cstat' # Poisson data
            print("Using C-statistics.")
        else:
            Fit.statMethod = 'chi' # Gaussian data
            print("Using chi squared statistics.")

        # Next, we need an initial guess for the flux (flux_guess), to be used in the spectral fitting routine. 
        # We will get an estimate of the flux using a simple powerlaw fit.
        # Define the spectral model -- absorbed power law
        # 1: nH (10^{22} atoms/cm^2); 2: photon index; 3: lower peg energy range (keV); 4: upper peg energy range (keV); 5: norm i.e. flux (10^{-12} ergs/cm^2/s)
        AllModels.clear()
        mod1 = Model('tbabs * pegpwrlw')
        mod1.setPars({1:'0.2 -1', 2:'2.0 -1', 3:1.0, 4:10.0, 5:'1000.0'})
        Fit.delta = 1e-2 # controls the convergence criterion for the fitting algorithm
        try: 
            Fit.perform()
            flux_guess = mod1.pegpwrlw.norm.values[0] # flux (in units of 10^{-12} ergs/cm^2/s)
        except: flux_guess = 1 # faint data... i.e. assume flux = 1 * 10^{-12} ergs/cm^2/s


        ## Iterate through the model versions, fitting each type, and appending the results to an output dictionary
        mod_2 = False # whether the second model succeeded in the 2-step model approach
        for mod_name in models:

            print('\n', mod_name)

            AllModels.clear() # AllModels represents the collection of all currently defined models in the XSPEC session

            # Initialise the model and systematics 
            mod, initial_pars = initialise_model(mod_name, flux_guess, fix_names, fix_values)
            mod_obj = Model(mod)
            mod_obj.setPars(initial_pars) # set the parameters to the initial parameters specified in the initialisation
            if add_systematic: AllModels.systematic = 0.03 # add systematic error
            AllModels.show() # for checking

            # Do not fit if tot_counts is less than count_threshold plus the fit parameters haven't been fixed
            fit = True
            if tot_counts[k] < count_threshold:
                if "nh" not in fix_names:
                    fit = False
                elif mod_name in {'powerlaw', 'pegged_powerlaw', 'powerlaw+diskbb', 'pegged_powerlaw+diskbb'} and "gamma" not in fix_names:
                    fit = False
                elif mod_name == 'diskbb' and "Tin" not in fix_names:
                    fit = False

            if not fit:
                print("Fitting will not be performed because #counts is below the threshold. Fix the fit parameters for this observation.")

            if fit: # Perform fitting
            
                try:
                    # Adjust the fit precision, and perform the fit of the spectrum. 
                    # The first pass allows the fit to reach a close approximation of the best-fit parameters quickly. 
                    # Each subsequent pass with a smaller delta refines the fit further.
                    for delt in [1e-2, 1e-3, 1e-4, 1e-5]:
                        Fit.delta = delt
                        Fit.perform()

                    # 2-step powerlaw+diskbb
                    # Note that there is no 'editmod' in pyxspec
                    # I am fixing nH, gamma, Tin for the second run because I have noticed that when the powerlaw component is so much less dominant than the diskbb one, the cflux convolution does not retain similar values for the parameters.
                    # I checked that fixing them doesn't make much difference to the flux uncertainty
                    ## TODO: In the future, if I am interested in the nH, gamma, and Tin unc, I cannot fix them, as I do not get the uncertainty. It seemed like normalising the normalisation may fix the issue?
                    if mod_name == "powerlaw+diskbb":
                        # Get best-fit parameters and pass to initialise_model
                        parameters = {
                        "nh": f'{mod_obj.TBabs.nH.values[0]} -1',   # I'm fixing it because of issues with cflux in some cases
                        "gamma": f'{mod_obj.powerlaw.PhoIndex.values[0]} -1',  # I'm fixing it because of issues with cflux in some cases
                        "Tin": f'{mod_obj.diskbb.Tin.values[0]} -1', # I'm fixing it because of issues with cflux in some cases  
                        "norm1": f'{mod_obj.powerlaw.norm.values[0]} -1', # fix the first normalisation
                        "norm2": f'{mod_obj.diskbb.norm.values[0]} -1'
                        }
                        AllModels.clear()
                        mod, initial_pars = initialise_model("cflux_(powerlaw+diskbb)", flux_guess, fix_names, fix_values, parameters)
                        mod_obj = Model(mod)
                        mod_obj.setPars(initial_pars) # set the parameters to the initial parameters specified in the initialisation
                        if add_systematic: AllModels.systematic = 0.03 # add systematic error
                        mod_2 = True 
                        #AllModels.show() # for checking
                        Fit.delta = 1e-5
                        Fit.perform()

                    # Store the fit statistics
                    xrt_dict[mod_name]['chi2'].append(Fit.testStatistic) # test statistic value from the most recent fit
                    xrt_dict[mod_name]['dof'].append(Fit.dof) # the degrees of freedom from the fit
                    try: xrt_dict[mod_name]['redchi2'].append(Fit.testStatistic/Fit.dof)
                    except: xrt_dict[mod_name]['redchi2'].append(-1) # dof =0 
                    print("Fit succeeded.")
                
                except: # Fit performing failed
                    print("Fit failed.")
                    fit = False
            
            if not fit:
                xrt_dict[mod_name]['chi2'].append(-1) 
                xrt_dict[mod_name]['dof'].append(-1) 
                xrt_dict[mod_name]['redchi2'].append(-1)

            
            if fit: # Get the parameter errors

                try:
                    # Determine the confidence intervals of a fit (FitManager), i.e. get the parameter errors
                    # Note, the default is to estimate the 90% confidence ranges for all parameters. 
                    # If we put "1.0", this indicates that we want 68% uncertainty.
                    # The maximum keyword ensures that error will not be run if the reduced chi-squared of the best fit exceeds <redchi>. The default value for <redchi> is 2.0. We set it a bit higher.
                    n_params = mod_obj.nParameters
                    Fit.error(f'maximum 3.0 nonew 1.0 1-{n_params}') 
                    print("Fit error calculation succeeded.")

                    for comp_name in mod_obj.componentNames: # Iterate through components, e.g. diskbb
                        comp = getattr(mod_obj, comp_name)
                        for par_name in comp.parameterNames: # Iterate through parameters, e.g. Tin
                            if par_name not in ["Emin", "Emax", "eMin", "eMax"]: # Skip fixed parameters (like norm for diskbb, and Emin/Emax)
                                if par_name == 'norm' and comp_name == 'diskbb': pass ###TODO: Could remove?
                                else:
                                    param = getattr(comp, par_name)
                                    xrt_dict[mod_name].setdefault(par_name, []).append(param.values[0]) # best value
                                    # error[0] is lower bound and error[1] is upper bound
                                    if param.error[0]==0 and param.error[1]==0: neg_er, pos_er = 0, 0 # parameter was fixed during fitting
                                    else: neg_er, pos_er = abs(param.values[0] - param.error[0]), abs(param.error[1] - param.values[0]) 
                                    xrt_dict[mod_name].setdefault(par_name + '_neg', []).append(neg_er) 
                                    xrt_dict[mod_name].setdefault(par_name + '_pos', []).append(pos_er) 
                
                except: 
                    print("Fit error calculation failed.")
                    fit = False

            if not fit: # fit was not performed or failed
                for comp_name in mod_obj.componentNames:
                    comp = getattr(mod_obj, comp_name)
                    for par_name in comp.parameterNames:
                        if par_name not in ["Emin", "Emax", "eMin", "eMax"]:
                            if par_name == 'norm' and comp_name == 'diskbb': pass ###TODO: Could remove?
                            else: # append -1s to indicate failure
                                xrt_dict[mod_name].setdefault(par_name, []).append(-1)
                                xrt_dict[mod_name].setdefault(par_name + '_neg', []).append(-1)
                                xrt_dict[mod_name].setdefault(par_name + '_pos', []).append(-1)
                if mod_name=="powerlaw+diskbb" and mod_2==False: # For the 2-step diskbb+powerlaw, if the fit failed before the cflux model was initialised
                    xrt_dict[mod_name].setdefault('lg10Flux', []).append(-1)
                    xrt_dict[mod_name].setdefault('lg10Flux' + '_neg', []).append(-1)
                    xrt_dict[mod_name].setdefault('lg10Flux' + '_pos', []).append(-1)

            if fit: # Plot residuals from the last fit 
                plot_resid(spectrum_file[10:-8], mod_name, fixing)
            
            print()
        print("------------------------------------------------------------")

    # Save the results dictionary
    with open(folder+'xrt_spectral_dict.json', 'w') as j:
        json.dump(xrt_dict, j, indent = 4)      

    Xset.closeLog()
    
    


