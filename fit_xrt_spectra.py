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

from xspec import * 


def iso2mjd(iso_dates):
    # Convert ISO dates to MJD using astropy Time
    times = Time(iso_dates, format='isot', scale='utc')
    return times.mjd


## Initialise different spectral models, based on the input model_index.
##TODO: Normalisation when there is both diskbb and powerlaw
def initialise_model(mod_index, flux_guess, fix_names, fix_values):
    
    ## Initialise the parameter parameters
    # We can enter ranges for the parameter by passing in a string containing <param>,<delta>,<hard min>,<soft min>,<soft max>,<hard max> where:
    # param: trial param value used initially in fit
    # delta: step size used in the numerical determination of the derivatives during the fitting process. This value may be overriden for all parameters by the xset delta command option, which will apply a proportional rather than a fixed delta.
    # Or, if we specify a value and then -1, it is fixes to this value.

    # Emin and Emax are set to the energy range over which we want the flux to be calculated.
    Emin = 1.0
    Emax = 10.0

    # Defaults for starting parameters and ranges
    parameters = {
        "nh": '0.2,,0.1,0.1,10,10',  # absorption
        "gamma": '1.7,,0.0,0.0,4.0,4.0',  # spectral index
        "Tin": '1.0,,0.05,0.05,5.0,5.0'  # temperature
    }
    
    # Set parameters to the specified values
    for name, value in zip(fix_names, fix_values):
        if name in parameters: parameters[name] = f"{value} -1"
    nh = parameters["nh"]
    gamma = parameters["gamma"]
    Tin = parameters["Tin"]


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
    if mod_index == 0:
        mod = 'tbabs * pegpwrlw' 
        initial_pars = {1:nh, 
                        2: gamma, 3:Emin, 4:Emax, 5:f'{flux_guess}'}

    ## Absorbed disk blackbody
    # par1: nH; equivalent hydrogen column, in units of 10^{22} atoms/cm^2 
    # par2: Emin; Minimum energy over which flux is calculated.
    # par3: Emax; Maximum energy over which flux is calculated.
    # par4: lg10Flux; log (base 10) flux in erg/cm^2/s
    # par5: temperature of the inner disk radius (keV)
    # par6: norm; (Rin/ D10)^2 cos0, where where Rin is an apparent inner disk radius in km, D10 the distance to the source in units of 10 kpc, and the angle of the disk ( = 0 is face-on). On the correction factor between the apparent inner disk radius and the realistic radius, see e.g. Kubota et al. 1998.
    elif mod_index == 1:
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
    elif mod_index == 2:
        mod = 'tbabs * (pegpwrlw + cflux *  diskbb)'
        bb_flux = np.log10(flux_guess * 0.9 * 1e-12)
        pwrlw_flux = flux_guess*0.1
        initial_pars = {1:nh, 
                        2: gamma, 3:Emin, 4:Emax, 5:f'{pwrlw_flux}',
                        6:Emin, 7:Emax, 8:f'{bb_flux}', 9: Tin, 10:'1.0 -1'}

    else:
        sys.exit('Incorrect index: Please input either 0, 1, or 2')

    return mod, initial_pars # return the model and initial guesses


## Helper method to plot the residuals for the fits
def plot_resid(spectrum_name, mod_name):

    mpl.rcParams['xtick.labelbottom'] = True
    plots_dir = "./spectral_fit_residuals/"
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)

    Plot('data') # data from most recent Fit that was performed

    # Get coordinates from plot:
    x = Plot.x()
    rates = Plot.y()
    yer = Plot.yErr()
    folded = Plot.model()
    resids = np.array(rates) - np.array(folded)

    # Create the Matplotlib plot
    fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={'height_ratios': [3, 1]})

    # Plot using Matplotlib:
    ax[0].errorbar(x, rates, yerr=yer, fmt='ro',ms=1, label="Data", elinewidth=0.2)
    #ax[0].plot(x, rates, 'ro', label="data", ms=1)
    ax[0].plot(x, folded, 'b', label="model", linewidth=1)
    ax[0].set_ylabel(r'counts/cm$^2$/sec/keV')
    ax[0].legend(fontsize=11)

    ax[1].plot(x, resids, 'g', label="residuals (data-model)", linewidth=1)
    ax[1].set_xlabel('Energy [keV]')
    ax[1].legend(fontsize=11)

    ax[0].set_title("Spectrum: " + spectrum_name + " & model: " + mod_name)

    plt.savefig('./spectral_fit_residuals/'+ spectrum_name+"_"+mod_name)

    ax[0].set_yscale("log")
    plt.savefig('./spectral_fit_residuals/'+ spectrum_name+"_"+mod_name+"_log")



##TODO: Add input for the normalisation of bb to powerlaw
# The parameter fix specifies which parameters need to be fixed during the fit and for which spectral files
def run_spectral_fit(fix={}):
    
    # Initialise XSPEC
    # Xset: storage class for XSPEC settings
    Xset.abund = "wilm"    # set the abundance table used in the plasma emission and photoelectric absorption models; 'wilm' is an XSPEC built-in abundance table
    Xset.xsect = "vern"    # set the photoelectric absorption cross-sections
    Xset.chatter = 0 # set the console chatter level... change this to 10 for more chatter
    Xset.logChatter = 20
    Xset.openLog("xspec.log") 
    Xset.parallel.error = 10 # use up to 10 parallel processes during the Fit.error() command runs

    # Fit: manager class for setting properties and running a fit
    Fit.query = "yes"      # when nIterations is reached, continue the fit without stopping to query
    Fit.nIterations = 100  # fit iterations for each attemp
    Fit.statTest = "chi"

    # Plot: manager class for performing XSPEC plots
    Plot.xAxis = "KeV"     # x-axis for plotting set to energy instead of channel
    Plot.device = '/null'
    Plot.yLog = True

    # Load in the spectral files
    spectral_files = glob.glob('./spectra/*final.pi') # the spectra after grouping script has been run
    n_spectra = len(spectral_files) 
    
    # Initialise and fill variables
    obs_isots  = [] # to hold the observation dates
    bin_counts = [] # to hold number of counts in grouped bin -- in order to decide whether chi^2 or C stats should be used
    for spectrum_file in spectral_files:
        header = fits.getheader(spectrum_file, ext = 1)
        obs_isots.append(header['DATE-OBS'])
        bin_counts.append(header['COUNTGRP'])

    # Convert dates to MJDs
    obs_mjds = iso2mjd(np.array(obs_isots))    
    
    # Re-order lists by time
    sort_index = np.argsort(obs_isots)
    obs_isots = np.array(obs_isots)[sort_index].tolist()
    spectral_files = np.array(spectral_files)[sort_index].tolist()
    bin_counts = np.array(bin_counts)[sort_index].tolist()
    obs_mjds = obs_mjds[sort_index].tolist()

    xrt_dict = {model: {"obs_isot": obs_isots, "obs_mjds": obs_mjds, "cstat?":[], "chi2/cstat": [], "dof": [], "redchi2/redcstat": []} for model in ['powerlaw', 'diskbb', 'both']}


    ## Iterate through each spectrum (i.e. each observation), getting the fit parameters
    # AllData: container for all loaded data sets (objects of class Spectrum).
    # AllModels: container for all Model objects.
    # XSCPEC models are defined by creating a Model object.
    ## Note: getattr(object, 'x') is completely equivalent to object.x, but used in the case we don't know exactly what 'x' is.
    for k, spectrum_file in enumerate(spectral_files):

        print('\n\n', spectrum_file)

        fix_names = []
        fix_values = []
        for param, details in fix.items():
            indices = details["indices"]
            value = details["value"]
            if indices == "all" or k in indices or k in np.array(indices)+n_spectra: # the last check is in case some of the values were negative
                fix_names.append(param)
                fix_values.append(value)

        # Load the spectrum into the AllData object, removing any previously loaded data sets. Could also run AllData.clear() before just to be sure.
        # If the rmf and arf files are located in the same folder, these will automatically be loaded.
        AllData(spectrum_file)
        # AllData.show() # check current state of the AllData container

        # Define range of energy to use -- i.e. ignore certain channels.
        # The floating point values are assumed to be in keV, as this is what was set for the Plot.xAxis value above.
        # This is done for all the loaded spectra.
        AllData.notice('*:0.5-10.0') # so that edge points are used
        AllData.ignore('*:10.0-**')
        AllData.ignore('*:**-0.5')
        AllData.ignore('bad') # ignore data that is bad -- using the quality column

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


        ## Next, we will get an estimate of the flux (flux_guess) using a simple powerlaw fit
        # Define the spectral model -- absorbed power law
        # 1: nH (10^{22} atoms/cm^2); 2: photon index; 3: lower peg energy range (keV); 4: upper peg energy range (keV); 5: norm i.e. flux (10^{-12} ergs/cm^2/s)
        mod1 = Model('tbabs * pegpwrlw')
        mod1.setPars({1:'0.2 -1', 2:'2.0 -1', 3:0.5, 4:10.0, 5:'1000.0'})
        Fit.delta = 1e-2 # controls the convergence criterion for the fitting algorithm
        Fit.perform()
        flux_guess = mod1.pegpwrlw.norm.values[0] # flux (in units of 10^{-12} ergs/cm^2/s)


        ## Iterate through the model versions, fitting each type and appending the results to an output dictionary
        for mod_index, mod_name in enumerate(['powerlaw', 'diskbb', 'both']):

            print('\n', mod_name)

            if bin_counts[k] == 1: xrt_dict[mod_name]['cstat?'].append(True)
            else: xrt_dict[mod_name]['cstat?'].append(False)

            AllModels.clear() # AllModels represents the collection of all currently defined models in the XSPEC session

            # Initialise the model and systematics 
            mod, initial_pars = initialise_model(mod_index, flux_guess, fix_names, fix_values)
            mod_obj = Model(mod)
            mod_obj.setPars(initial_pars) # set the parameters to the initial parameters specified in the initialisation
            AllModels.systematic = 0.03 # add systematic error
            AllModels.show()

            try:
                # Adjust the fit precision, and perform the fit of the spectrum. The factor will multiply the parameter value to produce a fit delta.
                # The first pass allows the fit to reach a close approximation of the best-fit parameters quickly. 
                # Each subsequent pass with a smaller delta refines the fit further.
                for delt in [1e-2, 1e-3, 1e-4, 1e-5]:
                    Fit.delta = delt
                    Fit.perform()

                # Store the fit statistics
                xrt_dict[mod_name]['chi2/cstat'].append(Fit.testStatistic) # test statistic value from the most recent fit
                xrt_dict[mod_name]['dof'].append(Fit.dof) # the degrees of freedom from the fit
                xrt_dict[mod_name]['redchi2/redcstat'].append(Fit.testStatistic/Fit.dof)
            
            except: # fit performing failed
                xrt_dict[mod_name]['chi2/cstat'].append(-1) # test statistic value from the most recent fit
                xrt_dict[mod_name]['dof'].append(-1) # the degrees of freedom from the fit
                xrt_dict[mod_name]['redchi2/redcstat'].append(-1)


            try:
                # Determine the confidence intervals of a fit (FitManager), i.e. get the parameter errors
                # Estimate the 90% confidence ranges for all parameters
                n_params = mod_obj.nParameters
                Fit.error(f'1-{n_params}')
                print("Fit succeeded.")

                for comp_name in mod_obj.componentNames: # Iterate through components
                    comp = getattr(mod_obj, comp_name)
                    for par_name in comp.parameterNames: # Iterate through parameters
                        if par_name not in ["Emin", "Emax", "eMin", "eMax"]: # Skip fixed parameters (like norm for diskbb, and Emin/Emax)
                            if par_name == 'norm' and comp_name == 'diskbb':
                                pass
                            else:
                                param = getattr(comp, par_name)
                                xrt_dict[mod_name].setdefault(par_name, []).append(param.values[0]) # best value
                                if param.error[0]==0 and param.error[1]==0: neg_er, pos_er = 0, 0 # parameter was fixed during fitting
                                else: neg_er, pos_er = abs(param.values[0] - param.error[0]), abs(param.error[1] - param.values[0]) # ... error[0] is lower bound and error[1] is upper bound
                                xrt_dict[mod_name].setdefault(par_name + '_neg', []).append(neg_er) 
                                xrt_dict[mod_name].setdefault(par_name + '_pos', []).append(pos_er) 
            
            except: # fit failed
                print("Fit failed.")
                for comp_name in mod_obj.componentNames:
                    comp = getattr(mod_obj, comp_name)
                    for par_name in comp.parameterNames:
                        if par_name not in ["Emin", "Emax", "eMin", "eMax"]:
                            if par_name == 'norm' and comp_name == 'diskbb':
                                pass
                            else: # append -1s to indicate failure
                                xrt_dict[mod_name].setdefault(par_name, []).append(-1)
                                xrt_dict[mod_name].setdefault(par_name + '_neg', []).append(-1)
                                xrt_dict[mod_name].setdefault(par_name + '_pos', []).append(-1)
                
            plot_resid(spectrum_file[10:-8], mod_name)
            
            print()
        print("------------------------------------------------------------")

    # Save the final XRT dictionary
    if not os.path.exists('./spectral_fit_results/'):
        os.makedirs('./spectral_fit_results/')

    with open('./spectral_fit_results/xrt_spectral_dict.json', 'w') as j:
        json.dump(xrt_dict, j, indent = 4)      

    Xset.closeLog()





