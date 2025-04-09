This is a pipeline to do spectral fitting of Swift XRT data to extract the 1-10 keV fluxes. Note that this pipeline uses simple spectral models, so is not appropriate if a more thorough spectral fit is required. However, it is sufficient if we would just like to extract fluxes (for example, for compiling the radio:X-ray plane). 

The initial code was obtained from Andrew Hughes. 

All parameters that need changing in order to run the pipeline are in input_parameters.py. 

To run the code, you will need swifttools and PyXspec. For the upper limit analysis, HeaSoft and the Swift XRT CALDB are required. 


The strategy is as follows:


# STEP 1: 
- Use `get_data` to get the light curve data and spectra using the [swifttools.ukssdc.xrt_prods](https://www.swift.ac.uk/API/ukssdc/xrt_prods.md) module.

For this, we require the following information:
- The source and observation IDs, names, and source position. These can be found [here](https://www.swift.psu.edu/operations/obsSchedule.php). The naming convention is 000 + targetID + segment number.


# STEP 2: 
- Use `group_spectra` to group the spectral data such that there are `min_counts_chi` (recommended is either 20 or 25) counts per bin if the total counts are more than or equal to 300. Else, set group to 1 count per bin.
- Also specify the lower energy bound for the fitting (which must be less than or equal to 1.0).


# STEP 3:
- Use `plot_lightcurve_and_hr` to plot the light curve and hardness ratio, as output by the swifttools.ukssdc pipeline.
- The code also outputs which observations are found to be upper limits by the swifttools.ukssdc pipeline. We will not fit these spectral files, and rather obtain the flux upper limits by the methods described in STEP 7. Use the website to get the IDs for these files. These will then be used as input to the following spectral fitting routines, to ensure that they are not fitted. 


# STEP 4:
- Run the initial spectral model fitting, using `unconstrained_fit`. 
- Specify the models to be fit.
- In this routine, we do not fit the observations identified as upper limits in STEP 3.
- We also choose to ignore PC spectra with fewer than 3 counts, and WT spectra with fewer than 15 counts because these cannot reliably be fit. 
- Lastly, we do not fit observations with counts less than `count_threshold`, which is set to 50 by default. These will need to have their spectral parameters fixed (except the flux) before fitting -- see STEP 6.
- All the models are fitted to all the remaining observations. 
- For the two-component model consisting of a powerlaw and disk blackbody, we have the option to do a two-step fit, rather than using `cflux` from the start. We have found that the former is more robust, compared to the latter where the model sometimes converges to local minima.
- The errors are calculated using the Xspec `error` command. This is only conducted if the chi^2/dof is less than or equal to 3.0. 
- The results are output to a folder called spectral_fit_results. The fit parameters are shown in fit_outputs.txt. If the fitting failed or the chi^2/dof was too high, '-1' is output.
- The residuals are output to a folder called spectral_fit_residuals.


# STEP 5:
- Using the previous results and any prior knowledge, decide on which observations should be fit with which models. The observations are indexed as per the file spectral_fit_results/fit_outputs.txt
- Input these indexes to a function called `get_param_averages`. This calculates the average and weighted averages of the parameter nH -- and gamma and Tin, if relevant.
- The parameters as a function of time, and distribution of the parameters are shown in the spectral_fit_results folder. 


# STEP 6:
- Re-run spectral fitting using `constrained_fit`, but this time fitting the parameters for some observations -- using the results of STEP 5 and/or any prior knowledge (e.g. nH values in papers).
- The results are shown in the folders spectral_fit_results_fixing and spectral_fit_residuals_fixing.


# STEP 7:
- Using `get_uplims`, get the upper limit flux values of the observations previously identified as upper limits. 
- For this, an nH and gamma value need to be supplied.  

The steps for this are as follows:
- Use the 3-sigma upper limit count rate supplied by the swifttools.ukssdc pipeline (STEP 3) as an initial count rate estimate.
- Create extraction regions for the source and background following [Evans et al. 2009](https://ui.adsabs.harvard.edu/abs/2009MNRAS.397.1177E/abstract).
- Using the extraction regions, get the number of total counts for the source and background by using `xselect` and extracting a curve/spectra. 
- Record the BACKSCAL for the source and background. Use the ratio of the BACKSCAL parameters to correct for the different area of the extraction regions.
- Adopt 3 times the upper Gehrels error (divided by the exposure time) as your count rate.
- To convert from a count rate to a flux, first make a spectrum, arf, and rmf file. 
- Load in the spectrum and the arf (using `arf [ARFFILE]`) and rmf (using `resp [RMFFILE]`) into Xspec. 
- Then pick a simple `tbabs*pegpwrlw` model (with the known nH), and predict a count rate using `show rates`. 
- Record the model-predicted rate for the assumed model flux (by default, just use 1e-12), and convert the count rate to flux.


Note that we re-calculate the 3-sigma count rate for the non-detections instead of using the ones in the pipeline because it was found that the region-size in the swifttools.ukssdc pipeline algorithm for non-detections is sometimes wrong, so we instead wanted more control over that. We could also have used WebPIMMS to convert count rates to fluxes, however we choose not to do this because we think the ARF/RMF files may cause very slight variations in time that WebPIMMS does not account for but doing it in HEASoft would.



# STEP 8:
- Plot and tabulate the final results, using `get_final_results`.  
- The final results and residuals and shown in final_results and final_residuals, respectively.



# TODO: 
- Clean up and check all code. 
- Especially check the upper limits routine.
- Scrape of observation IDs from the observation schedule website
- Scrape target coordinates and names from the observation schedule website
- Staggered data retrieval, to avoid overload
- Get the IDs for the upper limits automatically, using website.
- Minimise the need for user input by calling from other functions to fetch data e.g. Make it such that the fitting functions fetch the upper limit IDs, so that it does not require this user input.