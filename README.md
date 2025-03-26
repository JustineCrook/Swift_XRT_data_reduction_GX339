This is a pipeline to do spectral fitting of Swift XRT data to extract the 1-10 keV fluxes. Note that this pipeline uses simple spectral models, so is not appropriate if a more thorough spectral fit is required.

The initial code was obtained from Andrew Hughes. 

All parameters that need changing in order to run the pipeline are in input_parameters.py. 

To run the code, you will need swifttools and pyxspec. For the upper limit analysis, the heasoft and caldb are required. 


The strategy is as follows:


# STEP 1: 
- Use `get_data` to get the light curve data and spectra from swifttools.ukssdc

For this, we require the following information:
- The source and observation IDs, and source position. These can be found [here](https://www.swift.psu.edu/operations/obsSchedule.php). The naming convention is 000 + targetID + segment number.


# STEP 2: 
- Use `group_spectra` to group the spectral data such that there are `min_counts_chi` (recommended is either 20 or 25) counts per bin if the total counts are more than or equal to 300. Else, set group to 1 count per bin.
- Also specify the lower energy bound for the fitting (which must be less than or equal to 1.0).


# STEP 3:
- Use `plot_lightcurve_and_hr` to plot the light curve and hardness ratio
- The code outputs which observations are found to be upper limits by the swifttools.ukssdc software.
- We will not fit these spectral files, and rather obtain the flux upper limits by the methods described in STEP 7.
- Use the website to get the IDs for these files. These will then be used as input to the following spectral fitting routine, to make sure they are not fitted. 


# STEP 4:
- Run the initial spectral model fitting, using `unconstrained_fit`. 
- Specify the models to be fit.
- In this routine, we do not fit the observations identified as upper limits in step 3.
- We also choose not to fit PC spectra with fewer than 3 counts, and WT spectra with fewer than 15 counts. 
- All the models are fitted to all the remaining observations. 
- For the two-component model consisting of a powerlaw and disk blackbody, we have the option to do a two-step fit, rather than using cflux from the start. The former is more robust, compared to the latter where the model sometimes converges to local minima.
- The results are output to a json file. If the fitting failed, '-1' is output.
- The errors are calculated using the xspec `error` command. This is only conducted if the chi^2/dof is less than or equal to 3.0. If the errors are not calculated as a result, the values are filled with '-1's.


# STEP 5:
- Decide on which observations should be fit with which models. The observations are indexed as per the file fit_outputs.txt
- Input these indexes to a function called `get_param_averages`. This calculates the average and weighted averages of the parameter nh -- and gamma and Tin, if relevant.
- The parameters as a function of time, and distribution of the parameters are shown. 


# STEP 6:
- Re-run spectral fitting using `constrained_fit`, but this time fitting the parameters for some observations -- using the results of step 5 and/or any prior knowledge (e.g. nH values in papers).
- The results are shown in folders ending with '_fixing'


# STEP 7:
- Using `get_uplims`, get the upper limit flux values of the observations previously identified as upper limits. 
- For this, an nh and gamma value needs to be supplied.  

The steps for this are as follows:
- Use the 3-sigma upper limit count rate supplied by the swifttools.ukssdc pipeline (step 1) as an initial count rate estimate.
- Create extraction regions for the source and background following Evans 2009.
- Using the extraction regions, get the number of total counts for the source and background by using `xselect' and extracting a curve/spectra 
- Record the BACKSCAL for the source and background. Use the ratio of the BACKSCAL parameters to correct for the different area of the extraction regions.
- Then, once again calculate a net number of counts (and a count rate using the exposure time). 
- Adopt 3x the upper Gehrels error as your count rate.
- To convert from a count rate to a flux, first make a spectrum, arf, and rmf file. 
- Load in the spectrum and the arf (using arf [ARFFILE]) and rmf (using resp [RMFFILE]) into xspec. 
- Then pick a simple tbabs*pegpwrlw model (with the known NH), and predict a count rate using show rates. 
- Record the model-predicted rate for the assumed model flux (by default just use 1e-12), and convert the count rate to flux.


# STEP 8:
- Plot and tabulate the final results, using `get_final_results`.  