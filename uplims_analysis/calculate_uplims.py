import numpy as np
import scipy.stats 
from astropy.io import fits
import gzip

"""
For some data points, we only have upper limits. 
We therefore cannot use spectral fitting to obtain the flux, so we need to rely on alternative methods. 

Run uplims_script.py, and then use the output to run get_uplims_fluxes. The results can then be used in the final plots.
"""


## Calculate a 3-sigma upper limit on the flux of a source using a Poisson distribution. 
# Compute the Poisson confidence interval limits (lambda_l and lambda_u) for a count n with a given confidence interval (CI).
# n: The observed count rate, which should be an array.
# CI: confidence interval, with a default value that is roughly 68% for 1-sigma.
def gehrels_poisson(n, CI=scipy.special.erf(1/np.sqrt(2))):
  if n.any() < 0:
      print("n ({}) cannot be less than 0.".format(n))
      return()
  elif n.any() > 0:
    CL = (1+CI)/2 # confidence level adjusted for two-sided intervals
    lambda_l = scipy.special.gammaincinv(n,1-CL) # lower confidence limit
  else: # all n values are 0
    CL = CI
    lambda_l = 0

  lambda_u = scipy.special.gammaincinv(n+1,CL) # upper confidence limit

  return(lambda_l, lambda_u)




## Method to read "EXPOSURE" in header
def read_exposure_from_evt(file_path):
    """
    Reads the EXPOSURE value from the header of a .evt.gz FITS file.

    Parameters:
    file_path (str): Path to the .evt.gz file.

    Returns:
    float: The value of the EXPOSURE keyword, or None if not found.
    """
    # Open the gzip file, then use fits.open to read it as a FITS file
    with gzip.open(file_path, 'rb') as gz_file:
        with fits.open(gz_file) as hdul:
            # Check the primary header (extension 0) or other extensions as needed
            for hdu in hdul:
                if 'EXPOSURE' in hdu.header:
                    return hdu.header['EXPOSURE']
    return None  # Return None if EXPOSURE keyword is not found



## Calculate a 3-sigma upper limit on the flux of a source using a Poisson distribution. 
## Input (all should be arrays):
# n_src: source counts in a circle with pre-determined radius
# n_bkd: background counts in an annulus with pre-determined radii
# backscal_src and backscal_bkg are backscaling factors, used to account for differences in source and background area size
# exposure is the exposure time in seconds
# count_rate_factors is a part of the scaling factor between count rate and flux
## To determine the count_rate_factors:
# Fix the flux at 1e-12, Gamma at 1.7, and the known NH. 
# Then (don't fit), use the Xspec command "show rates" to get the denominator. 
# Therefore the conversion is 1e-12/[answer]
def get_uplim_fluxes(n_src= [2, 4], n_bkg=[29, 64], backscal_src=[8e-5,8e-5], backscal_bkg=[2.67e-2, 2.67e-2], exposure=[1.952880358167567E+03, 7.566808800853420E+02], count_rate_factors = [1.536e-03, 5.286e-03]):

    # Use the appropriate source region size from Evans et al. 2009
    # xselect on the src (or background) region, and it will give the correct number of 'good' counts
    n_src = np.array(n_src) # circle
    n_bkg = np.array(n_bkg) # annulus
    
    # Read these parameters (i.e., backscal and exposure from header files)
    backscal_src= np.array(backscal_src) 
    backscal_bkg = np.array(backscal_bkg) 
    backscal_factor = backscal_bkg/backscal_src

    exposure = np.array(exposure)

    # Error calculations
    src_err = (np.array(gehrels_poisson(n_src))-n_src)[1]
    bkg_err = (np.array(gehrels_poisson(n_bkg))-n_bkg)[1]

    # Calculate the background counts and error on the background counts, adjusted by backscal_factor
    bkg_in_aperture     = n_bkg/backscal_factor 
    bkg_in_aperture_err = bkg_err/backscal_factor 

    # Calculate net source counts (source counts minus background counts in aperture).
    net = n_src - bkg_in_aperture
    net_err = np.sqrt(src_err*src_err+bkg_in_aperture_err*bkg_in_aperture_err) # net error -- just add errors in quadrature, as a simple scheme

    # Scale the net error for a 3-sigma upper limit
    count_rate = 3 * net_err / exposure

    flux_conversion = 1e-12/np.array(count_rate_factors)

    # Convert the count rate to flux units
    flux_uplims = flux_conversion * count_rate

    print('3-sigma upper limit on the flux is: ', flux_uplims)

    return flux_uplims


