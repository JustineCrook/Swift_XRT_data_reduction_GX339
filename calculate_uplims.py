import numpy as np
import scipy.stats 
from astropy.io import fits
import gzip
from astropy.io import fits
import re
import os
import numpy as np
import subprocess
from xspec import * 


"""
For some data points, we only have upper limits. 
We therefore cannot use spectral fitting to obtain the flux, so we need to rely on alternative methods. 
"""


##TODO: Changing model/parameters for this fit
##TODO: Functionality for when it is in the soft state, and the diskbb is therefore what needs to be used.


# Note: The upper limits will always be when the mode is PC (instead of WC), as this is the mode that is used when the lowest count rates are observed.

##########################################################################################################################################################################################################################


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



## Make .reg files for the source and background
def generate_region_files(ra, dec, count_rate):

    # From Adams 2009 Table 1
    if count_rate>0.5: npix =30
    elif (count_rate>0.1 and count_rate<=0.5): npix=25
    elif (count_rate>0.05 and count_rate<=0.1): npix=20 
    elif (count_rate>0.01 and count_rate<=0.05): npix=15 
    elif (count_rate>0.005 and count_rate<=0.01): npix=12
    elif (count_rate>0.001 and count_rate<=0.005): npix=9
    elif (count_rate>0.0005 and count_rate<=0.001): npix=7
    elif count_rate<=0.0005: npix=5   
    else: print("Invalid count rate.")

    print("Source npix: ", npix)

    f = 6.548089e-04  # degrees...  size of a Swift XRT pixel

    # Calculate size in arcseconds
    size_asec = np.round(int(npix) * f * 3600, 3)  # Convert to arcseconds and round to 3 decimal places

    # Write to bkg.reg
    with open("./uplims_analysis/bkg.reg", "w") as bkg_file:
        bkg_file.write("icrs\n")
        bkg_file.write(f"annulus({ra},{dec},141.439\",259.304\")\n")

    # Write to src_{npix}.reg
    src_filename = f"./uplims_analysis/src_{npix}pix.reg"
    with open(src_filename, "w") as src_file:
        src_file.write("icrs\n")
        src_file.write(f"circle({ra},{dec},{size_asec}\")\n")

    print(f"Files 'bkg.reg' and '{src_filename}' have been created with the specified parameters.")

    return npix


def run_command(command, input_data=None):
    """Run a shell command with optional input and return output."""
    result = subprocess.run(command, input=input_data, text=True, shell=True, capture_output=True)
    return result.stdout


def extract_curve_info(type, obs_id, npix=None):
    """Extracts n_src or n_bkg value after running 'extract curve'."""

    if type=="src": name= f"./uplims_analysis/src_{npix}pix"
    elif type=="bkg": name="./uplims_analysis/bkg"

    xselect_cmds = f"""
    xselect << EOF
xsel
no
read events ./uplims_analysis/sw{obs_id}xpcw3po_cl.evt
./
yes
filter region {name}.reg
filter pha_cutoff 100 1000
extract curve
EOF
    """
    output = run_command(xselect_cmds)
    index=0
    lines =output.splitlines()
    for i, element in enumerate(lines):
        if element.startswith("    Grand"):
            index=i
    n_value = lines[index+1].split()[1]
    return n_value


def extract_backscal_exposure(type, obs_id):
    """Extracts BACKSCAL and EXPOSURE from a PHA file header."""

    if type=="src": filename= f"./uplims_analysis/src_{obs_id}.pha"
    elif type=="bkg": filename=f"./uplims_analysis/bkg_{obs_id}.pha"

    with fits.open(filename) as hdul:
        header = hdul[1].header
        backscal = header.get("BACKSCAL")
        exposure = header.get("EXPOSURE")
    return backscal, exposure


def extract_spectrum(type, obs_id, npix=None):

    if type=="src": name= f"./uplims_analysis/src_{npix}pix"
    elif type=="bkg": name="./uplims_analysis/bkg"

    xselect_spectrum_cmds = f"""
    xselect << EOF
xsel
no
read events ./uplims_analysis/sw{obs_id}xpcw3po_cl.evt
./
yes
filter region {name}.reg
extract spectrum
save spectrum ./uplims_analysis/{type}_{obs_id}.pha
EOF
    """

    output = run_command(xselect_spectrum_cmds)
    


## TODO:
# Checking output of xrtmkarf to get location rmf file to use
def run_xspec(obs_id, nH):
    """Run XSPEC with the specified model and get count_rate_factor."""

    AllData.clear()
    AllModels.clear()

    # Set the abundance table used in the plasma emission and photoelectric absorption models; 'wilm' is an XSPEC built-in abundance table
    Xset.abund = "wilm"    
    Xset.xsect = "vern"  

    Xset.openLog("xspec.log")

    spectrum_file = f"./uplims_analysis/src_{obs_id}.pha"
    AllData(spectrum_file) # Load the spectrum into the AllData object
    AllData(1).response = "/usr/local/shared/caldb/data/swift/xrt/cpf/rmf/swxpc0to12s6_20210101v016.rmf"
    AllData(1).response.arf = f"./uplims_analysis/src_{obs_id}.arf"
    AllData.ignore("**-1.0,10.0-**")
    AllData.notice("1.0-10.0")
    AllData.ignore("bad")
    
    # Set up the model
    Model("tbabs*pegpwrlw", setPars={1:str(nH)+" -1",2:"1.7 -1",3:1.0, 4:10.0, 5:1.0})

    print("---------------DATA:---------------")
    AllData(1).show()

    # Show rates for the model applied to the spectrum
    count_rate_factor = AllData(1).rate[3] # this is the model predicted rate

    Xset.closeLog()
    
    return count_rate_factor


def get_values_for_uplims(obs_id, ra, dec, count_rate, nH): 
    
    # Delete pre-existing results
    file_paths = [f"./uplims_analysis/bkg_{obs_id}.pha", f"./uplims_analysis/src_{obs_id}.pha", f"./uplims_analysis/src_{obs_id}.arf"]
    for file_path in file_paths:
        if os.path.exists(file_path): os.remove(file_path)

    # Create the region files for the source and background
    npix = generate_region_files(ra, dec, count_rate)

    ## Source:
    # Get the number of good counts
    n_src = extract_curve_info("src", obs_id, npix) 
    # Create the .pha file
    extract_spectrum("src", obs_id, npix) 
    # Run xrtmkarf to generate ARF file
    xrtmkarf_cmd = f"xrtmkarf phafile=./uplims_analysis/src_{obs_id}.pha outfile=./uplims_analysis/src_{obs_id}.arf expofile=./uplims_analysis/sw{obs_id}xpcw3po_ex.img srcx=-1 srcy=-1 psfflag=yes"
    run_command(xrtmkarf_cmd)
    # Get backscal_src and exposure
    backscal_src, exposure = extract_backscal_exposure("src", obs_id)
    
    ## Background:
    # Get the number of good counts
    n_bkg = extract_curve_info("bkg", obs_id)
    # Create the .pha file
    extract_spectrum("bkg", obs_id)
    # Get backscal_bkg
    backscal_bkg, _ = extract_backscal_exposure("bkg", obs_id)
    
    ## Run XSPEC to get count_rate_factor:
    count_rate_factor = run_xspec(obs_id, nH)
    
    ## Print results:
    print()
    print(f"RESULTS for {obs_id}:")
    print(f"n_src: {n_src}")
    print(f"backscal_src: {backscal_src}")
    print(f"exposure: {exposure}")
    print(f"n_bkg: {n_bkg}")
    print(f"backscal_bkg: {backscal_bkg}")
    print(f"count_rate_factor: {count_rate_factor}")
    print()

    return int(n_src), float(backscal_src), float(exposure), int(n_bkg), float(backscal_bkg), float(count_rate_factor)



##########################################################################################################################################################################################################################


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
def get_uplim_fluxes(n_src= [2, 4], n_bkg=[29, 64], backscal_src=[8e-5,8e-5], backscal_bkg=[2.67e-2, 2.67e-2], exposure=[1.952880358167567E+03, 7.566808800853420E+02], count_rate_factor = [1.536e-03, 5.286e-03]):

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

    flux_conversion = 1e-12/np.array(count_rate_factor)

    # Convert the count rate to flux units
    flux_uplims = flux_conversion * count_rate

    print('3-sigma upper limit on the flux is: ', flux_uplims)

    return flux_uplims


##########################################################################################################################################################################################################################


##TODO: Try put xrtpipeline_runner.sh file into this python file?
# I can set get_data=False, when I am re-running this and the data has already been loaded
def uplims_runner(obs_id_ar, count_rate_ar, target_coords, nH, get_data=True):

    ra, dec = str(target_coords[0]), str(target_coords[1])

    if get_data:

        sh_file = "xrtpipeline_runner.sh"
        if not os.path.exists("./uplims_analysis"): os.makedirs("./uplims_analysis")
        if not os.path.exists("./uplims_analysis/data"): os.makedirs("./uplims_analysis/data")
        if not os.path.exists("./uplims_analysis/xrtpipeline_output"): os.makedirs("./uplims_analysis/xrtpipeline_output")
        
        ## Get the raw event files for the observations of interest: https://www.swift.ac.uk/swift_portal/. To do this, use the wget option in the /data sub-folder.
        ## Then, in /uplims_analysis, run xrtpipeline for the observations of interest.
        ## Copy the required data to /uplims_analysis -- the event files and exposure files.
        try: # Run the shell script with RA, DEC, and obs_id_ar as arguments
            print("Running bash script...")
            result = subprocess.run(['bash', sh_file, ra, dec] + obs_id_ar, check=True, text=True, capture_output=True)
            #print("Output:\n", result.stdout)
        except subprocess.CalledProcessError as e:
            print("Error:\n", e.stderr)
  
    ## Get the values needed to calculate the 3-sigma flux upper limits, using the event and exposure files
    print("Getting values...")
    all_n_src, all_backscal_src, all_exposure, all_n_bkg, all_backscal_bkg, all_count_rate_factor = [], [], [], [], [], []
    for obs_id, count_rate in zip(obs_id_ar, count_rate_ar): 
      n_src, backscal_src, exposure, n_bkg, backscal_bkg, count_rate_factor = get_values_for_uplims(obs_id, ra, dec, count_rate, nH)
      all_n_src.append(n_src), all_backscal_src.append(backscal_src), all_exposure.append(exposure), all_n_bkg.append(n_bkg), all_backscal_bkg.append(backscal_bkg), all_count_rate_factor.append(count_rate_factor)

    ## Calculate the 3-sigma flux upper limits.
    print("Calculating fluxes...")
    flux_uplims = get_uplim_fluxes(all_n_src, all_n_bkg, all_backscal_src, all_backscal_bkg, all_exposure, all_count_rate_factor)

    return flux_uplims
    