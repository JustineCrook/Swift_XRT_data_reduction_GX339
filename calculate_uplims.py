import numpy as np
import scipy.stats 
from astropy.io import fits
import gzip
import re
import os
import numpy as np
import subprocess
from xspec import * 

"""
For some data points, we only have upper limits. 
We therefore cannot use spectral fitting to obtain the flux, so we need to rely on alternative methods. 
Before running this file, we need: 
module load heasoft/6.33
source /mnt/users/crookmansourj/swiftenv/bin/activate
"""


## NOTE
# We may need to look at the background region manually to check that it does not extend beyond the edge of the detector window. 
# For the WT upper limits, the strategy to reset the BACKSCAL values is only appropriate if the source is located within the WT window and not near the bad columns.


##########################################################################################################################################################################################################################


print("Current working directory:", os.getcwd())  # Print the current directory for debugging

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
# For the source region, we use the recommendations of Adams 2009 Table 1
# Note that the Adams table refers to uncorrected count rates whereas count_rate is the corrected value, but it it makes extremely little difference in the low-count regime considered here.
# For the background region, we also use the recommendation of Adams 2009 -- an annulus with an inner radius of 60 pixels (142 arcsec) and an outer radius of 110 pixels (260 arcsec)
def generate_region_files_pc(ra, dec, count_rate):

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

    # Write to bkg_pc.reg
    bkg_filename = "./uplims_analysis/bkg_pc.reg"
    if not os.path.exists(bkg_filename):
        with open(bkg_filename, "w") as bkg_file:
            bkg_file.write("icrs\n")
            bkg_file.write(f"annulus({ra},{dec},141.439\",259.304\")\n")
        print(f"File '{bkg_filename}' has been created.")

    # Write to src_pc_{npix}.reg
    src_filename = f"./uplims_analysis/src_pc_{npix}pix.reg"
    if not os.path.exists(src_filename):
        with open(src_filename, "w") as src_file:
            src_file.write("icrs\n")
            src_file.write(f"circle({ra},{dec},{size_asec}\")\n")
        print(f"File '{src_filename}' have been created with the specified parameters.")

    return npix

##TODO: May want to make this a box instead, but then need to be careful about the ARF file
# 20-pixel circle centred on the source position for the source file 
# 80-120 annulus for the background
# This should work okay if the source is quite well if the source is more or less in the centre of the window.
def generate_region_files_wt(ra, dec):
    
    # Write to bkg_pc.reg
    bkg_filename = "./uplims_analysis/bkg_wt.reg"
    if not os.path.exists(bkg_filename):
        with open(bkg_filename, "w") as bkg_file:
            bkg_file.write("icrs\n")
            bkg_file.write(f"annulus({ra},{dec},188.585\",282.877\")\n")
        print(f"File '{bkg_filename}' has been created.")

    npix = 20 
    print("Source npix: ", npix)

    f = 6.548089e-04  # degrees...  size of a Swift XRT pixel

    # Calculate size in arcseconds
    size_asec = np.round(int(npix) * f * 3600, 3)  # Convert to arcseconds and round to 3 decimal places

    # Write to src_pc_{npix}.reg
    src_filename = f"./uplims_analysis/src_wt_{npix}pix.reg"
    if not os.path.exists(src_filename):
        with open(src_filename, "w") as src_file:
            src_file.write("icrs\n")
            src_file.write(f"circle({ra},{dec},{size_asec}\")\n")
        print(f"File '{src_filename}' have been created with the specified parameters.")

    return npix


def run_command(command, input_data=None):
    """Run a shell command with optional input and return output."""
    result = subprocess.run(command, input=input_data, text=True, shell=True, capture_output=True)
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
        raise RuntimeError(f"Command failed: {command}")
    return result.stdout



def extract_curve_info(type, obs_id, npix=None):
    """Extracts n_src or n_bkg value after running 'extract curve'."""

    mode = obs_id[-2:]
    ID = obs_id[:-2]
    num = 3 if mode == "pc" else 2

    if type=="src": name= f"./uplims_analysis/src_{mode}_{npix}pix"
    elif type=="bkg": name=f"./uplims_analysis/bkg_{mode}"

    xselect_cmds = f"""
    xselect << EOF
xsel
no
read events ./uplims_analysis/sw{ID}x{mode}w{num}po_cl.evt
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



def extract_spectrum(type, obs_id, npix=None):
    """Create a spectrum file."""

    mode = obs_id[-2:]
    ID = obs_id[:-2]
    num = 3 if mode == "pc" else 2

    if type=="src": name= f"./uplims_analysis/src_{mode}_{npix}pix"
    elif type=="bkg": name=f"./uplims_analysis/bkg_{mode}"

    
    xselect_spectrum_cmds = f"""
    xselect << EOF
xsel
no
read events ./uplims_analysis/sw{ID}x{mode}w{num}po_cl.evt
./
yes
filter region {name}.reg
extract spectrum
save spectrum ./uplims_analysis/{type}_{obs_id}.pha
EOF
    """

    output = run_command(xselect_spectrum_cmds)


# In WT mode, the backscal factor needs editing, as per https://www.swift.ac.uk/analysis/xrt/backscal.php
# Note that the unit doesn't matter, since all we will do with it is calculate the ratio of the source and background backscal factors -- both of which are changed using this function
def change_backscal(type, obs_id, npix=20, r1=80, r2=120):

    filename = f'./uplims_analysis/{type}_{obs_id}.pha'

    if type=="src": N = 2*npix
    elif type=="bkg": N = r2 - r1 -1

    with fits.open(filename, mode='update') as hdul:
        # Inspect the primary HDU or the appropriate extension
        # For .pha files, BACKSCAL is usually in the first or second extension
        for hdu in hdul:
            if "BACKSCAL" in hdu.header:
                hdu.header["BACKSCAL"] = N
                print(f"Updated BACKSCAL to {N}.")
        
        # Save changes back to the file
        hdul.flush()

    

def extract_backscal_exposure(type, obs_id):
    """Extracts BACKSCAL and EXPOSURE from a PHA file header."""

    if type=="src": filename= f"./uplims_analysis/src_{obs_id}.pha"
    elif type=="bkg": filename=f"./uplims_analysis/bkg_{obs_id}.pha"

    with fits.open(filename) as hdul:
        header = hdul[1].header
        backscal = header.get("BACKSCAL")
        exposure = header.get("EXPOSURE")
    return backscal, exposure



## TODO: Checking output of xrtmkarf to get location rmf file to use
def run_xspec(obs_id, nH=0.2, gamma=1.7):
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
    # Fix the norm to 1.0 -- which is actually 1e-12 ergs/s/cm
    Model("tbabs*pegpwrlw", setPars={1:f"{nH} -1",2:f"{gamma} -1",3:1.0, 4:10.0, 5:1.0})

    AllData(1).show()

    # Show rates for the model applied to the spectrum
    count_rate_factor = AllData(1).rate[3] # this is the model predicted rate

    Xset.closeLog()
    
    return count_rate_factor


def get_values_for_uplims(obs_id, ra, dec, count_rate, nH, gamma): 

    mode = obs_id[-2:]
    if mode != "pc" and mode != "wt": raise ValueError(f"Invalid mode: {mode}. Expected 'pc' or 'wt'.")
    num = 3 if mode == "pc" else 2
    ID = obs_id[:-2]
    
    # Delete pre-existing results
    file_paths = [f"./uplims_analysis/bkg_{obs_id}.pha", f"./uplims_analysis/src_{obs_id}.pha", f"./uplims_analysis/src_{obs_id}.arf"]
    for file_path in file_paths:
        if os.path.exists(file_path): os.remove(file_path)

    # Create the region files for the source and background
    if mode=="pc": npix = generate_region_files_pc(ra, dec, count_rate)
    elif mode=="wt": npix = generate_region_files_wt(ra, dec)


    ## Source:
    # Get the number of good counts
    n_src = extract_curve_info("src", obs_id, npix) 
    # Create the .pha file
    extract_spectrum("src", obs_id, npix) 
    # If wt mode, change the backscal of source file
    if mode=="wt": change_backscal("src", obs_id, npix)

    phafile = f"./uplims_analysis/src_{obs_id}.pha"
    expofile = f"./uplims_analysis/sw{ID}x{mode}w{num}po_ex.img"
    if not os.path.exists(phafile):
        raise FileNotFoundError(f"pha file does not exist: {phafile}")
    if not os.path.exists(expofile):
        raise FileNotFoundError(f"expo file does not exist: {expofile}")

    # Run xrtmkarf to generate ARF file
    xrtmkarf_cmd = f"xrtmkarf phafile={phafile} outfile=./uplims_analysis/src_{obs_id}.arf expofile={expofile} srcx=-1 srcy=-1 psfflag=yes"
    print(f"./uplims_analysis/sw{ID}x{mode}w{num}po_ex.img")
    print()
    output = run_command(xrtmkarf_cmd)
    # Get backscal_src and exposure
    backscal_src, exposure = extract_backscal_exposure("src", obs_id)
    
    ## Background:
    # Get the number of good counts
    n_bkg = extract_curve_info("bkg", obs_id)
    # Create the .pha file
    extract_spectrum("bkg", obs_id)
    # If wt mode, change the backscal of source file
    if mode=="wt": change_backscal("bkg", obs_id)
    # Get backscal_bkg
    backscal_bkg, _ = extract_backscal_exposure("bkg", obs_id)
    
    ## Run XSPEC to get count_rate_factor:
    count_rate_factor = run_xspec(obs_id, nH, gamma)
    
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


# xrtmkarf phafile=./uplims_analysis/src_00010627132wt.pha outfile=./uplims_analysis/src_00010627132wt.arf expofile=./uplims_analysis/sw00010627132xwtw2po_ex.img srcx=-1 srcy=-1 psfflag=yes

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
# n_bkg: background counts in an annulus with pre-determined radii
# backscal_src and backscal_bkg are backscaling factors, used to account for differences in source and background area size
# exposure is the exposure time in seconds
# count_rate_factors is a part of the scaling factor between count rate and flux
## To determine the count_rate_factors:
# Fix the flux at 1e-12, gamma (default 1.7), and the known NH. 
# Then (don't fit), use the Xspec command "show rates" to get the denominator. 
# Therefore the conversion is 1e-12/[answer]
def get_uplim_fluxes(n_src= [], n_bkg=[], backscal_src=[], backscal_bkg=[], exposure=[], count_rate_factor = []):

    # Use the appropriate source region size from Evans et al. 2009
    # xselect on the src (or background) region, and it will give the correct number of 'good' counts
    n_src = np.array(n_src) # circle
    n_bkg = np.array(n_bkg) # annulus

    # Read these parameters (i.e., backscal and exposure from header files)
    backscal_src= np.array(backscal_src) 
    backscal_bkg = np.array(backscal_bkg) 
    backscal_factor = backscal_bkg/backscal_src

    # Error calculations
    src_err = (np.array(gehrels_poisson(n_src))-n_src)[1] # element 1 is the upper error
    bkg_err = (np.array(gehrels_poisson(n_bkg))-n_bkg)[1] # element 1 is the upper error

    # Calculate the background counts and error on the background counts, adjusted by backscal_factor
    bkg_in_aperture     = n_bkg/backscal_factor 
    bkg_in_aperture_err = bkg_err/backscal_factor 

    # Calculate net source counts (source counts minus background counts in aperture).
    net = n_src - bkg_in_aperture
    net_err = np.sqrt(src_err*src_err + bkg_in_aperture_err*bkg_in_aperture_err) # net error -- just add errors in quadrature, as a simple scheme

    # Scale the net error for a 3-sigma upper limit on the count rate
    exposure = np.array(exposure)
    count_rate = 3 * net_err / exposure

    # The norm when obtaining count_rate_factor was fixed to flux = 1e-12 ergs/s/cm
    flux_conversion = 1e-12/np.array(count_rate_factor)

    # Convert the count rate to flux units
    flux_uplims = flux_conversion * count_rate

    print('3-sigma upper limit on the flux is: \n' + np.array2string(flux_uplims, separator=",") )

    return flux_uplims


##########################################################################################################################################################################################################################


# Set get_data=False when you are re-running this and the data has already been loaded
def uplims_runner(uplims_files, count_rate_ar, target_coords, nH, gamma):

    ra, dec = str(target_coords[0]), str(target_coords[1])

    ## Get the raw event files for the observations of interest, using the xrtpipeline_runner.sh script
    # This gets data from https://www.swift.ac.uk/swift_portal/ using the wget option and puts it in /data.
    # Then, in /uplims_analysis, its runs xrtpipeline.
    # The required data (event files and exposure files) is then copied to /uplims_analysis.
    if not os.path.exists("./uplims_analysis"):
    
        print("Getting data...")

        os.makedirs("./uplims_analysis")
        os.makedirs("./uplims_analysis/data")
        os.makedirs("./uplims_analysis/xrtpipeline_output")
        
        sh_file = "xrtpipeline_runner.sh"
        try: # Run the shell script with RA, DEC, and obs_id_ar as arguments
            print("Running bash script...")
            result = subprocess.run(['bash', sh_file, ra, dec] + uplims_files, check=True, text=True, capture_output=True)
            #print("Output:\n", result.stdout)
        except subprocess.CalledProcessError as e:
            print("Error:\n", e.stderr)
  
    
    ## Get the values needed to calculate the 3-sigma flux upper limits, using the event and exposure files
    print("Getting values...")
    all_n_src, all_backscal_src, all_exposure, all_n_bkg, all_backscal_bkg, all_count_rate_factor = [], [], [], [], [], []
    for obs_id, count_rate in zip(uplims_files, count_rate_ar): 
      n_src, backscal_src, exposure, n_bkg, backscal_bkg, count_rate_factor = get_values_for_uplims(obs_id, ra, dec, count_rate, nH, gamma)
      all_n_src.append(n_src), all_backscal_src.append(backscal_src), all_exposure.append(exposure), all_n_bkg.append(n_bkg), all_backscal_bkg.append(backscal_bkg), all_count_rate_factor.append(count_rate_factor)

    ## Calculate the 3-sigma flux upper limits.
    print("Calculating fluxes...")
    flux_uplims = get_uplim_fluxes(all_n_src, all_n_bkg, all_backscal_src, all_backscal_bkg, all_exposure, all_count_rate_factor)

    return flux_uplims

    