import subprocess
from xspec import * 
from astropy.io import fits
import re
import os
import numpy as np

"""
Before running this script, the following must be done in the uplims_analysis directory:

(1) 
- mkdir data 
- Get the raw event files for the observations of interest: https://www.swift.ac.uk/swift_portal/
- To do this, use the wget option in the /data sub-folder.

(2)
- mkdir xrtpipeline_output
- In /uplims_analysis, run xrtpipeline for the observations of interest
e.g. xrtpipeline indir=./data/reproc/00016584014 outdir=./xrtpipeline_output/00016584014 steminputs=sw00016584014 srcra=OBJECT srcdec=OBJECT clobber=yes 
- However, in the above command, use the actual source coordinates in place of 'OBJECT'.

(3)
- Copy the required data to /uplims_analysis

For the particular observations of interest:
 - event files (original names)
 - exposure files (original names)

 
"""


## TODO:
# Checking output of xrtmkarf to get location rmf file to use


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
    with open("bkg.reg", "w") as bkg_file:
        bkg_file.write("icrs\n")
        bkg_file.write(f"annulus({ra},{dec},141.439\",259.304\")\n")

    # Write to src_{npix}.reg
    src_filename = f"src_{npix}pix.reg"
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

    if type=="src": name= f"src_{npix}pix"
    elif type=="bkg": name="bkg"

    xselect_cmds = f"""
    xselect << EOF
xsel
no
read events sw{obs_id}xpcw3po_cl.evt
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

    if type=="src": filename= f"src_{obs_id}.pha"
    elif type=="bkg": filename=f"bkg_{obs_id}.pha"

    with fits.open(filename) as hdul:
        header = hdul[1].header
        backscal = header.get("BACKSCAL")
        exposure = header.get("EXPOSURE")
    return backscal, exposure


def extract_spectrum(type, obs_id, npix=None):

    if type=="src": name= f"src_{npix}pix"
    elif type=="bkg": name="bkg"

    xselect_spectrum_cmds = f"""
    xselect << EOF
xsel
no
read events sw{obs_id}xpcw3po_cl.evt
./
yes
filter region {name}.reg
extract spectrum
save spectrum {type}_{obs_id}.pha
EOF
    """

    output = run_command(xselect_spectrum_cmds)
    



def run_xspec(obs_id, nH):
    """Run XSPEC with the specified model and get count_rate_factor."""

    AllData.clear()
    AllModels.clear()

    # Set the abundance table used in the plasma emission and photoelectric absorption models; 'wilm' is an XSPEC built-in abundance table
    Xset.abund = "wilm"    
    Xset.xsect = "vern"  

    Xset.openLog("xspec.log")

    spectrum_file = f"src_{obs_id}.pha"
    AllData(spectrum_file) # Load the spectrum into the AllData object
    AllData(1).response = "/usr/local/shared/caldb/data/swift/xrt/cpf/rmf/swxpc0to12s6_20210101v016.rmf"
    AllData(1).response.arf = f"src_{obs_id}.arf"
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


def main(obs_id, ra, dec, count_rate, nH): 

    # Create the region files for the source and background
    npix = generate_region_files(ra, dec, count_rate)

    ## Source:
    # Get the number of good counts
    n_src = extract_curve_info("src", obs_id, npix) 
    # Create the .pha file
    extract_spectrum("src", obs_id, npix) 
    # Run xrtmkarf to generate ARF file
    xrtmkarf_cmd = f"xrtmkarf phafile=src_{obs_id}.pha outfile=src_{obs_id}.arf expofile=sw{obs_id}xpcw3po_ex.img srcx=-1 srcy=-1 psfflag=yes"
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
    print("RESULTS:")
    print(f"n_src: {n_src}")
    print(f"backscal_src: {backscal_src}")
    print(f"exposure: {exposure}")
    print(f"n_bkg: {n_bkg}")
    print(f"backscal_bkg: {backscal_bkg}")
    print(f"count_rate_factor: {count_rate_factor}")


if __name__ == "__main__":

    #obs_id = input("Enter observation ID: ")
    #ra = input("Source RA: ")
    #dec = input("Source dec: ")
    # count_rate = input("Enter the number of counts/s for the observation: ")
    #nH = input("nH value [10^22 cm^-2]: ")
    obs_id = "00016584014"
    ra = "261.930583"
    dec="-16.205322"
    count_rate = 0.01155
    nH="0.2"

    file_paths = [f"bkg_{obs_id}.pha", f"src_{obs_id}.pha", f"src_{obs_id}.arf"]
    for file_path in file_paths:
        if os.path.exists(file_path): os.remove(file_path)
    
    main(obs_id, ra, dec, count_rate, nH)