#!/usr/bin/python

import os, optparse, sys
import numpy as np
from astropy.io import fits

np.set_printoptions(threshold=np.inf) 

## Parse the arguments to the python file
def parse_args( argv):

    desc="""%prog creates a grouped pha spectral file suitable for Xspec over a specified energy range. \
    User must specify input and ouput file. User can specifiy minimum binning based on spectral resolution and minimum number of counts per bin."""
    
    parser = optparse.OptionParser(description=desc)

    parser.add_option('-i', '--input', \
                    help='input pi file', \
                    dest='input_pi')
    parser.add_option('-r', '--rmf', \
                        help='rmf file', \
                        dest='rmf')
    parser.add_option('-a', '--arf', \
                        help='arf file', \
                        dest='arf')
    parser.add_option('-b', '--background', \
                        help='background pi file', \
                        dest='background')
    parser.add_option('-o', '--output', \
                        help='output (grouped) pi file', \
                        dest='output_pi')
    parser.add_option('-e', '--energy_bin', \
                        help='minimum energy change per bin in eV or fraction [%default]', \
                        dest='bin_min_energy', \
                        type='float', \
                        default=100.0)
    parser.add_option('-f', '--fractional_energy_binning', \
                        help='switch to set fractional energy binning [%default]', \
                        dest='fraction_switch', \
                        default=False, \
                        action='store_true')
    parser.add_option('-l', '--lower_energy', \
                        help='minimum energy of first bin in eV [%default]', \
                        dest='lower_energy', \
                        type='float', \
                        default=300.0)
    parser.add_option('-u', '--upper_energy', \
                        help='maximum energy of last bin in eV [%default]', \
                        dest='upper_energy', \
                        type='float', \
                        default=10000.0)
    parser.add_option('-c', '--count_bin', \
                        help='minimum counts per bin [%default]', \
                        dest='bin_min_counts', \
                        type='long', \
                        default=25)

    if (argv == []):
        parser.print_help()
        exit(-1)

    return parser



## Check whether the arguments are valid
def valid_args_check(opts, parser):

    # Check the input file exists and is well specified
    if opts.input_pi is None:
        print("The input pi file is missing\n")
        parser.print_help()
        exit(-1)
    if os.path.isfile(opts.input_pi) == 0:
        print("The specified input pi file ("+opts.input_pi+") does not exist\n")
        exit(-1)

    input_fits   = fits.open(opts.input_pi)
    input_pi     = input_fits[1]
    #  input_pi_hdr = input_pi.header.ascardlist()
    input_pi_hdr = input_pi.header


    # Check the output file exists and is well specified
    if opts.output_pi is None:
        print("The output pi file is missing\n")
        parser.print_help()
        exit(-1)

    # Check the RMF file exists and is well specified
    if opts.rmf is None:
        if 'RESPFILE' in input_pi_hdr:
            opts.rmf = input_pi_hdr['RESPFILE']
            print("Found RESPFILE: {}".format(opts.rmf))
        if opts.rmf is None:
            print("No rmf file is specified in either input pi file or command line\n")
            parser.print_help()
            exit(-1)

    # Check the background file exists and is well specified
    if opts.background is None:
        if 'BACKFILE' in input_pi_hdr:
            opts.background = input_pi_hdr['BACKFILE']
            print("Found BACKFILE: {}".format(opts.background))
        if opts.background is None:
            print("*WARNING*: The grouped spectra has no background file\n")

    # Check the ARF file exists and is well specified
    if opts.arf is None:
        if 'ANCRFILE' in input_pi_hdr:
            opts.arf = input_pi_hdr['ANCRFILE']
            print("Found ANCRFILE: {}".format(opts.arf))
        if opts.arf is None:
            print("*WARNING*: The grouped spectra has no arf file\n")

    input_fits.close



## Grouping counts to ensure good statistics
def group_pha(opts, pi_data, ebounds):

    # print("Pi data shape: ", pi_data.shape) # (1024,)
    pi_rows = pi_data.shape[0] 
    print("Pi rows: ", pi_rows) # there are 1024 detector channels
    
    quality             = np.zeros(pi_rows)+5 # quality of data; set to 5 at the start, which means bad
    grouping            = np.zeros(pi_rows)+1 # to track the grouping; set all to 1 at the start

    channels            = ebounds.data.field('CHANNEL')
    e_min               = ebounds.data.field('E_MIN')*1000
    # print("Emin shape: ", e_min.shape) # (1024,)
    e_max               = ebounds.data.field('E_MAX')*1000
    # print("Emax shape: ", e_max.shape) # (1024,)
    energy              = np.sqrt(e_min*e_max)
    energy_width        = e_max-e_min # width in energy of each each channel
    min_width_at_energy = np.zeros(pi_rows) + opts.bin_min_energy # set minimum energy change per channel
    
    if (opts.fraction_switch): # using fractional energy binning
        min_width_at_energy = opts.bin_min_energy * energy # scale bin width by energy
        
    # Initiate values
    tot_counts          = 0
    row_counter         = 0
    last_new_row        = 0
    counts_in_bin       = opts.bin_min_counts
    energy_width_of_bin = np.sum(energy_width) # i.e. total energy of the whole spectrum

    # AKH added this to find the index that corresponded to counts
    count_index = np.where(np.array(pi_data.columns.names).astype(str) == 'COUNTS')[0][0] # = 1

    ## TODO: For row_counter, can just use enumerate

    all_counts=[]

    # For each channel in pi_data:
    for row in (pi_data): # row is type FITS_record
        
        pi_counts = row[count_index] # This was originally an index of 2 potentially a difference between chandra v. swift
        all_counts.append(pi_counts)

        pi_row = row[0] # This is also just the channels so I'm not sure what the difference is
        matching_index = (np.where(channels==pi_row))[0][0]

        # If the energy of the current channel falls in the specified energy range
        if ((e_min[matching_index] > opts.lower_energy) and (opts.upper_energy > e_max[matching_index])):
            quality[row_counter] = 0 # quality is marked as good (0)
            tot_counts += pi_counts # total counts for that channel are added to tot_counts

            # Check whether the counts in the bin and bin width meet the minimum requirements.
            # If so, a new group/bin is started.
            if((counts_in_bin >= opts.bin_min_counts) and (energy_width_of_bin >= min_width_at_energy[row_counter])):
                counts_in_bin         = pi_counts
                energy_width_of_bin   = energy_width[row_counter]
                last_new_row          = row_counter
                grouping[row_counter] =1
            # Otherwise, the counts and energy for the current channel are added to the existing group/bin.
            else:
                counts_in_bin         += pi_counts
                energy_width_of_bin   += energy_width[row_counter]
                grouping[row_counter] =-1 # i.e. this channel is part of an ongoing group/bin

        row_counter += 1
    

    if counts_in_bin < opts.bin_min_counts:
        grouping[last_new_row] =-1

    print()
    print("ALL COUNTS:")
    ar = np.array(all_counts)
    print(ar)
    print(len(all_counts))
    print(np.sum(all_counts))
    indices = np.where(ar == 1)
    print(indices[0])
    
    return (quality, grouping, tot_counts)



## Wrapper to group counts and output an updated FITS file
def fitsio_grppha(opts):

    # Open PHA FITS file containing the photon counts
    input_fits   = fits.open(opts.input_pi)

    # Open RMF which contains energy information
    ebounds = fits.open(opts.rmf)['EBOUNDS']

    # Use the group_pha function to group counts
    (quality, grouping, tot_counts) = group_pha(opts, input_fits[1].data, ebounds)

    print("Total counts: ", tot_counts)
    
    # Adjust the binning if counts are too low
    if tot_counts < 300: 
        print("Bin counts are too low. Adjusting minimum counts per bin.")
        opts.bin_min_counts = 1 # adjust the bin min counts to 1 (default is 25)
        (quality, grouping, tot_counts) = group_pha(opts, input_fits[1].data, ebounds)


    print(f'Final min counts per bin: {opts.bin_min_counts}')
    
    print("QUALITY: ", quality)
    print("GROUPING: ", grouping)
    
    # Create new columns for quality and grouping
    quality = fits.Column(name='QUALITY', format='I', array=quality)
    grouping = fits.Column(name='GROUPING', format='I', array=grouping)

    # Delete existing quality and grouping columns if present
    if 'QUALITY' in input_fits[1].columns.names:
        print('Deleting QUALITY column')
        input_fits[1].columns.del_col('QUALITY')
    if 'GROUPING' in input_fits[1].columns.names:
        print('Deleting GROUPING column')
        input_fits[1].columns.del_col('GROUPING')


    # Make a copy of the input FITS file for the output
    output_fits = input_fits.copy() # This is a reference not a unique object despite it claiming to be a unique object
    # Add the new quality and grouping columns to the first HDU
    output_fits[1] = fits.BinTableHDU.from_columns(input_fits[1].columns + quality + grouping, header=input_fits[1].header)
    # Update the FITS headers with additional information
    output_fits[1].header['BACKFILE'] = opts.background
    output_fits[1].header['ANCRFILE'] = opts.arf
    output_fits[1].header['COUNTGRP'] = opts.bin_min_counts
    # Write the modified FITS file
    output_fits.writeto(opts.output_pi, overwrite=True)
    output_fits.close()
    input_fits.close()



def main(argv):
    # Parse the arguments
    parser = parse_args(argv) # initialise the parser.
    (opts, args) = parser.parse_args() # opts: an object that contains the parsed options as attributes; args: a list of positional arguments not associated with options.
    # Check that the inputted arguments are valid
    valid_args_check(opts, parser) 
    # Run wrapper method to group counts and output file
    fitsio_grppha(opts) 

if __name__ == "__main__":
    main(sys.argv[1:])
