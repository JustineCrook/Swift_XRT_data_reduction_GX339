#!/bin/bash
# To run this: addqueue -c "4 hours" -n 1 runner.sh

## Load the necessary module
module load heasoft/6.33
source ~/.bashrc
export HEADASNOQUERY=1 # allows xrtpipeline to proceed without needing access to /dev/tty (i.e. prevents it from attempting to prompt the user during execution)


## RUN XRTPIPELINE

# Assign RA and DEC from the first two arguments
ra=$1
dec=$2

# Shift the first two arguments, so only the names remain
shift 2
names=("$@")

cd ./uplims_analysis

# Loop through each name in the array
for name in "${names[@]}"; do
  echo "Processing $name with RA=$ra and DEC=$dec..."

  cd ./data
  wget -e 'robots=off' -nv -w 2 -nH --cut-dirs=1 -r --no-parent --reject "index.html*" https://www.swift.ac.uk/archive/reproc/$name/xrt/
  wget -e 'robots=off' -nv -w 2 -nH --cut-dirs=1 -r --no-parent --reject "index.html*" https://www.swift.ac.uk/archive/reproc/$name/auxil/ # auxil files
  cd ..
  
  # Run the xrtpipeline command
  xrtpipeline indir=./data/reproc/$name outdir=./xrtpipeline_output/$name steminputs=sw${name} srcra=$ra srcdec=$dec clobber=yes
  
  # Optional: Check for errors or completion status
  if [ $? -ne 0 ]; then
    echo "Error processing $name"
  else
    echo "Completed $name successfully"
  fi

  # Copy the required files one level up in the directory structure
  cp ./xrtpipeline_output/$name/sw${name}xpcw3po_cl.evt ./
  cp ./xrtpipeline_output/$name/sw${name}xpcw3po_ex.img ./
  
  # Check if the copy was successful
  if [ $? -ne 0 ]; then
    echo "Error copying files for $name"
  else
    echo "Files copied successfully for $name"
  fi
  
done