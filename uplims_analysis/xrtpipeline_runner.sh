#!/bin/bash
# To run this: addqueue -c "4 hours" -n 1 runner.sh

## Load the necessary module
module load heasoft/6.33
source ~/.bashrc
export HEADASNOQUERY=1 # allows xrtpipeline to proceed without needing access to /dev/tty (i.e. prevents it from attempting to prompt the user during execution)


## RUN XRTPIPELINE

names=(00089766002 00089766003 00089766004 00089766005 00089766006 00089766007 00089766012 00016584002 00016584004 00016584005 00016584006 00016584007 00016584008 00016584009 00016584010 00016584011 00016584012 00016584014 00016584016 00016584017 00016584018 00016584019)

# Loop through each name in the array
for name in "${names[@]}"; do
  echo "Processing $name..."

  xrtpipeline indir=./data/reproc/$name outdir=./xrtpipeline_output/$name steminputs=sw${name} srcra=261.930583 srcdec=-16.205322 clobber=yes
  
  # Optional: Check for errors or completion status
  if [ $? -ne 0 ]; then
    echo "Error processing $name"
  else
    echo "Completed $name successfully"
  fi
done