#!/bin/bash
make -f simultaneous_fit_make_file
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sps/nemo/scratch/spratt/analysis/
./simultaneous_fit.exe
