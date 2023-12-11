#!/bin/bash
make -f fit_make_file
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sps/nemo/scratch/spratt/analysis/
./fitting_analysis.exe
