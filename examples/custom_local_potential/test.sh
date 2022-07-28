#!/bin/bash
 
# This script will compare the .h5 file

set -e
SEED=ones


# Replace the custom file
cp aux.cpp ../../Src/Hamiltonian/aux.cpp
cd ../../build
cmake .. && make 
cd ../examples/custom_local_potential

# Run the example with the modified local potential
python config.py > log_config
SEED=$S ../../build/KITEx config.h5 > log_KITEx
#python ../compare.py configREF.h5 $file config.h5 $file
rm -r __pycache__

# Put the file back
cp aux_default.cpp ../../Src/Hamiltonian/aux.cpp
cd ../../build
cmake .. && make 
cd ../examples/custom_local_potential

