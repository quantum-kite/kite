#!/bin/bash
 
# This script will compare the .h5 file


S=3  # the seed that will be used by KITE
file=/Calculation/dos/MU



if [[ "$1" == "redo" ]]; then
    # Use the existing configuration file
    cp configORIG.h5 config.h5
    chmod 755 config.h5
    SEED=$S ../KITEx config.h5 > log_KITEx
    python3 ../compare.py configREF.h5 $file config.h5 $file
fi

if [[ "$1" == "script" ]]; then
    # Create a configuration file from scratch
    python3 config.py > log_config
    SEED=$S ../KITEx config.h5 > log_KITEx
    python3 ../compare.py configREF.h5 $file config.h5 $file
    rm -r __pycache__
fi

