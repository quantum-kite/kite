#!/bin/bash

S=3
 
# This script will compare the .h5 file
if [[ "$1" == "redo" ]]; then
    # Recreate all the .h5 configuration files from scratch

    python config.py > log_config
    chmod 755 config.h5
    cp config.h5 configORIG.h5 

    SEED=$S ../KITEx config.h5 > log_KITEx
    cp config.h5 configREF.h5

    ../KITE-tools config.h5 > log_KITEx
    cp config.h5 configREF.h5
    cp dos.dat dosREF.dat

    python test.py
    rm -r __pycache__
fi

if [[ "$1" == "quick" ]]; then
    # Run KITEx immediately on the existing configuration file
    cp configORIG.h5 config.h5
    python test.py

fi

