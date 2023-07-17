#!/bin/bash
 
S=3

# This script will compare the .h5 file
if [[ "$1" == "redo" ]]; then
    # Recreate all the .h5 configuration files from scratch
    python config.py > log_config

    chmod 755 config.h5
    SEED=$S ../KITEx config.h5 > log_KITEx

    python test.py
    rm -r __pycache__
fi


if [[ "$1" == "quick" ]]; then
    python test.py

fi

