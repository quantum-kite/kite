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

    ../KITE-tools config.h5 --ARPES -K green 0.1 -E -10 10 2048 -F 100
    python ../../tools/process_arpes.py arpes.dat

    #python test.py
    rm -r __pycache__
fi


if [[ "$1" == "quick" ]]; then
    # Run KITEx immediately on the existing configuration file
    cp configORIG.h5 config.h5
    SEED=$S ../KITEx config.h5 > log_KITEx
    ../KITE-tools config.h5 --ARPES -K green 0.1 -E -10 10 2048 -F 100
    python ../../tools/process_arpes.py arpes.dat
    #python test.py

fi

