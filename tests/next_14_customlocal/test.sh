#!/bin/bash
 
# This script will compare the .h5 file


S=3  # the seed that will be used by KITE
file=/Calculation/dos/MU




if [[ "$1" == "redo" ]]; then
    
    # Replace the library in the lib directory
    rm ../../lib/libaux.so.1
    ln -s libaux.so ../tests/next_14_customlocal/libaux.so.1
    cp libaux.so.1 ../../lib/

    # Use the existing configuration file
    #cp configORIG.h5 config.h5
    #chmod 755 config.h5
    #SEED=$S ../KITEx config.h5 > log_KITEx
    #python ../compare.py configREF.h5 $file config.h5 $file

    ## Put the library back the way it was, so it does not change the other tests
    #rm ../../lib/libaux.so.1
    #ln -s ../../lib/libaux.so ../../lib/libaux.so.1
fi

if [[ "$1" == "script" ]]; then
    # Create a configuration file from scratch
    python config.py > log_config
    SEED=$S ../KITEx config.h5 > log_KITEx
    python ../compare.py configREF.h5 $file config.h5 $file
    rm -r __pycache__
fi

