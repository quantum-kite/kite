#!/bin/bash
 
# This script will compare the .h5 file


S=3  # the seed that will be used by KITE
file=/Calculation/dos/MU



dir=../tests/next_14_customlocal

if [[ "$1" == "redo" ]]; then
    
    # Replace the library in the lib directory
    cd ../../lib
    rm libaux.so.1
    ln -s $dir/libaux.so libaux.so.1
    cd $dir

    # Use the existing configuration file
    cp configORIG.h5 config.h5
    chmod 755 config.h5
    SEED=$S ../KITEx config.h5 > log_KITEx
    python ../compare.py configREF.h5 $file config.h5 $file

    ## Put the library back the way it was, so it does not change the other tests
    cd ../../lib
    rm libaux.so.1
    ln -s libaux.so libaux.so.1
    cd $dir
fi

if [[ "$1" == "script" ]]; then
    # Create a configuration file from scratch
    python config.py > log_config
    SEED=$S ../KITEx config.h5 > log_KITEx
    python ../compare.py configREF.h5 $file config.h5 $file
    rm -r __pycache__
fi

