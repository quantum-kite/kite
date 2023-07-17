#!/bin/bash
 
S=3
tools=../KITE-tools

# This script will compare the .h5 file
if [[ "$1" == "redo" ]]; then
    # Recreate all the .h5 configuration files from scratch
    python config.py > log_config
    python config_sym.py >> log_config

    chmod 755 config.h5
    chmod 755 config_sym.h5
    cp config.h5 configORIG.h5 
    cp config_sym.h5 config_symORIG.h5 

    SEED=$S ../KITEx config.h5 > log_KITEx
    SEED=$S ../KITEx config_sym.h5 >> log_KITEx
    cp config.h5 configREF.h5
    cp config_sym.h5 config_symREF.h5

    #echo Testing DoS for the same Hamiltonian with and without energy shift.
    $tools config_sym.h5 --DOS -N dos_sym.dat -X > log
    $tools config.h5 --DOS -N dos.dat -X >> log

    #python test.py
    rm -r __pycache__
fi


if [[ "$1" == "quick" ]]; then
    # Run KITEx immediately on the existing configuration file
    cp configORIG.h5 config.h5
    cp config_symORIG.h5 config_sym.h5

    #python test.py
    #echo Testing DoS for the same Hamiltonian with and without energy shift.
    $tools config_sym.h5 --DOS -N dos_sym.dat -X > log
    $tools config.h5 --DOS -N dos.dat -X >> log
    echo Done

fi

