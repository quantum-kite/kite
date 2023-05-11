#!/bin/bash
file="../large000_KITEx_sq_dos/config.h5"

if [[ -f "$file" ]]; then
    sleep 0
else
    echo "Cannot run KITE-tools without the configuration script"
fi

# This script will compare the .h5 file
if [[ "$1" == "redo" ]]; then
    # Recreate all the .h5 configuration files from scratch

    ../KITE-tools $file > log_KITEx
    cp dos.dat dosREF.dat

    python test.py
fi

if [[ "$1" == "quick" ]]; then
    # Run KITE-tools immediately on the existing configuration file

    python test.py

fi

