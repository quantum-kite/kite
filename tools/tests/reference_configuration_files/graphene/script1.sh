#!/bin/bash



echo "Generating reference hdf files"


PYTHONPATH=../../../../
export PYTHONPATH
name=REF
program=../../../../pp
python generate_hdf.py 2 $name 512 512


$program ${name}dos.h5
