#!/bin/bash

tools=../../build/KITE-tools

echo Density of states of the square lattice.
$tools configREF.h5 --DOS -X -N dos.dat > log
python ../compare.py dos.dat dosREF.dat
