#!/bin/bash

tools=../../build/KITE-tools

echo DC conductivity of the square lattice.
$tools configREF.h5 --CondDC -X
python3 ../compare.py condDC.dat condDCREF.dat
