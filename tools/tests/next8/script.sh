#!/bin/bash
cat description.txt

tools=../../build/KITE-tools
$tools configREF.h5 --ARPES -M 1024 -K green 0.01 -E -4 -3 100 > log
python process_arpes.py arpes.dat >> log
mv dos_k-0.0.dat dos_k.dat 
rm arpes.dat
#python ../compare.py dos.dat dosREF.dat
