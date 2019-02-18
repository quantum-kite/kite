#!/bin/bash

tools=../../build/KITE-tools

echo Testing DoS for the same Hamiltonian with and without energy shift.
$tools config_sym.h5 --DOS -N dos_sym.dat -X > log1
$tools config.h5 --DOS -N dos.dat -X > log2

head -n 900 dos.dat > dos1.dat
tail -n 900 dos_sym.dat > dos2.dat

paste dos1.dat dos2.dat | awk '
BEGIN {max=0; ORS=" "} 
{sum += ($2-$5)^2} 
sqrt(($2-$5)^2)>max {max = sqrt(($2-$5)^2)} 
END { print max; print sqrt(sum); ORS="\n"; print ""}'

