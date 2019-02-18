#!/bin/bash

tools=../../build/KITE-tools

echo Testing CondDC for the same Hamiltonian with and without energy shift.
$tools config_sym.h5 --CondDC -N condDC_sym.dat -X -S 0.1 > log1
$tools config.h5 --CondDC -N condDC.dat -X -S 0.1 > log2


paste condDC_sym.dat condDC.dat | awk '
BEGIN {max=0; ORS=" "} 
{sum += ($2-$5)^2} 
sqrt(($2-$5)^2)>max {max = sqrt(($2-$5)^2)} 
END { print max; print sqrt(sum); ORS="\n"; print ""}'

