#!/bin/bash

tools=../../build/KITE-tools

echo Comparing the optical conductivity for different block divisions.
$tools configREF.h5 --CondOpt -O 0 10 16 -E 16 -C 4 4 -X -N opt_44.dat > log1
$tools configREF.h5 --CondOpt -O 0 10 16 -E 16 -C 2 2 -X -N opt_22.dat > log1
$tools configREF.h5 --CondOpt -O 0 10 16 -E 16 -C 1 1 -X -N opt_11.dat > log2
$tools configREF.h5 --CondOpt -O 0 10 16 -E 16 -C 4 2 -X -N opt_42.dat > log4
$tools configREF.h5 --CondOpt -O 0 10 16 -E 16 -X -N opt.dat > log4

paste g64_d64_opt_* opt.dat | awk 'BEGIN {sum=0} 
    {sum+=($3-$15)^2 + ($6-$15)^2 + ($9-$15)^2 + ($12-$15)^2} 
    END {print sum}'

