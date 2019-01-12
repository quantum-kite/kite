#!/bin/bash

# Test the 2nd-order optical conductivity
../../build/KITEx config.h5

echo "dif      max      norm1    norm2" ; for i in 0 1 2 3; do file=/Calculation/conductivity_optical_nonlinear/Gamma${i}yyy; python ../compare.py configREF.h5 $file config.h5 $file; done

