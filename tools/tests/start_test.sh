#!/bin/bash

tools=../../build/KITE-tools

for i in test*; do
  cd $i
  filename=$(cat function)
  par=$(cat params)
  cat description
  $tools configREF.h5 $par > /dev/null
  norm=$(python ../compare.py ${filename}.dat ${filename}REF.dat)
  echo $norm
  cd ..

done
