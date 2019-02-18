#!/bin/bash

tools=../../build/KITE-tools

for i in test*; do
  cd $i
  ./script.sh
  cd ..

done
