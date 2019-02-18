#!/bin/bash

tools=../../build/KITE-tools

echo Linear optical conductivity of graphene.
$tools configREF.h5 --CondOpt -X > log
python ../compare.py optcond.dat optcondREF.dat
