#!/bin/bash

python config.py 
../../../build/KITEx config.h5 
../../../tools/build/KITE-tools config.h5 --DOS -E 8192


