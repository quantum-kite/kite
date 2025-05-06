!!! Info

    This file describes the contents of the `#!bash kite/tools`-folder. Note that the main post-processing tool is located in the `#!bash kite/build`-folder.

# `kite/tools`-folder
This folder contains bespoke tools:

## 1. `process_arpes.py`
to process the output of a spectral function calculation.

## 2. `process_single_shot.py`
to extract the output data from a single-shot calculation of the longitudinal DC-conductivity.
This tool is a simple python script that outputs a data file that contains three columns of data that
are read out of the processed HDF5, namely:

| Column 1       | Column 2                      | Column 3                           |
|----------------|-------------------------------|------------------------------------|
| Fermi Energies | Real Part of the Conductivity | Imaginary Part of the Conductivity |
