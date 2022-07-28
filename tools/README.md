# `kite/tools`-folder
This folder contains two distinct tools:

## 1. `KITE-tools`
The main KITE-tools post-processing program that is used to reconstruct calculated observable
out of the data contained in the processed HDF5 file.
This tool needs to be compiled into a binary executable before being used.

## 2. `process_single_shot.py`
to extract the output data from a single-shot calculation of the longitudinal DC-conductivity.
This tool is a simple python script that outputs a data file that contains three columns of data that
are read out of the processed HDF5, namely:

| Column 1       | Column 2                      | Column 3                           |
|----------------|-------------------------------|------------------------------------|
| Fermi Energies | Real Part of the Conductivity | Imaginary Part of the Conductivity |
