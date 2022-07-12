## KITE Examples
This folder contains various examples showing the functionalities of the KITE library. To execute an example, run

    python3 example.py
    ../build/KITEx example-output.h5
    ../tools/build/KITE-tools example-output.h5

The following examples are given:
* A simple DOS calculation of various systems, like
    * Checkerboard lattice
    * Mixed disordered lattice
    * On-site disorder
    * Simple square lattice
    * Vacancies
* Optical conductivity
* XX/YY conductivity

The examples start with the initialization of the Pybinding Lattice in the Python-script.
The type of calculation is then described, with some default parameters associate with this calculation.
The python-script returns a HDF5-file.
This file can then be passed to the KITEx-executable to do the core calculation.
Finally, the raw calculations are analyzed using the KITE-tools-executable.
Depending on the type of calculation, various output files will be created in the folder where the code was executed from.