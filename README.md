Lattice building script is used for making and exporting the lattice. Recent changes, and the features implemented in
the scripts that are used for making and exporting different system from Pybinding to the C++ code are shown here.
# Building Lattice
  the building is defined in 3 scripts,

  * lattice_building.py - user changes this, defines the data he needs, lattice,
  disorder and field modifications, calculation type, size, periodicity...

  * define_lattice.py - script with available lattices, user can use one from
  here, or define similar by it's own.

  * export_lattice.py - user doesn't need to think about this, changes the format
  of pybinding lattice the one in C++ code, and exports everything needed.

  in addition:

  * Building_script - jupyter notebook script -  imports the define_lattice.py and
  export_lattice.py and does the user specified building same as in
  lattice_building.py

# Implemented things


# Changing the structure of the lattice
The pybinding format is successfully convert to the format of the C++ code.

### ~~Add disorder:~~
 - Disorder is added to the plain lattice, defined with respect to the structure.
 (Previously this point was called onsite potential, it's the same.)

### ~~Add structural disorder:~~
 - Structural disorder is added to the plain lattice, affects only a certain area of the structure and can be defined as
 bond disorder or onsite disorder.

### ~~Magnetic field On/Off:~~
 - Add magnetic field support. For now only on and off due to PBC.

### ~~Chose between the functions that you want to calculate:~~
 - Other strings that the labels should result in errors.   
 - labels are DOS, CondXX, CondXY, SpinCond
 - DOS (number of disorder, num of moments, num of random vectors),
 - Conductivity (number of disorder, num of moments, num of random vectors),
 - Optical Conductivity (number of disorder, num of moments, num of random vectors),
 - Spin Dependent Conductivity (number of disorder, num of moments, num of random vectors).

### ~~Pauli matrices.~~
  - These matrices can be found in pybinding.constants.pauli.

## Define ploting from the resulted data

# A list of things that could/should be implemented

### Define plotting from the resulted data
### Add the choice of magnetic field to the user.
  - In case of periodic boundary conditions give available values.