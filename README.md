Lattice building script is used for making and exporting the lattice.
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

### ~~Convert pybinding format to common format used in C++~~

### ~~Add  onsite disorder:~~
 - For now only constant values as onsite potential and the labels for the
 rectangular/uniform and gaussian distribution are supported.
 Both distributions have mean value and width.

### ~~Magnetic field:~~
 - Add magnetic field support. For now only on and off due to PBC.

### ~~Chose between the functions that you want to calculate:~~
 - Other strings that the labels result in an error.   
  - labels are DOS, CondXX, CondXY, OptCond, SpinCond
  - DOS (number of disorder, num of moments, num of random vectors) 'DOS' == 1,
  - Conductivity (number of disorder, num of moments, num of random vectors)
  'CondXX' == 2, 'CondXY' ==   3,
  - Optical Conductivity (number of disorder, num of moments, num of random
    vectors), 'OptCond' == 4
  - Spin Dependent Conductivity (number of disorder, num of moments, num of
    random vectors), 'SpinCond' == 5.

# A list of things that could/should be implemented

### Define plotting from the resulted data
### Add vacancy disorder modifications
### Add the choice of magnetic field to the user.
  - In case of periodic boundary conditions give available values.
