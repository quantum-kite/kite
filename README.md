Lattice building script is used for making and exporting the lattice.
# Building Lattice
  the builing is defined in 3 scripts,

  lattice_building.py - user changes this, defines the data he needs, lattice,
  disorder and field modifications, calculation type, size, periodicity..

  define_lattice.py - script with available lattices, user can use one from
  here, or define similar by it's own.

  export_lattice.py - user doesn't need to think about this, changes the format
  of pybinding lattice to our, and exports everything needed.

  Building_script - jupyter notebook script -  imports the define_lattice.py and
  export_lattice.py and does the user specified building same as in
  lattice_building.py

# Implemented things

### ~~Convert pybinding format to common format used in C++~~

### ~~Add  onsite disorder:~~
 ~~- For now only constant values exported from python script are supported.
 Add the labels for the rectangular and gaussian distribution.~~

### ~~Magnetic field:~~
 ~~- Add magnetic field support. For now only on and off due to PBC.~~

### ~~Chose between the functions that you want to calculate:~~
 ~~- Other strings that the labels should result in errors.   
 - labels are DOS, CondXX, CondXY, SpinCond
 - DOS (number of disorder, num of moments, num of random vectors) 'DOS' == 1,
 - Conductivity (number of disorder, num of moments, num of random vectors) 'CondXX' == 2, 'CondXY' == 3,
 - Optical Conductivity (number of disorder, num of moments, num of random vectors), 'OptCond' == 4
 - Spin Dependent Conductivity (number of disorder, num of moments, num of random vectors),  'SpinCond' == 5.~~

# A list of things that should be implemented

### Define ploting from the resulted data
