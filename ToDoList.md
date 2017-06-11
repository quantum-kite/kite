
# A list of things that should be implemented

# Changing the structure of the lattice
The pybinding format is successfully convert to the format of the C++ code.

### ~~Add disorder:~~
 - Disorder is added to the plain lattice. (Previouslly this point was called onsite potential, it's the same.)

### ~~Magnetic field On/Off:~~
 - Add magnetic field support. For now only on and off due to PBC.
## Magnetic field arbitrary:
 - Implementation of arbitrary magnetic field.

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
