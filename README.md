Lattice building script is used for making and exporting the lattice.

# A list of things that should be implemented

### ~~Add  onsite disorder:~~
 ~~- For now only constant values exported from python script are supported.
 Add the labels for the rectangular and gaussian distribution.~~

### ~~Magnetic field:~~
 - Add magnetic field support. For now only on and off due to PBC.

### ~~Chose between the functions that you want to calculate:~~
 - Other strings that the labels should result in errors.   
 - labels are DOS, CondXX, CondXY, SpinCond
 - DOS (number of disorder, num of moments, num of random vectors) 'DOS' == 1,
 - Conductivity (number of disorder, num of moments, num of random vectors) 'CondXX' == 2, 'CondXY' == 3,
 - Optical Conductivity (number of disorder, num of moments, num of random vectors), 'OptCond' == 4
 - Spin Dependent Conductivity (number of disorder, num of moments, num of random vectors),  'SpinCond' == 5.

### Define ploting from the resulted data