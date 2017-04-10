Lattice building script is used for making and exporting the lattice.

# A list of things that should be implemented

### ~~Add labels for onsite potential:~~
 ~~- For now only constant values exported from python script are supported.
 Add the labels for the rectangular and gaussian distribution.~~
 Labels are added. Other strings than the labels results in errors.

### ~~Magnetic field:~~
 - Add magnetic field support. For now only on and off due to PBC.

### ~~Chose between the functions that you want to calculate:~~
 - Other strings that the labels should result in errors.   
 - labels are DOS, CondXX, CondXY, SpinCond
 - DOS (number of disorder, num of moments, num of random vectors),
 - Conductivity (number of disorder, num of moments, num of random vectors),
 - Optical Conductivity (number of disorder, num of moments, num of random vectors),
 - Spin Dependent Conductivity (number of disorder, num of moments, num of random vectors).

### Define ploting from the resulted data