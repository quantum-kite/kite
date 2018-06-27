# Post-processing tools

*KITE* calculates and stores the Chebyshev moments of a given expansion in the hdf file. In general, it is possible to calculate a quantity at different conditions with the same moments of an expansion, without the need to perform time consuming iterations. *KITE* post-processing tool (*KITE-tools*), based on the hdf file, automatically identifies which quantities should be calculated. By default, when defining the settings with the python script, the user already pre-defines parameters for the postprocessing tools. This is the case, for example, in the number of energy points and temperature in the calculation of the DC conductivity, or the frequencies in the optical conductivity calculation.

Before trying to use it, one needs to compile *KITE-tools*
``` bash
cd /tools/
make
```
It's usage is very simple:
``` bash
./tools//KITE-tools archive.h5
```
where archive.h5 is the hdf files that stores the output of the calculation. If *KITE-tools* does not find this output, it will return an error.

In the table below, we specify the name of the files that are created by *KITE-tools* according to the calculated quantity and the format of the data file.

| Quantity | File | Column 1 | Column 2 | Column 3 |
|:------------------------------------------------------------------------:|:------------------:| ------------ |:--------------:|:--------------:|
| Density of States | dos.dat | energy | DOS [Re] | DOS[Im] |
| Optical Conductivity | optical_cond.dat | Frequency | Opt. Cond [Re] | Opt. Cond [Im] |
| DC Conductivity | condDC.dat | Fermi energy | Cond [Re] | Cond[Im] |
| Non-linear Optical Cond. (experimental feature , currently only for hBN) | nonlinear_cond.dat | Frequency | NL Cond[Re] | NL Cond[Im] |

* All linear conductivities are in units of e2/h

* Both Planck's constant and electron charge are set to 1.

For more details on the type of calculations performed during post-processing, check [Resources][1] where we discuss our method.

The single shot DC conductivity does not need any post-processing as it is an energy dependent calculation where the conductivity is calculated *on the fly*. In this particular case, the data is extracted directly from the hdf file with the following python script:
``` python
import h5py #read h5 files
import numpy as np
file_name = 'archive.h5' #h5 file
file_input = h5py.File(file_name, 'r+')

# extract the single shot
single_shot = file_input['Calculation']['singleshot_conductivity_dc']['SingleShot']
np.savetxt('cond_singleshot.dat', np.c_[single_shot[:, 0], single_shot[:, 1]])
```

[1]: https://quantum-kite.com/resources/
