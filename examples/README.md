!!! Info

    This file describes the contents of the [`#!bash kite/examples/`][examples-folder-github]-folder.

## KITE Scripts
This folder contains various scripts with examples showing the functionalities of the KITE library. To execute a script, run

``` bash
python3 example.py
../build/KITEx example-output.h5
../build/KITE-tools example-output.h5
```

Several examples are given, including

* A simple DOS calculation of various systems, such as
    * Checkerboard lattice
    * Graphene
    * Cubic lattice
  
* Two types of disorder:
  * On-site (uncorrelated) disorder
  * Vacancies
  
* Optical conductivity
* XX/YY conductivity
* Weyl semi-metal calculations
* Fu-Kane-Mele model calculations

The examples start with the initialization of the Pybinding Lattice in the Python-script.
The type of calculation is then described in the *main()*-part, with some default parameters associate with this calculation.
The python-script returns a HDF5-file.
This file can then be passed to the KITEx-executable to do the core calculation.
Finally, the raw calculations are analyzed using the KITE-tools-executable.
Depending on the type of calculation, various output files will be created in the folder where the code was executed from.

If the system is installed in the right folder as mentioned in the tutorial, all the results can be generated automatically by running

``` bash
python3 run_all_examples.py
```

After running this command, all the examples will be executed. This can take several minutes.
Besides the outputfiles, like *name-dos.dat*, plots will be given for the *DOS*, *optical conductivity* and *DC conductivity*.

To clean up the folder after running all the examples, execute the following commands

``` bash
python3
>>> import run_all_examples
>>> run_all_examples.clean()
>>> exit()
```

[examples-folder-github]: https://github.com/quantum-kite/kite/tree/master/examples