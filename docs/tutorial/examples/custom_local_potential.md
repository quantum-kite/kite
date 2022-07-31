!!! Info

    This file describes the contents of the `#!bash kite/examples/custom_local_potential`-folder.

This example shows how to use the custom local potential functionality in KITE. 

Begin by defining the local energy function that you want to use.
This is a function that takes as arguments the position ($x$, $y$) (in 2D) or ($x$, $y$, $z$) (in 3D) and
orbital `#!python "orb"` and returns the value of the local energy at that position and orbital.
This function is defined by the user in a simple cpp file `#!bash aux.cpp`, which is then compiled into a
shared library to be used by KITE `#!bash libaux.so.1`.
To force KITE to use the function defined by the user, the file `#!bash /lib/libaux.so.1` should be replaced by
our custom shared library, and the flag `#!python custom_local=True` must be used in the python configuration script.
By default, this flag is set to `#!python False`. 

## 1. Compiling and replacing the default library
Begin by compiling the c++ file into a shared library by executing the following lines:

``` bash
g++ -Wall -fPIC -c aux.cpp
g++ -shared -Wl,-soname,libaux.so.1 -o libaux.so aux.o
```

The first line compiles the `#!bash aux.cpp` file and the second generates the shared library from the compiled object.
Next, replace the default shared library file in the `#!bash /lib` folder by the library that we just built

``` bash
rm ../../lib/libaux.so.1
cp libaux.so ../../lib/libaux.so.1
```

The file `#!bash libaux.so.1` is the file that KITE will search for.
Now, the next time that KITE runs, it will use the library we just built.

## 2. Running the example
The rest proceeds as usual.
Use python to generate the configuration file, run the KITE executable on that file and run KITE-tools on the same file again:

``` bash
python config.py
../../build/KITEx config.h5
../../tools/build/KITE-tools config.h5
```

There is an extra flag in the configuration file `#!python custom_local_print=True`, which tells KITE to generate the
files `#!bash local_potentialX.dat`.
By default, this flag is set to `#!python False`.
These files contain the potential as a function of position for all the lattice points used by KITE (one file per thread).
These files can be further processed to generate a color plot of the potential for each orbital:

```
python3 test_potential.py
```

## 3. Finalizing
As a result of this whole process, the default shared library used by KITE in `#!bash /lib` has been changed.
If the flag `#!python custom_local`" is unspecified, or set to `#!python False`, KITE will simply ignore this library,
but if it is set to `#!python True`, KITE will always use the shared library that was built in this example.
It is recommended to simply not specify this flag, and KITE will have the intended behavior.
If the user wants to restore the original library, a backup should be used,
or the following three lines can be executed inside the `#!bash /lib` folder

``` bash
g++ -Wall -fPIC -c aux.cpp 
g++ -shared -Wl,-soname,libaux.so.1 -o libaux.so aux.o
ln -sf libaux.so libaux.so.1
```
