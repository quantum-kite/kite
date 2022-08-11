## Usage


### Default usage

Its default usage is very simple:

``` bash
./KITE-tools archive.h5
```

where archive.h5 is the HDF file that stores the output of KITE. If KITE-tools does not find this output, it will return an error. The output of KITE-tools is a set of .dat files, one for each of the requested quantities. KITE-tools may be executed without any additional parameters; all the unspecified parameters required for the calculation will be set to sensible default values. At the moment, KITE-tools is able to compute the following quantities:

  * Local density of states (LDOS)
  * Angle-resolved photoemission spectroscopy (experimental) (ARPES)
  * Density of states (DOS)
  * DC conductivity (CondDC)
  * Optical conductivity (CondOpt)
  * Second-order optical conductivity (CondOpt2)

The SingleShot DC conductivity does not require the post-processing through KITE-tools.


### Advanced usage

KITE-tools supports a set of command-line instructions to force it to use user-specified parameters for each of the quantities mentioned in the previous section. The syntax is as following:

``` bash
./KITE-tools archive.h5 --quantity_to_compute1 -key_1 value_1 -key_2 value_2 --quantity_to_compute2 -key_3 value_3 ...
```

Each function to compute is specified after the double hyphens — and the parameters of each function is specified after the single hyphen -. The list of available commands is as follows:

| Function            | Parameter    | Description                                                                                         |
|---------------------|--------------|-----------------------------------------------------------------------------------------------------|
| `#!bash --LDOS`     | `#!bash -N`  | Name of the output file                                                                             |
| `#!bash --LDOS`     | `#!bash -M`  | Number of Chebyshev moments                                                                         |
| `#!bash --LDOS`     | `#!bash -K`  | Kernel to use (jackson/green). green requires broadening parameter. Example: `#!bash -K green 0.01` |
| `#!bash --LDOS`     | `#!bash -X`  | Exclusive. Only calculate this quantity                                                             |
| `#!bash --ARPES`    | `#!bash -N`  | Name of the output file                                                                             |
| `#!bash --ARPES`    | `#!bash -E`  | min max num Number of energy points                                                                 |
| `#!bash --ARPES`    | `#!bash -F`  | Fermi energy                                                                                        |
| `#!bash --ARPES`    | `#!bash -T`  | Temperature                                                                                         |
| `#!bash --ARPES`    | `#!bash -V`  | Wave vector of the incident wave                                                                    |
| `#!bash --ARPES`    | `#!bash -O`  | Frequency of the incident wave                                                                      |
| `#!bash --ARPES`    | `#!bash -X`  | Exclusive. Only calculate this quantity                                                             |
| `#!bash --DOS`      | `#!bash -N`  | Name of the output file                                                                             |
| `#!bash --DOS`      | `#!bash -E`  | Number of energy points                                                                             |
| `#!bash --DOS`      | `#!bash -M`  | Number of Chebyshev moments                                                                         |
| `#!bash --DOS`      | `#!bash -K`  | Kernel to use (jackson/green). green requires broadening parameter. Example: `#!bash -K green 0.01` |
| `#!bash --DOS`      | `#!bash -X`  | Exclusive. Only calculate this quantity                                                             |
| `#!bash --CondDC`   | `#!bash -N`  | Name of the output file                                                                             |
| `#!bash --CondDC`   | `#!bash -E`  | Number of energy points used in the integration                                                     |
| `#!bash --CondDC`   | `#!bash -M`  | Number of Chebyshev moments                                                                         |
| `#!bash --CondDC`   | `#!bash -T`  | Temperature                                                                                         |
| `#!bash --CondDC`   | `#!bash -S`  | Broadening parameter of the Green’s function                                                        |
| `#!bash --CondDC`   | `#!bash -d`  | Broadening parameter of the Dirac delta                                                             |
| `#!bash --CondDC`   | `#!bash -F`  | min max numRange of Fermi energies. min and max may be omitted if only one is required              |
| `#!bash --CondDC`   | `#!bash -t`  | Number of threads                                                                                   |
| `#!bash --CondDC`   | `#!bash -I`  | If `#!bash 0`, CondDC uses the DOS to estimate the integration range                                |
| `#!bash --CondDC`   | `#!bash -X`  | Exclusive. Only calculate this quantity                                                             |
| `#!bash --CondOpt`  | `#!bash -N`  | Name of the output file                                                                             |
| `#!bash --CondOpt`  | `#!bash -E`  | Number of energy points used in the integration                                                     |
| `#!bash --CondOpt`  | `#!bash -M`  | Number of Chebyshev moments                                                                         |
| `#!bash --CondOpt`  | `#!bash -T`  | Temperature                                                                                         |
| `#!bash --CondOpt`  | `#!bash -F`  | Fermi energy                                                                                        |
| `#!bash --CondOpt`  | `#!bash -S`  | Broadening parameter of the Green’s function                                                        |
| `#!bash --CondOpt`  | `#!bash -O`  | min max num Range of frequencies                                                                    |
| `#!bash --CondOpt2` | `#!bash -N`  | Name of the output file                                                                             |
| `#!bash --CondOpt2` | `#!bash -E`  | Number of energy points used in the integration                                                     |
| `#!bash --CondOpt2` | `#!bash -M`  | Number of Chebyshev moments                                                                         |
| `#!bash --CondOpt2` | `#!bash -R`  | Ratio of the second frequency relative to the first one                                             |
| `#!bash --CondOpt2` | `#!bash -P`  | If set to 1: writes all the different contributions to separate files                               |
| `#!bash --CondOpt2` | `#!bash -T`  | Temperature                                                                                         |
| `#!bash --CondOpt2` | `#!bash -F`  | Fermi energy                                                                                        |
| `#!bash --CondOpt2` | `#!bash -S`  | Broadening parameter of the Green’s function                                                        |
| `#!bash --CondOpt2` | `#!bash -O`  | min max num Range of frequencies                                                                    |

All the values specified in this way are assumed to be in the same units as the ones used in the configuration file. All quantities are double-precision numbers except for the ones representing integers, such as the number of points. This list may be found in KITE-tools, run `KITE-tools --help`:

## Output

In the table below, we specify the name of the files that are created by KITE-tools according to the calculated quantity and the format of the data file.

| Quantity                          | File                         | Column 1          | Column 2         | Column 3         |
|-----------------------------------|------------------------------|-------------------|------------------|------------------|
| Local Density of States           | `#!bash ldos{E}.dat`         | lattice position  | LDOS ($Re$)      |                  |
| ARPES                             | `#!bash arpes.dat`           | k-vector          | ARPES ($Re$)     |                  |	
| Density of States                 | `#!bash dos.dat`             | energy            | DOS ($Re$)       | DOS ($Im$)       |
| Optical Conductivity              | `#!bash optical_cond.dat`    | Frequency         | Opt. Cond ($Re$) | Opt. Cond ($Im$) |
| DC Conductivity                   | `#!bash condDC.dat`          | Fermi energy      | Cond ($Re$)      | Cond ($Im$)      |
| Second-order optical conductivity | `#!bash nonlinear_cond.dat`  | Frequency         | NL Cond ($Re$)   | NL Cond ($Im$)   |
| Single-shot DC Conductivity       | `#!bash [HDF5-filename].dat` | Fermi energy      | Cond ($Re$)      | Cond ($Im$)      |

* All linear conductivities are in units of $e^2/h$
* Both Planck’s constant and electron charge are set to 1.
* LDOS outputs one file for each requested energy. The energy is in the E in the file name.

For more details on the type of calculations performed during post-processing, check [Resources][resources] where we discuss our method.

!!! info "Processing the single-shot DC conductivity"

    The single shot DC conductivity does not need any post-processing as it is an energy dependent calculation where the conductivity is calculated on the fly.
    In this particular case, the data is extracted directly from the hdf file with the following python script in the `#!bash kite/tools/`-folder:
    
    ``` bash
    python3 process_single_shot.py output.h5
    ```
    
    The output of this script will be a data-file with name `#!bash output.dat` or similar, structured as given in the table above.

## Examples


### Example 1

``` bash
./KITE-tools h5_file.h5 --DOS -E 1024
```

Processes the .h5 file as usual but ignores the number of energy points in the density of states present there. Instead, KITE-tools will use the value 1024 as specified in the example.

### Example 2

``` bash
./KITE-tools h5_file.h5 --CondDC -E 552 -S 0.01
```

Processes the .h5 file but uses 552 points for the energy integration and a broadening parameter of 0.01.

### Example 3

``` bash
./KITE-tools h5_file.h5 --CondDC -T 0.4 -F 500
```

Calculates the DC conductivity using a temperature of 0.4 and 500 equidistant Fermi energies spanning the spectrum of the Hamiltonian.

### Example 4

``` bash
./KITE-tools h5_file.h --CondDC -F -1.2 2.5 30 --CondOpt -T 93
```

Calculates the DC conductivity using 30 equidistant Fermi energies in the range `#!python [-1.2, 2.5]` and the optical conductivity using a temperature of 93.

@@include[kite_tools_readme.md](kite_tools_readme.md)

[resources]: ../background/index.md
 