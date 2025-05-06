KITE is written in C++ with code optimized for large systems and optimal multithreading performance. 
[Pybinding][pybinding] is KITE's default interface, 
which is primarily used to build the configuration (HDF5) file for [KITEx].

The KITE team endeavours to assist researchers run KITE on UNIX-based systems, such as GNU/Linux and Mac OS X.
Thus, feel free to contact any of our team members if you have any queries (contacts can be found at the bottom of the landing page). 

In what follows, we provide detailed installation instructions and additional tips for both Linux and MAC users. 

## 1. Download KITE

First download the source code from our official repository on GitHub [repository][repository]:

``` bash
git clone https://github.com/quantum-kite/kite.git
```

!!! info
   
    Git's installation process for Mac users is outlined in section 2.2.


## 2. Get dependencies

* [Pybinding][pybinding]
* [Eigen3][eigen3] (version 3.3.7 or newer)
* [HDF5][hdf5] (version 1.8.13 or newer)
* [CMake][cmake] (version 3.9 or newer)
* [gcc][gcc] (version 4.8.1 or newer)
* [h5py][h5py]

The compiler **must** support *C++17* (or newer) features and [*OpenMP*][openmp] parallelization.


To enable KITE's [Gaussian wavepacket propagation][calculation-gaussian_wave_packet] functionality, compile the source code with a recent gcc version
(gcc 8.0.0 or newer).
To check the gcc version, you can use the following command in the terminal:

``` bash
g++ --version
```

### 2.1 For Ubuntu users

Install *Eigen3* for various linear algebra tools:

``` bash
sudo apt-get install libeigen3-dev
```

Make sure you retrieve the [latest stable release][eigen3] of Eigen3. Unzip the file and copy the *Eigen* directory to */usr/include/*.

Hierarchical Data Format (*HDF5*) is used to store the inputs/outputs of the program:

``` bash
sudo apt-get install h5utils
sudo apt-get install libhdf5-dev
```

Calculations on KITE are configured using a python script which interfaces with Pybinding.
Pybinding requires CMake and pip:

``` bash
sudo apt-get install cmake
sudo apt-get install pip3
```

Pybinding also requires the SciPy packages but pip will resolve all the SciPy dependencies automatically:

``` bash
pip3 install pybinding
```

!!! info

    Due to outdated dependencies, the installation of pybinding might fail.
    You can install the development version of Pybinding with:
    ``` bash
    pip3 install git+https://github.com/BertJorissen/pybinding
    ```

To construct the HDF5-files, KITE requires *h5py*: 

``` bash
pip3 install h5py
```

### 2.2 For Mac OS X users

The installation of KITE's dependencies on Apple machines is slightly more evolved. We provide below a recipe that has been tested on some Mac OS X systems, but users are encouraged to contact the KITE team shall they encounter any difficulties.  

The *Xcode* command-line tools from Apple Developer are required.  Install these using the terminal:

``` bash
xcode-select --install
```

KITE requires an open-source software package management system like [Homebrew][homebrew] or [MacPorts][ports]. We provide here step-by-step instructions for Homebrew (pointers for MacPorts users are given below). To install HomeBrew, run the following command in the terminal and follow the subsequent instructions provided by software:

``` bash
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

Install an up-to-date C++ compiler via Homebrew:

``` bash
brew install gcc
```

Now **close** the terminal window, and open a **new terminal** window.

!!! info
    
    1.  The default directory for Homebrew is */usr/local/bin/*.
        Correct this path if Homebrew was installed in a different directory
    
    2.  In the following sections, replace **n** with the version of gcc installed by Homebrew as given by `#!bash brew info gcc`.


The hierarchical Data Format (*HDF5*) is used to store the inputs/outputs of the program. Install *HDF5* from source, _whilst enforcing the C++17 standard_, using:

``` bash
HOMEBREW_CC=gcc-n HOMEBREW_CXX=g++-n HOMEBREW_CXXFLAGS="-std=c++17" brew install hdf5 --build-from-source
```

!!! info
    
    [MacPorts][ports] users can use the following command:

    ``` bash
    sudo port -v install hdf5 +gcc-n +cxx +hl configure.ldflags="-stdlib=libstdc++" configure.cxx_stdlib="libstdc++" configure.cxxflags="-std=c++17" 
    ```

Install *Eigen3* for various linear algebra tools, CMake and Python:

``` bash
brew install eigen python Cmake git
```

Calculations on KITE are configured using a python script which interfaces with Pybinding.
Pybinding also requires the SciPy packages but pip will resolve all the SciPy dependencies automatically:

!!! warning

    To install the pyhton requirements, you **must** run the Homebrew-python version.
    You can find the Homebrew-python binary at `#!bash /opt/homebrew/bin/python3`.

``` bash
/usr/local/bin/python3 -m pip install numpy h5py pybinding
```

!!! info

    Due to outdated dependencies, the installation of pybinding might fail.
    You can install the development version of Pybinding with:
    ``` bash
    pip3 install git+https://github.com/BertJorissen/pybinding
    ```
    Alternatively, you might prefer to follow the instructions on the [Pybinding][pybinding] webpage.

Next, download the source code by the command given in section 1.
Edit *CMakeLists.txt* in the `#!bash kite/`-directory:

* locate the following statements
  ```
  set(CMAKE C COMPILER "gcc")
  set(CMAKE CXX COMPILER "g++")
  ```

* replace by
  ```
  set(CMAKE C COMPILER "gcc-n")
  set(CMAKE CXX COMPILER "g++-n")
  ```

where **n** is the version number as used previously.

## 3. KITEx & KITE-tools
From within the `#!bash kite/` directory (containing *CMakeLists.txt* and [*kite.py*][kitepython]), run the following commands:

``` bash
mkdir build
cd build
cmake ..
make
```

!!! info

    Any warnings appearing during the compilation process can typically be ignored.

If these commands have run successfully, you will now find [KITEx][kitex] and [KITE-tools][kitetools]  in the `#!bash kite/build/` directory, which are now ready to use!


## 4. Test KITE

To generate an input file using [KITE's python-interface][kitepython], try one of our examples in the `#!bash kite/examples/` directory:

``` bash
python dos_graphene.py
```

It creates a file named *graphene_lattice-output.h5* that is used as an input for [KITEx][kitex]:

``` bash
../build/KITEx graphene_lattice-output.h5
```

This first example calculates the density of states (DOS) of pristine graphene.
To obtain the file with the DOS-data, you need to [post-process][kitetools] the output with  

``` bash
../build/KITE-tools graphene_lattice-output.h5
```

which generates the appropriate data file. For more details refer to the [tutorial][tutorial].

!!! info

    The three command above were run from the `#!bash kite/examples/`-directory. If you didn't build [KITEx][kitex] or [KITE-tools][kitetools] in the
    `#!bash kite/build/` and `#!bash kite/build/` directories respectively, the commands won't work.

## 5. Common issues

### 5.1 Finding Eigen3 with CMake

If experiencing difficulties running `#!bash cmake..` due to CMake not being able to locate Eigen3, you may opt for the following solution. Extract [Eigen3][eigen3] from source and then run
the following command from Eigen3's directory

``` bash
mkdir build
cd build
cmake ..
make
```
You may now wish to try and re-compile KITEx and KITE-tools following the steps of Sec. 3.

### 5.2 Apple Silicon (ARM) architecture 

On machines from the *Apple Silicon lineup* (M1, M1Max, etc.) you may have to use the *Rosetta* translator if experiencing
architecture compatibility issues. *Rosetta* simulates an *Intel-x64* system and translates existing software for use with *Apple Silicon*. 
To load *Rosetta*, run `#!bash arch -x86_64 zsh` when starting a **new terminal**.


[repository]: https://github.com/quantum-kite/kite
[eigen3]: https://eigen.tuxfamily.org/
[cmake]: https://cmake.org/
[gcc]: https://gcc.gnu.org/
[h5py]: https://www.h5py.org/
[calculation-gaussian_wave_packet]: api/kite.md#calculation-gaussian_wave_packet
[hdf5]: https://github.com/HDFGroup/
[openmp]: https://gcc.gnu.org/onlinedocs/libgomp/
[homebrew]: https://brew.sh/
[ports]: https://www.macports.org 
[pybinding]: https://docs.pybinding.site/en/stable/install/quick.html
[tutorial]: documentation/index.md
[kitepython]: api/kite.md
[kitex]: api/kitex.md
[kitetools]: api/kite-tools.md


