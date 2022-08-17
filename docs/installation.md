!!! warning

    KITE runs exclusively on UNIX-based systems, such as GNU/Linux and Mac OS X.
    Currently, there is **no** Windows support.

KITE is written in C++ with code optimisation, including multithreading performance.
The Python package [Pybinding][pybinding] is used during pre- and post-processing steps.

!!! tip

    To run the pre-processing tools using Python, run the Python-script from the `#!bash kite/`-directory.
    By default, [KITE's python interface][kitepython] is only avaible within the `#!bash kite/`-directory.


## 1. Download KITE

First download the source code from our official repository on GitHub [repository][repository]:

``` bash
git clone https://github.com/quantum-kite/kite.git
```

!!! info

    If git is not installed on your Mac OS X device, come back later to this step.
    A description of the installation process for Git will be given in the next section.


## 2. Get dependencies

KITE dependends on:

* [Pybinding][pybinding]
* [Eigen3][eigen3] (version 3.3.7 or newer)
* [HDF5][hdf5] (version 1.8.13 or newer)
* CMake (version 3.9 or newer)
* gcc (version 4.8.1 or newer)
* h5py

The compiler **must** support *C++11* features and *OpenMP* parallelization.


To enable KITE's Gaussian wavepacket propagation functionality, compile the source code with a recent gcc version
(gcc 8.0.0 or newer).
To check the gcc version, you can use the following command in the terminal:

``` bash
g++ --version
```

### 2.1 For Ubuntu

Install *Eigen3* for various linear algebra tools:

``` bash
sudo apt-get install libeigen3-dev
```

Ubuntu may not provide the required version of Eigen3, retrieve the [latest stable release][eigen3] of Eigen3.
Unzip the file and copy the *Eigen* directory to */usr/include/*.

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

Alternatively, you might prefer to follow the instructions on the [Pybinding][pybinding] webpage.

To construct the HDF5-files, KITE requires *h5py*: 

``` bash
pip3 install h5py
```

### 2.2 For Mac OS X

Install the *Xcode* command-line tools from Apple Developer, install these using the terminal:

``` bash
xcode-select --install
```


!!! info

    On machines from the *Apple Silicon lineup* (M1, M1Max, etc.) you may have to use the *Rosetta* translator while
    a few remaining incompatibility issues are resolved.
    *Rosetta* simulates an *Intel-x64* system and translates existing software for use with *Apple Silicon*.
    
    To load *Rosetta*, run `#!bash arch -x86_64 zsh` **each time** when starting a **new terminal**.

Kite requires an open-source software package management system like [Homebrew][homebrew].
Run the following command in the terminal and follow the subsequent instructions provided by Homebrew.

``` bash
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

Install a C++ compiler via Homebrew:

``` bash
brew install gcc
```

Now **close** the terminal window, and open a **new terminal** window.

Hierarchical Data Format (*HDF5*) is used to store the inputs/outputs of the program.

!!! info
    
    In the following sections, replace **n** with the version of gcc installed by Homebrew as given by `#!bash brew info gcc`.

Install *HDF5* from source:

``` bash
brew install hdf5 --build-from-source --cc=gcc-n
```

!!! info
    
    The default directory for Homebrew is */usr/local/bin/*.
    Change this path to the right location if Homebrew was installed in a different directory

Install *Eigen3* for various linear algebra tools, CMake and Python:

``` bash
CXX=g++-n brew install eigen python Cmake git
```

Calculations on KITE are configured using a python script which interfaces with Pybinding.
Pybinding also requires the SciPy packages but pip will resolve all the SciPy dependencies automatically:

!!! warning

    To install the pyhton requirements, you **must** run the Homebrew-python version.
    You can find the Homebrew-python binary at `#!bash /opt/homebrew/bin/python3.

``` bash
CXX=g++-n /usr/local/bin/python3 -m pip install numpy h5py pybinding
```

Alternatively, you might prefer to follow the instructions on the [Pybinding][pybinding] webpage.

Next, download the source code by the command given in section 1.
Edit *CMakeLists.txt* in the `#!bash kite/`-directory:
* locate the following statements
  ```
  set(CMAKE C COMPILER "gcc")
  set(CMAKE C COMPILER "g++")
  ```

* replace by
  
  ```
  set(CMAKE CXX COMPILER "gcc-n")
  set(CMAKE CXX COMPILER "g++-n")
  ```

where **n** is the version number as used previously.

## 3. KITEx & KITE-tools
From within the `#!bash kite/` directory (containing *CMakeLists.txt* and [*kite.py*][kitepython]), run the following commands for [KITEx][kitex]:

``` bash
mkdir build
cd build
cmake ..
make
```

!!! info

    Any warnings appearing during the compilation process can typically be ignored.

For [KITE-tools][kitetools], run the following commands from the `#!bash kite/tools/` directory

``` bash
mkdir build
cd build
cmake ..
make
```

If these commands have run successfully, you will now find [KITEx][kitex] in the `#!bash kite/build/` directory and [KITE-tools][kitetools in the
`#!bash kite/tools/build/` directory, which are now ready to use!


## 4. Test KITE

To generate an input file using [KITE's python-interface][kitepython], try one of our examples in the `#!bash kite/examples/` directory:

``` bash
python dos_graphene.py
```

It creates a file named *graphene_lattice-output.h5* that is used as an input for [KITEx][kitex]:

``` bash
../build/KITEx graphene_lattice-output.h5
```

This first example calculates the density of states (DOS) of pure graphene.
To obtain the file with the DOS-data, you need to [post-process][kitetools] the output with a tool 

``` bash
../tools/build/KITE-tools graphene_lattice-output.h5
```

that generates the appropriata data file. For more details refer to the [tutorial][tutorial].

!!! info

    The three command above were run from the `#!bash kite/examples/`-directory. If you didn't build [KITEx][kitex] or [KITE-tools][kitetools] in the
    `#!bash kite/build/` and `#!bash kite/tools/build/` directories respectively, the commands won't work.

[repository]: https://github.com/quantum-kite/kite
[eigen3]: https://eigen.tuxfamily.org/
[hdf5]: https://www.hdfgroup.org/
[homebrew]: https://brew.sh/
[pybinding]: https://docs.pybinding.site/en/stable/install/quick.html
[tutorial]: documentation/index.md
[kitepython]: api/kite.md
[kitex]: api/kitex.md
[kitetools]: api/kite-tools.md


