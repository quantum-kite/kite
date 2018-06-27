Quantum KITE is an open-source software suite for accurate real space evaluations of electronic structure and response functions of large-scale tight-binding (TB) models of complex molecules and crystals with up to multi billions of atoms (orbitals).

KITE is written in C++ with advanced code optimization at the design level (including automated multithreaded lattice decomposition) and allow users to change defaults to enable the best possible use of resources. KITE's user-friendly interface and processing tools are based on [Pybinding][1], a scientific Python package for TB calculations.

KITE is also a quantum transport simulator. Multi-orbital bond disorder can be defined at the interface level and added to the system according to pre-defined probability distributions, allowing to simulate the behavior of realistic disordered systems.

*This is our pre-release BETA version*. The fully-debugged official release is scheduled for 01/11/2018.

For currently available functionalities and to-do-list please refer to the README file. Please share your feedback and bug reports. We would love to hear your ideas for future improvements (contact email: support at quantum-kite.com).

# Download KITE

You can download a zip file with our code **[here][2]**

Or you can get the code from our **[repository][3]**

# Install KITE

KITE runs on UNIX®-based systems (including Mac OS X) and requires [Pybinding][1], [Eigen][4] C++ template library and [HDF5][5] support for multi-dimensional datasets. These packages are available from public domain servers (see below).

* [Linux ][6]

* [Ubuntu Installation][7]

* [Compiling Libraries From Source Code][8]

* [ Mac OS X ][9]

## Linux:


### Ubuntu installation


The instructions below were tested with Ubuntu LTS realese 16.04; for other Linux distributions please refer to *Compiling libraries from Source Code*.

In order to compile our source code, the compiler must be up to date (e.g. GCC 4.8.1 or newer). The required dependencies are:

* [Eigen3][10]

* [HDF5][5] (version 1.8.13 or newer)

* [Pybinding][1]

Eigen (Eigen3) features various linear algebra tools. To install it, run:

sudo apt-get install libeigen3-dev
Hierarchical Data Format (HDF5) is used to store the inputs/outputs of the program. To install HDF5, run:

sudo apt-get install h5utils
sudo apt-get install libhdf5-dev
Calculations on KITE are configured using a python script which interfaces with Pybinding. To install Pybinding, you will need **CMake** and **pip**:

sudo apt-get install cmake
sudo apt-get install pip3
Pybinding also requires the SciPy packages but pip will resolve all the SciPy dependencies automatically:

pip3 install pybinding
Alternativelly, you might prefer to follow the instructions on [Pybinding][1] webpages.

**IMPORTANT**: the last version of matlibplot is having issues with pybinding. Until this is resolved, use:

pip3 install matplotlib==2.1.1
**After successfully installing these libraries, you will be ready to compile KITE.**

Fetch the source code from the Git Hub [repository][3]

git clone https://github.com/quantum-kite/kite.git
Execute the Makefile inside the KITE folder to compile the code

make
That’s it! To compile the post-processing program run

cd tools/KITE
make
To generate the input file, try one of our examples. In the KITE folder, run

python example_initial.py
It creates a file names example.h5 that is used as an input for KITE:

./KITEx example_initial.h5
This first example calculates the density of states of pure graphene. To obtain the data file, you need to postprocess the output:

./tools/KITE/KITEx-tools example_initial.h5
For more details refer to the KITE [Documentation][11].

### Compiling Libraries From Source Code


In order to compile our source code, the compiler must be up to date (e.g. GCC 4.8.1 or newer). The required dependencies are:

* [Eigen3][10]

* [HDF5][5] (version 1.8.13 or newer)

* [Pybinding][1]

The detailed installation steps are described on the respective websites. The GCC (g++) compiler must support C++11 features and OpenMP parallelization (version 4.8.1 or newer). To check the version installed on your computer you can use:

g++ --version
You can retrieve the latest stable release of Eigen3 from [eigen.tuxfamily.org][12]. Unzip the file and copy the Eigen folder to /usr/include/. We suggest the installation of [Pybinding][13] using Miniconda.

## Mac OS X:


In order to install KITE on a MAC OS X system, you will need Xcode command-line tools from Apple Developer. You can download Xcode from the Apple Store. Alternatively, you can install the command tools directly from the terminal:

xcode-select --install
You will also need an open-source software package management system. Here, we provide detailed instructions using [Homebrew][14]. To install Homebrew run the following command on the terminal:

/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
Follow the instructions provided. Next, install the C++ compiler:

brew install gcc@6
KITE has the following dependencies:

* [Eigen3][10]

* [HDF5][5] (version 1.8.13 or newer)

* [Pybinding][1]

Eigen (Eigen3) features various linear algebra tools. To install it, run:

brew install eigen
Hierarchical Data Format (HDF5) is used to store the inputs/outputs of the program. To install HDF5, run:

brew install hdf5 --cc=gcc-6
Calculations on KITE are configured using a python script which interfaces with Pybinding. To install the Pybinding package, we suggest using Homebrew. (Alternatively, you can proceed with the [installation suggested by Pybinding][13], with the use of Miniconda.)

brew install cmake
brew install python
Last, install Pybinding with pip:

pip3 install pybinding
**After successfully installing these libraries, you are now ready to compile KITE.**

**IMPORTANT**: the last version of matlibplot is having issues with pybinding. Until this is resolved, use:

pip3 install matplotlib==2.1.1
Fetch the source code from the Git Hub repository

 git clone https://github.com/quantum-kite/kite.git
Execute the Makefile inside the KITE folder to compile the code

make
That’s it! To compile the post-processing program run

cd tools/KITE/
make
To generate the input file, try one of our examples. In the KITE folder, run

python example_initial.py
It creates a file names example.h5 that is used as an input for KITE:

./KITEx example_initial.h5
This first example calculates the density of states of pure graphene. To obtain the data file, you need to postprocess the output:

./tools/KITE/KITE-tools example_initial.h5
For more details refer to the KITE [Documentation][11].

[1]: http://docs.pybinding.site/en/stable/
[2]: https://github.com/quantum-kite/kite/archive/master.zip
[3]: https://github.com/quantum-kite/kite.git
[4]: eigen.tuxfamily.org
[5]: https://www.hdfgroup.org/
[6]: #linux
[7]: #ubuntu
[8]: #source
[9]: #macosx
[10]: http://eigen.tuxfamily.org
[11]: https://quantum-kite.com/Documentation/
[12]: http://eigen.tuxfamily.org/
[13]: http://docs.pybinding.site/en/stable/install/quick.html
[14]: https://brew.sh/
