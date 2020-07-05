Quantum KITE is an open-source software suite for accurate real space evaluations of electronic structure and response functions of large-scale tight-binding (TB) models of complex molecules and crystals with up to multi billions of atoms (orbitals).

KITE is written in C++ with advanced code optimization at the design level (including automated multithreaded lattice decomposition) and allow users to change defaults to enable the best possible use of resources. KITE's user-friendly interface and processing tools are based on [Pybinding][1], a scientific Python package for TB calculations.

KITE is also a quantum transport simulator. Multi-orbital bond disorder can be defined at the interface level and added to the system according to pre-defined probability distributions, allowing to simulate the behavior of realistic disordered systems.

For currently available functionalities and to-do-list please refer to the README file. Please share your feedback and bug reports. We would love to hear your ideas for future improvements (contact email: support at quantum-kite.com).

We are currently working on an installation package to automate the whole installation process. While this is not finished, the manual installation can be done in three simple steps which are detailed below: downloading the source code, making sure the dependencies are met, and compiling the software from the source code.

KITE only runs on UNIX-based systems such as GNU/Linux and Mac OS X.

# 1. Download KITE

The first step is to download the source code. You can find the latest version of KITE at our official repository on GitHub **[repository][3]**. Fetch the source code from the Git Hub [repository][3] using the command
``` bash
git clone https://github.com/quantum-kite/kite.git
```
which will create a folder called kite with the source code inside kite/Src and kite/tools/src.

# 2. Dependencies

The second step is to make sure that all of KITE's dependencies have been met. KITE requires:
* [Pybinding][1]
* [Eigen3][4] (version 3.3.7 or newer)
* [HDF5][5] (version 1.8.13 or newer)
* CMake (version 3.9 or newer)
* gcc (version 4.8.1 or newer)

These packages are available from public domain servers. The compiler must support C++11 features and OpenMP parallelization. Some of KITE's features such as the Gaussian Wavepacket propagation require a more recent version of gcc (8.0.0), but are not required for the program to have its core functionality. To check the version installed in your computer, you can use the following command in the terminal:
``` bash
g++ --version
```
In the next few paragraphs, we provide some convenient instructions on how to install these libraries in Ubuntu and Mac OS X. Instructions on how to install these dependencies for other Linux distributions may be found on the respective websites of the libraries.

## 2.1 Dependencies in Ubuntu


The instructions below were tested with Ubuntu LTS release 16.04; for other Linux distributions please refer to *Compiling libraries from Source Code*.

Eigen (Eigen3) features various linear algebra tools. To install it, run:
``` bash
sudo apt-get install libeigen3-dev
```
Ubuntu may not provide the required version of Eigen3. You can retrieve the latest stable release of Eigen3 from [eigen.tuxfamily.org][12]. Unzip the file and copy the Eigen folder to /usr/include/. 

Hierarchical Data Format (HDF5) is used to store the inputs/outputs of the program. To install HDF5, run:
``` bash
sudo apt-get install h5utils
sudo apt-get install libhdf5-dev
```
Calculations on KITE are configured using a python script which interfaces with Pybinding. To install Pybinding, you will need **CMake** and **pip**:
``` bash
sudo apt-get install cmake
sudo apt-get install pip3
```
Pybinding also requires the SciPy packages but pip will resolve all the SciPy dependencies automatically:
``` bash
pip3 install pybinding
```
Alternativelly, you might prefer to follow the instructions on the [Pybinding][1] webpage. 

**IMPORTANT**: the last version of matlibplot is having issues with pybinding. Until this is resolved, use:
``` bash
pip3 install matplotlib==2.1.1
```


## 2.2 Dependencies in Mac OS X

In order to install the dependencies on a MAC OS X system, you will need Xcode command-line tools from Apple Developer. You can download Xcode from the Apple Store. Alternatively, you can install the command tools directly from the terminal:
``` bash
xcode-select --install
```
You will also need an open-source software package management system. Here, we provide detailed instructions using [Homebrew][14]. To install Homebrew run the following command on the terminal:
``` bash
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```
Follow the instructions provided. Next, install the C++ compiler:
``` bash
brew install gcc@6
```

Eigen (Eigen3) features various linear algebra tools. To install it, run:
``` bash
brew install eigen
```
Hierarchical Data Format (HDF5) is used to store the inputs/outputs of the program. To install HDF5, run:
``` bash
brew install hdf5 --cc=gcc-6
```
Calculations on KITE are configured using a python script which interfaces with Pybinding. To install the Pybinding package, we suggest using Homebrew. (Alternatively, you can proceed with the [installation suggested by Pybinding][13], with the use of Miniconda.)
``` bash
brew install cmake
brew install python
```
Last, install Pybinding with pip:
``` bash
pip3 install pybinding
```
**IMPORTANT**: the last version of matlibplot is having issues with pybinding. Until this is resolved, use:
``` bash
pip3 install matplotlib==2.1.1
```

# 3. Install KITE

The last step is to build KITE from source. The following set of commands will create a build directory for both KITE and KITE-tools, and compile the programs from their respective source codes. From within the kite directory (the one containing CMakeLists.txt and kite.py), run the following commands to compile the main program KITEx:
``` bash
mkdir build
cd build
cmake ..
make 
```
In order to compile the suite of post-processing utilities, KITE-tools, run the following commands from the kite/tools directory

```
mkdir build
cmake ..
make
```
If these commands ran successfully, you will now find KITEx in the build directory and KITE-tools in the tools/build directory, which are now ready to use!


# Testing the program

To generate the input file, try one of our examples. In the kite/examples folder, run
``` bash
python example_initial.py
```
It creates a file named example_initial.h5 that is used as an input for KITE:
``` bash
../build/KITEx example_initial.h5
```
This first example calculates the density of states of pure graphene. To obtain the data file, you need to postprocess the output:
``` bash
../tools/build/KITEx-tools example_initial.h5
```
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
