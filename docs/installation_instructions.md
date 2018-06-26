Quantum Kite is a multithreaded C++ package for efficient evaluation of spectral properties of large-scale tight-binding (TB) hamiltonians.

This is our pre-release, which means that KITE is still under development and debug. More information about the to-do list for the first release is available in the README file. Suggestions for improvements and bug reports please send us an e-mail to support at quantum-kite.com.

KITE is written in C++ but has an interface in python, based on pybinding, a tight-binding package for python. The instructions bellow include the installation of our software and the python interface.

- **Linux**
	- <a href="#ubuntu">Ubuntu Installation</a>
	- <a href="#source">Compiling Libraries From Source Code</a>
		
- <a href="#macosx"> Mac OS X </a>	

##Linux:

<a name="linux"></a>
<a name="ubuntu"></a>
###Ubuntu installation


These instructions are geared towards Ubuntu LTS realese 16.04 although some of the instructions work for other distributions as well (for other Linux distributions, check the section *Compiling libraries from Source Code*)

In order to compile our source code, some libraries have to be present and the compiler must be up to date. The required dependencies are:

* Eigen3

* HDF5 (version 1.8.13 and onward) 

* Pybinding

Eigen (Eigen3) features various linear algebra tools. To install it, run:
~~~bash
sudo apt-get install libeigen3-dev
~~~

Hierarchical Data Format (HDF5) is used to store the outputs of the program in an compact way. To install HDF5, run:
~~~bash
sudo apt-get install h5utils
sudo apt-get install libhdf5-dev
~~~

The calculations are configured using a python script which interfaces with Pybinding, To install Pybinding, you will need **CMake**  and **pip**:
~~~bash
sudo apt-get install cmake
sudo apt-get install pip3
~~~
Pybinding also requires the SciPy packages but  pip will resolve all the SciPy dependencies automatically:
~~~bash
pip3 install pybinding
~~~

Alternativelly, you might prefer to follow the instructions in pybinding webpage directly:
http://docs.pybinding.site/en/stable/install/index.html

**After successfully installing these libraries, you are now ready to compile KITE.**


1. Download the source code from () 
Alternativelly, fetch the source code from the git repository in (…)
~~~bash
git clone (...)
git checkout develop
~~~
2. Execute the Makefile inside the kpm_transport folder to compile the code
~~~bash
make
~~~
3. That’s it! If you want to compile the post-processing program as well, do
~~~bash
cd tools/kpm_calculate
make
~~~
####Compiling Libraries From Source Code
<a name="source"></a>
Instead of installing the dependencies with **apt-get** you can install the libraries manually. For KITE, the compiler must be up to date. The required dependencies are:

* Eigen3

* HDF5 (version 1.8.13 and onward) 

* Pybinding


The detailed installation instructions are described in the respective websites. This is a list of the requirements and some instructions on how to fulfill them:

The g++ compiler must support C++11 features and OpenMP parallelization, so make sure your g++ version is at least 4.8.1. To check the version you may use:
~~~bash
g++ --version
~~~
Hierarchical Data Format (HDF5) is used to store the outputs of the program in an organized fashion. The code requires some functions that are only available from version 1.8.13 and onwards, so make sure you’re installing the correct version.

For more information on this library, visit https://support.hdfgroup.org/HDF5/release/obtainsrc518.html#conf.

 

Eigen (Eigen3) features various linear algebra tools.  You can obtain the latest stable release from http://eigen.tuxfamily.org/index.php?title=Main_Page#Download, unzip the folder and copy the Eigen folder into /usr/include/

The calculations are configured using a python script which interfaces with Pybinding, so make sure that you have python 3 installed with SciPy and Pybinding. We suggest to proceed with the installation suggested by Pybinding, with the use of Miniconda:
http://docs.pybinding.site/en/stable/install/quick.html.

##Mac OS X:
<a name="macosx"></a>
In order to compile our source code, some libraries have to be present and the compiler must be up to date. There are several ways of installing the libraries and c++ in the OS X operational system. Our detailed instructions are focused on an installation in High Sierra OS X with [Homebrew](https://brew.sh/)
To install Homebrew, run the following command in a terminal:
~~~bash
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
~~~

The script explains what it will do and then pauses before it does it.

** Homebrew requires command-line tools from xcode. You can either download xcode from Mac App Store (basic users) or just download the command line tools from the developer site of apple (advanced users)

Next, you need to use Homebrew to install g++:

~~~bash
brew install gcc@6
~~~

The required dependencies are:

* Eigen3

* HDF5 (version 1.8.13 and onward) 

* Pybinding

Eigen (Eigen3) features various linear algebra tools. To install it, run:

~~~bash
brew install eigen
~~~

Hierarchical Data Format (HDF5) is used to store the outputs of the program in an compact way. To install HDF5, run:
~~~bash
brew install hdf5 --cc=gcc-6
~~~
The calculations are configured using a python script which interfaces with Pybinding, so make sure that you have python 3 installed with SciPy and Pybinding. We suggest using Hombrew but you can also proceed with the  [installation suggested by Pybinding](http://docs.pybinding.site/en/stable/install/quick.html), with the use of Miniconda.

~~~bash
brew install cmake
~~~

~~~bash
brew install python
~~~
Now you are ready to install pybiding with pip:
~~~bash
pip3 install pybinding
~~~
**After successfully installing these libraries, you are now ready to compile KITE.**

**IMPORTANT: the last version of matlibplot is having issues with pybinding. Until this is resolved, use:

~~~bash
pip3 install matplotlib==2.1.1 
~~~

1. Download the source code from () 
Alternativelly, fetch the source code from the git repository in (…)
~~~bash
git clone (...)
git checkout develop
~~~
2. Execute the Makefile inside the kpm_transport folder to compile the code
~~~bash
make
~~~
3. That’s it! If you want to compile the post-processing program as well, do
~~~bash
cd tools/postprocessing
make
~~~
To generate the input file, try one of our examples:
~~~bash
python3 example_initial.py
~~~
It creates a file names example.h5 that is used as an input for KITE:
~~~bash
./KITEx example_initial.h5
~~~
This first example calculates the density of states of pure graphene. To obtain the data file, you need to postprocess the output:

~~~bash
./tools/postprocessing/KITEpos example_initial.h5
~~~

For more details on how to create Python scripts and postprocess the data, refer to [Documentation](https://quantum-kite.com/Documentation/)