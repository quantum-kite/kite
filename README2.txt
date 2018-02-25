___________________________ Compilation ___________________________

To use the KITE software, make sure the program pp exists in the kpm_transport folder. If not, compile it by running 

	make

If you only wish to obtain the matrices of Chebyshev moments, this suffices. However, if you also want to use our software to calculate density of states, conductivity, etc, you also need to make sure another program exists, calc_kpm in the folder tools/kpm_calculate/. If it does not exist, move to that folder 

	cd tools/kpm_calculate/

and run

	make

If the compilation process was successful, you're set to go. If not, follow the installation instructions in the KITE website.






___________________________ Some documentation on how to use the KITE software ___________________________

These simple examples illustrate how to use the KITE software. You first create a configuration script in python, run it and feed the .h5 file into the pp program to obtain the Chebyshev matrices. Then feed that .h5 again into the calc_kpm program to calculate the quantities. Make sure you have pybinding installed and that the export_lattice script is present in the same folder as the python script.

Some examples are provided in the main folder. Here we show a concrete example on how to use our software to calculate the density of states of graphene. 






---------------- Example 1 ----------------
1. First of all, run the configuration script with python to generate the configuration file. In this script you can define everything you need about the simulation: lattice geometry, lattice size, parallelization parameters, tight-binding hoppings, structural disorder, Anderson disorder, magnetic field and the objects you want to calculate (density of states, conductivity). This generates a single .h5 file with all the parameters needed for the simulation.

	python example1.py

2. You should now have a file called 'example1.h5'. This file should be fed into the main program, here called 'pp'. No further arguments are needed because everything the program needs to know is inside the .h5 file. This program calculates and stores the matrices of Chebyshev moments which may be used to calculate various quantities such as the density of states and the conductivity. In this case, we'll obtain the 1D matrix of Chebyshev moments of the density of states.

	./pp example1.h5

After running the program, the configuration file should now be filled with the objects you asked it to calculate (MU is the 1D matrix of Chebyshev moments relative to the Density of states, CondXX is the one relative to the first order XX conductivity). In this case, we asked for the density of states, so the .h5 file should have a 'MU' entry. 

3. If you want to go further, you can feed this filled .h5 file to the calc_kpm program, which takes the existing Chebyshev moment matrices and tries to calculate the quantities which use those objects (DOS and CondXX).

	tools/kpm_calculate/calc_kpm example1.h5

This will generate a file called DOS.dat with two columns: the energy and the density of states. You can now plot this file and get the density of states of graphene!






---------------- Example 2 ----------------
The procedure is precisely the same as in the previous example, but this script calculates the longitudinal (XX) optical conductivity of graphene. This will take considerably longer because the matrix being calculated is now two-dimensional. In our tests, this took 12 minutes, but if you don't want to wait, you can reduce the lattice size and the number of Chebyshev polynomials being used. The resulting file is called optical_cond.dat and has three columns: the frequency, the real part of the conductivity and the imaginary part of the conductivity.

	python example2.py
	./pp example2.h5
	tools/kpm_calculate/calc_kpm example2.h5






---------------- Example 3 ----------------
Again with the same procedure, this example uses a different lattice geometry: the square lattice and calculates the density of states. It also uses the domain decomposition into 4 sublattices to speed up the calculation.

	python example3.py
	./pp example3.h5
	tools/kpm_calculate/calc_kpm example3.h5






---------------- Example 4 ----------------
This example shows how to include structural disorder. The density of states of graphene is calculated, but this time there is one selected bond missing from each unit cell. This is addressed by treating it as a disorder that is added with 100% probability to each cell.

	python example4.py
	./pp example4.h5
	tools/kpm_calculate/calc_kpm example4.h5



