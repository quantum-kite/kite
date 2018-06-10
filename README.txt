___________________________ Compilation ___________________________

To use the KITE software, make sure the program pp exists in the kpm_transport folder. If not, compile it by running 

	make

If you only wish to obtain the matrices of Chebyshev moments, this suffices. However, if you also want to use our software to calculate density of states, conductivity, etc, you also need to make sure another program exists, calc_kpm in the folder tools/kpm_calculate/. If it does not exist, move to that folder 

	cd tools/

and run

	make

If the compilation process was successful, you're set to go. If not, follow the installation instructions in the KITE website.






___________________________ Some documentation on how to use the KITE software ___________________________

These simple examples illustrate how to use the KITE software. You first create a configuration script in python, run it and feed the .h5 file into the KITEx program to obtain the Chebyshev matrices. Then feed that .h5 again into the calc_kpm program to calculate the quantities. Make sure you have pybinding installed and that the export_lattice script is present in the same folder as the python script.

Some examples are provided in the main folder. Here we show a concrete example on how to use our software to calculate the density of states of graphene. 
