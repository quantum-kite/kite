#include "headers.hpp"
#include "ComplexTraits.hpp"
#include "H5Cpp.h"
#include "tensor.hpp"
#include "myHDF5.hpp"
#include "info.hpp"
#include "calculate.hpp"

#define debug 1

//void calc_dos(double, int);
//void calc_optical_cond(int, int, int, double, double);
//void single_shot(Eigen::Array<double, -1, 1> energies);
//void single_shot_matrix(Eigen::Array<double, -1, 1> energies);

//#################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

//TO DO:
// Clarify what each quantity to be calculated means. Cond, SingleCond, OptCond????

//https://support.hdfgroup.org/HDF5/doc/cpplus_RM/readdata_8cpp-example.html



void choose_simulation_type(char *name){
	/* The type of data used for the simulation is not known beforehand. It may be 
	 * float, double or long double. Each of those may or may not be complex. To simplify
	 * we assume everything is complex. If they are not, they are cast into their complex form.
	 * Furthermore, the dimension is also not known. Therefore, each and every one of those
	 * combinations of cases is compiled individually, so the correct one will be
	 * available at runtime. This function reads from the hdf5 file which are the
	 * parameters and runs the correct version of the code to deal with it.
	 */
	
	H5::H5File file;
	file = H5::H5File(name, H5F_ACC_RDONLY);
	int precision = 1, dim, complex;

	get_hdf5(&complex, &file, (char *) "/IS_COMPLEX");
	get_hdf5(&precision,  &file, (char *) "/PRECISION");
	get_hdf5(&dim,        &file, (char *) "/DIM");
	
	 if(dim < 1 || dim > 3)
		exit(0);
  
	if(precision < 0 || precision > 2)
		exit(0);

	if(complex < 0 || complex > 1)
		exit(0);

	int index = dim - 1 + 3 * precision;
	switch (index ) {
		case 0:
		{
			if(debug)std::cout << "The program is using data type: float.\nDimension: 1.\n" << std::flush;
			calculate<float, 1u>(name);
			break;
		}
		case 1:
		{
			if(debug)std::cout << "The program is using data type: float.\nDimension: 2.\n" << std::flush;
			calculate<float, 2u>(name);
			break;
		}
		case 2:
		{
			if(debug)std::cout << "The program is using data type: float.\nDimension: 3.\n" << std::flush;
			calculate<float, 3u>(name);
			break;			
		}
		case 3:
		{
			if(debug)std::cout << "The program is using data type: double.\nDimension: 1.\n" << std::flush;
			calculate<double, 1u>(name);
			break;
		}
		case 4:
		{
			if(debug)std::cout << "The program is using data type: double.\nDimension: 2.\n" << std::flush;
			calculate<double, 2u>(name);
			break;
		}
		case 5:
		{
			if(debug)std::cout << "The program is using data type: double.\nDimension: 3.\n" << std::flush;
			calculate<double, 3u>(name);
			break;
		}
		case 6:
		{
			if(debug)std::cout << "The program is using data type: long double.\nDimension: 1.\n" << std::flush;
			calculate<long double, 1u>(name);
			break;
		}
		case 7:
		{
			if(debug)std::cout << "The program is using data type: long double.\nDimension: 2.\n" << std::flush;
			calculate<long double, 2u>(name);
			break;
		}
		case 8:
		{
			if(debug)std::cout << "The program is using data type: long double.\nDimension: 3.\n" << std::flush;
			calculate<long double, 3u>(name);
			break;
		}
		default:
		{
			if(debug)std::cout << "Please enter a valid data type: float, double or long double and a valid dimension: 1, 2, \n" << std::flush;
			exit(0);
			break;
		}
	} 
}


int main(int argc, char *argv[]){
	choose_simulation_type(argv[1]);
	return 0;
}
