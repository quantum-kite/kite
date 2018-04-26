#include "headers.hpp"
#include "ComplexTraits.hpp"
#include "H5Cpp.h"
#include "tensor.hpp"
#include "myHDF5.hpp"
#include "info.hpp"
#include "calculate_simple.hpp"
#include "parse_input.hpp"

//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################
//#####################################################################################

//TO DO:
// Clarify what each quantity to be calculate_simpled means. Cond, SingleCond, OptCond????

//https://support.hdfgroup.org/HDF5/doc/cpplus_RM/readdata_8cpp-example.html

#define debug 1

void choose_simulation_type(char *name){
	debug_message("Entered choose_simulation.\n");
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
	
	file.close();
	
	 if(dim < 1 || dim > 3)
		exit(0);
  
	if(precision < 0 || precision > 2)
		exit(0);

	if(complex < 0 || complex > 1)
		exit(0);

	int index = dim - 1 + 3 * precision;
	switch (index ) {
#ifdef debug1
		case 0:
		{
			verbose_message("The program is using data type: float.\nDimension: 1.\n");
			calculate_simple<float, 1u>(name);
			break;
		}
		case 1:
		{
			verbose_message("The program is using data type: float.\nDimension: 2.\n");
			calculate_simple<float, 2u>(name);
			break;
		}
		case 2:
		{
			verbose_message("The program is using data type: float.\nDimension: 3.\n");
			calculate_simple<float, 3u>(name);
			break;			
		}
		case 3:
		{
			verbose_message("The program is using data type: double.\nDimension: 1.\n");
			calculate_simple<double, 1u>(name);
			break;
		}
#endif
		case 4:
		{
			verbose_message("The program is using data type: double.\nDimension: 2.\n");
			calculate_simple<double, 2u>(name);
			break;
		}
		case 5:
		{
			verbose_message("The program is using data type: double.\nDimension: 3.\n");
			calculate_simple<double, 3u>(name);
			break;
		}
#ifdef debug1
		case 6:
		{
			verbose_message("The program is using data type: long double.\nDimension: 1.\n");
			calculate_simple<long double, 1u>(name);
			break;
		}
		case 7:
		{
			verbose_message("The program is using data type: long double.\nDimension: 2.\n");
			calculate_simple<long double, 2u>(name);
			break;
		}
		case 8:
		{
			verbose_message("The program is using data type: long double.\nDimension: 3.\n");
			calculate_simple<long double, 3u>(name);
			break;
		}
#endif
		default:
		{
			std::cout << "Please enter a valid data type: float, double or long double and a valid dimension: 2. Exiting.\n" << std::flush;
			exit(1);
			break;
		}
	} 
	debug_message("Left choose_simulation.\n");
}


int main(int argc, char *argv[]){
	//parser(argc, argv);
	choose_simulation_type(argv[1]);
	return 0;
}
