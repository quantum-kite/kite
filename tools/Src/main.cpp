/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

#include "headers.hpp"
#include "ComplexTraits.hpp"
#include "H5Cpp.h"
#include "tensor.hpp"
#include "myHDF5.hpp"
#include "systemInfo.hpp"
#include "conductivity_dc.hpp"
#include "conductivity_optical.hpp"
#include "conductivity_nonlinear.hpp"
#include "dos.hpp"
#include "calculate_simple.hpp"
#include "parse_input.hpp"

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
	
	if(dim < 1 || dim > 2){
    std::cout << "Invalid value for the dimension of the system. Exiting.\n";
		exit(0);
  }
  
	if(precision < 0 || precision > 2){
    std::cout << "Invalid value for the precision of the data values. Exiting. \n";
		exit(0);
  }

	if(complex < 0 || complex > 1)
		exit(0);

	int index = dim - 1 + 3 * precision;
	switch (index ) {
#ifdef debug1
		case 0:
		{
			debug_message("The program is using data type: float.\nDimension: 1.\n");
			calculate_simple<float, 1u>(name);
			break;
		}
		case 1:
		{
			debug_message("The program is using data type: float.\nDimension: 2.\n");
			calculate_simple<float, 2u>(name);
			break;
		}
		case 2:
		{
			debug_message("The program is using data type: float.\nDimension: 3.\n");
			calculate_simple<float, 3u>(name);
			break;			
		}
		case 3:
		{
			debug_message("The program is using data type: double.\nDimension: 1.\n");
			calculate_simple<double, 1u>(name);
			break;
		}
#endif
		case 4:
		{
			debug_message("The program is using data type: double.\nDimension: 2.\n");
			calculate_simple<double, 2u>(name);
			break;
		}
		case 5:
		{
			debug_message("The program is using data type: double.\nDimension: 3.\n");
			calculate_simple<double, 3u>(name);
			break;
		}
#ifdef debug1
		case 6:
		{
			debug_message("The program is using data type: long double.\nDimension: 1.\n");
			calculate_simple<long double, 1u>(name);
			break;
		}
		case 7:
		{
			debug_message("The program is using data type: long double.\nDimension: 2.\n");
			calculate_simple<long double, 2u>(name);
			break;
		}
		case 8:
		{
			debug_message("The program is using data type: long double.\nDimension: 3.\n");
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
  verbose_message(
      "\n+------------------------------------------------------------------------+\n"
      "|            Chebyshev Polynomial Green's Function Approach              | \n"
      "|             to Real-Space Quantum Transport Simulations                | \n"    
      "|                                                                        | \n"                                              
      "|                     KITE | Pre-Release Beta 0.1                        | \n"         
      "|                     Kite home: quantum-kite.com                        | \n"
      "|                                                                        | \n"
      "|    Created by Simao M. Joao, Joao V. Lopes (Universidade do Porto),    | \n"
      "|      Tatiana G. Rappoport (Universidade Federal Rio de Janeiro),       | \n"
      "|        Misa Andelkovic, Lucian Covaci (University of Antwerp)          | \n"
      "|                and Aires Ferreira (University of York)                 | \n"
      "|                                                                        | \n"                                            
      "|            Funded by The Royal Society| royalsociety.org               | \n"
      "|                                                                        | \n"
      "|  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,            | \n"
      "|                      S. M. Joao, J. V. Lopes, T. G. Rappoport          | \n"
      "|                                                                        | \n"
      "|  This program is free software: you can redistribute it and/or modify  | \n"
      "|  it under the terms of the GNU General Public License as published by  | \n"
      "|  the Free Software Foundation, either version 3 of the License, or     | \n"
      "|  (at your option) any later version.                                   | \n"
      "|                                                                        | \n"
      "|  This program is distributed in the hope that it will be useful,       | \n"
      "|  but WITHOUT ANY WARRANTY; without even the implied warranty of        | \n"
      "|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                  | \n"
      "|  See the GNU General Public License for more details.                  | \n"
      "+------------------------------------------------------------------------+\n"
      );

  verbose_message(
      "\n------------------------------ INFORMATION ------------------------------\n"
      "Linear response functions in units of e^2/h.                              \n"
      "To stop these messages, set VERBOSE to 0 in the Makefile.                 \n"
      "To see debug messages, set DEBUG to 1 in the Makefile.                    \n"
      "------------------------------------------------------------------------- \n\n"
      );

  verbose_message("------- FLAGS SET -------\n");
  verbose_message("DEBUG: "); verbose_message(DEBUG); verbose_message("\n");
  verbose_message("VERBOSE: "); verbose_message(VERBOSE); verbose_message("\n");
  verbose_message("-------------------------\n\n");


  verbose_message("\nStarting program...\n\n");


	//parser(argc, argv);
	choose_simulation_type(argv[1]);
	return 0;
  verbose_message("Done.\n");
}
