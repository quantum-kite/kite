/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/


#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <string>
#include <omp.h>

#include "H5Cpp.h"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"

#include "parse_input.hpp"
#include "systemInfo.hpp"
#include "dos.hpp"
#include "ldos.hpp"
#include "arpes.hpp"
#include "cond_dc/conductivity_dc.hpp"
#include "conductivity_optical.hpp"
//#include "functions.hpp"
#include "cond_2order/conductivity_2order.hpp"
#include "calculate.hpp"

#include "macros.hpp"

template <typename U, unsigned DIM>
void calculate(char *name, shell_input & variables){
	debug_message("Entered calculate_simple.\n");
  // Attempts to calculate all quantities implemented.
  // For each of those functions, KITE will attempt to retrieve the relevant quantities
  // from the configuration file. If they cannot be retrieved, that means that
  // they were not requested, and the program moves on to the next quantity.
  
  system_info<U, DIM> info;
  info = system_info<U, DIM>(std::string(name));
  info.read();
  info.print_info();

  verbose_message("----------------- CALCULATIONS ----------------- \n");
  arpes<U, DIM>                   arpes(info, variables); 
  ldos<U, DIM>                    lDOS(info, variables); 
  dos<U, DIM>                     DOS(info, variables); 
  conductivity_dc<U, DIM>         condDC(info, variables);
  conductivity_optical<U, DIM>    condOpt(info, variables);
  conductivity_nonlinear<U, DIM>  condOpt2(info, variables);
  verbose_message("------------------------------------------------ \n\n");


	debug_message("Left calculate_simple.\n");
}



void choose_simulation_type(char *name, shell_input & variables){
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
    debug_message("IS_COMPLEX "); debug_message(complex); debug_message("\n");
    debug_message("precision "); debug_message(precision); debug_message("\n");
    debug_message("dim "); debug_message(dim); debug_message("\n");

    file.close();
    
    if(precision < 0 || precision > 2){
    std::cout << "Invalid value for the precision of the data values. Exiting. \n";
        exit(0);
  }

    if(complex < 0 || complex > 1)
        exit(0);

    int index = dim - 1 + 3 * precision;
    switch (index ) {
        case 0:
        {
            debug_message("case ");debug_message(0); debug_message("\n");
            calculate<float, 1u>(name, variables);
            break;
        }
        case 1:
        {
            debug_message("case ");debug_message(1); debug_message("\n");
            calculate<float, 2u>(name, variables);
            break;
        }
        case 2:
        {
            debug_message("case ");debug_message(2); debug_message("\n");
            calculate<float, 3u>(name, variables);
            debug_message("left case ");debug_message(2); debug_message("\n");
            break;			
        }
        case 3:
        {
            debug_message("case ");debug_message(3); debug_message("\n");
            calculate<double, 1u>(name, variables);
            break;
        }
        case 4:
        {
            debug_message("case ");debug_message(4); debug_message("\n");
            calculate<double, 2u>(name, variables);
            break;
        }
        case 5:
        {
            debug_message("case ");debug_message(5); debug_message("\n");
            calculate<double, 3u>(name, variables);
            break;
        }
        case 6:
        {
            debug_message("case ");debug_message(6); debug_message("\n");
            calculate<long double, 1u>(name, variables);
            break;
        }
        case 7:
        {
            debug_message("case ");debug_message(7); debug_message("\n");
            calculate<long double, 2u>(name, variables);
            break;
        }
        case 8:
        {
            debug_message("case ");debug_message(8); debug_message("\n");
            calculate<long double, 3u>(name, variables);
            break;
        }
        default:
        {
            debug_message("case default"); debug_message("\n");
            std::cout << "Please enter a valid data type: float, double or long double and a valid dimension: 2. Exiting.\n" << std::flush;
            exit(1);
            break;
        }
    } 
    debug_message("Left choose_simulation.\n");
}
