/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

#include "systemInfo.hpp"
#include "dos.hpp"
#include "ldos.hpp"
#include "conductivity_dc.hpp"
#include "conductivity_optical.hpp"
#include "functions.hpp"
#include "cond_2order/conductivity_2order.hpp"


void choose_simulation_type(char, shell_input &);


template <typename U, unsigned DIM>
void calculate_conductivity_nonlinear(system_info<U, DIM>& sysinfo, shell_input & variables){
  conductivity_nonlinear<U, DIM> info(sysinfo, variables);
  if(info.isRequired and variables.CondOpt2_is_required){
    verbose_message(
        "Retrieving nonlinear conductivity...\n "
        );
    variables.printOpt2();
    info.set_default_parameters();
    info.fetch_parameters();
    info.override_parameters();
    info.calculate();
  }
};

template <typename U, unsigned DIM>
void calculate_conductivity_optical(system_info<U, DIM>& sysinfo, shell_input & variables){
  conductivity_optical<U, DIM> info(sysinfo, variables);
  if(info.isRequired and variables.CondOpt_is_required){
    verbose_message("Retrieving optical conductivity...\n");
    variables.printOpt();
    info.fetch_parameters();
    info.override_parameters();
    info.calculate_efficient();
  }
};

template <typename U, unsigned DIM>
void calculate_dos(system_info<U, DIM>& sysinfo, shell_input & variables){
  dos<U, DIM> info(sysinfo, variables); // The constructor checks whether this quantity is required
  if(info.isRequired and variables.DOS_is_required){
    info.calculate();
  }
};

template <typename U, unsigned DIM>
void calculate_ldos(system_info<U, DIM>& sysinfo, shell_input & variables){
  ldos<U, DIM> info(sysinfo, variables); // The constructor checks whether this quantity is required
  if(info.isRequired and variables.lDOS_is_required){
    info.calculate();
  }
};

template <typename U, unsigned DIM>
void calculate_simple(char *name, shell_input & variables){
	debug_message("Entered calculate_simple.\n");
  // Attempts to calculate all quantities implemented.
  // For each of those functions, KITE will attempt to retrieve the relevant quantities
  // from the configuration file. If they cannot be retrieved, that means that
  // they were not requested, and the program moves on to the next quantity.
  
  system_info<U, DIM> info;
  info = system_info<U, DIM>(std::string(name));
  info.read();

  verbose_message("----------------- CALCULATIONS ----------------- \n");
  calculate_dos<U, DIM>(info, variables);
  calculate_ldos<U, DIM>(info, variables);

  conductivity_dc<U, DIM> condDC(info, variables);
  if(condDC.isRequired and variables.CondDC_is_required) condDC.calculate();

  calculate_conductivity_optical<U, DIM>(info, variables);
  calculate_conductivity_nonlinear<U, DIM>(info, variables);
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
            calculate_simple<float, 1u>(name, variables);
            break;
        }
        case 1:
        {
            calculate_simple<float, 2u>(name, variables);
            break;
        }
        case 2:
        {
            calculate_simple<float, 3u>(name, variables);
            break;			
        }
        case 3:
        {
            calculate_simple<double, 1u>(name, variables);
            break;
        }
#endif
        case 4:
        {
            calculate_simple<double, 2u>(name, variables);
            break;
        }
        case 5:
        {
            calculate_simple<double, 3u>(name, variables);
            break;
        }
#ifdef debug1
        case 6:
        {
            calculate_simple<long double, 1u>(name, variables);
            break;
        }
        case 7:
        {
            calculate_simple<long double, 2u>(name, variables);
            break;
        }
        case 8:
        {
            calculate_simple<long double, 3u>(name, variables);
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
