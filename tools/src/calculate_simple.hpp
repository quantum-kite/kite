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
void calculate_conductivity_dc(system_info<U, DIM>& sysinfo, shell_input & variables){
  conductivity_dc<U, DIM> info(sysinfo, variables);
  if(info.isRequired and variables.CondDC_is_required){
    //variables.printDC();
    info.calculate();
  }
};

template <typename U, unsigned DIM>
void calculate_dos(system_info<U, DIM>& sysinfo, shell_input & variables){
  dos<U, DIM> info(sysinfo, variables); // The constructor checks whether this quantity is required
  if(info.isRequired and variables.DOS_is_required){
    verbose_message("Retrieving DOS...\n");
    variables.printDOS();
    info.fetch_parameters();
    info.override_parameters();
    info.calculate();
  }
};

template <typename U, unsigned DIM>
void calculate_ldos(system_info<U, DIM>& sysinfo, shell_input & variables){
  ldos<U, DIM> info(sysinfo, variables); // The constructor checks whether this quantity is required
  if(info.isRequired and variables.lDOS_is_required){
    verbose_message("Retrieving lDOS...\n");
    variables.printlDOS();
    info.fetch_parameters();
    info.override_parameters();
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
  calculate_conductivity_dc<U, DIM>(info, variables);
  calculate_conductivity_optical<U, DIM>(info, variables);
  calculate_conductivity_nonlinear<U, DIM>(info, variables);
  verbose_message("------------------------------------------------ \n\n");


	debug_message("Left calculate_simple.\n");
}



