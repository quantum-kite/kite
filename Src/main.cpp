/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

#include <stdio.h>
#include <omp.h>
#include <time.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "H5Cpp.h"
#include <H5Group.h>
#include <complex>
#include <random>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <typeinfo>
#include <type_traits>
#include <complex>
#include <iostream>
#include <chrono>
#include <thread>
#include <cmath>
#include <math.h>
#include <initializer_list>

// Set of compilation parameters chosen in the Makefile
// MEMORY is the number of KPM vectors stored in the memory while calculating Gamma2D
// STRIDE is the size of the memory blocks used in the program
// COMPILE_MAIN is a flag to prevent compilation of unnecessary parts of the code when testing
#ifndef MEMORY
#define MEMORY 128
#endif

#ifndef STRIDE
#define STRIDE 64
#endif

#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef VERBOSE
#define VERBOSE 1
#endif

#ifndef COMPILE_MAIN
#define COMPILE_MAIN 1
#endif

#ifndef ESTIMATE_TIME
#define ESTIMATE_TIME 1
#endif

// other compilation parameters not set in the Makefile
// NGHOSTS is the extra length in each direction, to be used with the blocks of size STRIDE
#define PATTERNS  4
#define NGHOSTS   2
#define VVERBOSE 0
#define NUM_GHOST_CORR 0
#define SSPRINT 8

// These are the verbose and debug messages
#define outcol "\033[1;31m"
#define outres "\033[0m"

#ifdef VERBOSE
#if VERBOSE==1
#define verbose_message(VAR)              \
  _Pragma("omp master")                   \
  {                                       \
    std::cout<<VAR<<std::flush;           \
  }                                       \
  _Pragma("omp barrier")
#else
#define verbose_message(VAR) 
#endif
#else
#define verbose_message(VAR) 
#endif


#ifdef VVERBOSE
#if VVERBOSE==1
#define vverbose_message(VAR)              \
  _Pragma("omp master")                   \
  {                                       \
    std::cout<<VAR<<std::flush;           \
  }                                       \
  _Pragma("omp barrier")
#else
#define vverbose_message(VAR) 
#endif
#else
#define vverbose_message(VAR) 
#endif

#ifdef DEBUG
#if DEBUG==1
#define debug_message(VAR) std::cout<<VAR<<std::flush; 
#else
#define debug_message(VAR) 
#endif
#else
#define debug_message(VAR) 
#endif


template<typename T, unsigned D>
class Simulation;
#include "Global.hpp"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"
#include "Random.hpp"
#include "LatticeStructure.hpp"
#include "Hamiltonian.hpp"
#include "KPM_Vector.hpp"
#include "KPM_Vector2D.hpp"
#include "Simulation.hpp"

typedef int indextype;



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
      "To estimate the calculation time, set ESTIMATE_TIME to 1 in the Makefile. \n"
      "------------------------------------------------------------------------- \n\n"
      );

  verbose_message("------- FLAGS SET -------\n");
  verbose_message("Flags set at compilation:\n");
  verbose_message("DEBUG: "); verbose_message(DEBUG); verbose_message("\n");
  verbose_message("VERBOSE: "); verbose_message(VERBOSE); verbose_message("\n");
  verbose_message("ESTIMATE_TIME: "); verbose_message(ESTIMATE_TIME); verbose_message("\n");
  verbose_message("-------------------------\n");

  verbose_message("\nStarting program...\n\n");
  debug_message("Starting program. The messages in red are debug messages. They may be turned off by setting DEBUG 0 in main.cpp\n");

  /* Define General characteristics of the data */  
  int precision = 1, dim, is_complex;

  H5::H5File *file = new H5::H5File(argv[1], H5F_ACC_RDONLY);
  get_hdf5(&is_complex, file, (char *) "/IS_COMPLEX");
  get_hdf5(&precision,  file, (char *) "/PRECISION");
  get_hdf5(&dim,        file, (char *) "/DIM");
  
  
  file->close();
  
  // Verify if the values passed to the program are valid. If they aren't
  // the program should notify the user and exit with error 1.
  if(dim < 1 || dim > 3){
    std::cout << "Invalid number of dimensions. The code is only valid for 2D. Exiting.\n";
    exit(1);
  }
  if(precision < 0 || precision > 2){
    std::cout << "Please use a valid value for the numerical precision. Accepted values: 0, 1, 2. Exiting.\n";
    exit(1);
  }
  if(is_complex < 0 || is_complex > 1){
    std::cout << "Bad complex flag. It has to be either 0 or 1. Exiting.\n";
    exit(1);
  }


  // Decide which version of the program should run. This depends on the
  // precision, the dimension and whether or not we want complex functions.
  // When COMPILE_MAIN is 1, only a small portion of the code is compiled.
  // This speeds up the compilation time considerably, while only support
  int index =   dim - 1 + 3 * precision + is_complex * 3 * 3; 
  switch (index ) {
#if COMPILE_MAIN==1
  case 0:
    {
      class GlobalSimulation <float, 1u> h(argv[1]); // float real 1D
      break;
    }
#endif
  case 1:
    {
      class GlobalSimulation <float, 2u> h(argv[1]); // float real 2D
      break;
    }
#if COMPILE_MAIN==1
  case 2:
    {
      class GlobalSimulation <float, 3u> h(argv[1]); // float real 3D
      break;
    }
  case 3:
      {
      class GlobalSimulation <double, 1u> h(argv[1]); // double real 1D
      break;
      }
#endif
  case 4:
      {
      class GlobalSimulation <double, 2u> h(argv[1]); //double real 2D. You get the picture.
      break;
      }
#if COMPILE_MAIN==1  
  case 5:
      {
      class GlobalSimulation <double, 3u> h(argv[1]);
      break;
      }
  case 6:
      {
      class GlobalSimulation <long double, 1u> h(argv[1]);
      break;
      }
  case 7:
      {
      class GlobalSimulation <long double, 2u> h(argv[1]);
      break;
      }
  case 8:
      {
      class GlobalSimulation <long double, 3u> h(argv[1]);
      break;
      }
  case 9:
      {
      class GlobalSimulation <std::complex<float>, 1u> h(argv[1]);
      break;
      }
  case 10:
      {
      class GlobalSimulation <std::complex<float>, 2u> h(argv[1]);
      break;
      }
  case 11:
      {
      class GlobalSimulation <std::complex<float>, 3u> h(argv[1]);
      break;
      }
  case 12:
      {
      class GlobalSimulation <std::complex<double>, 1u> h(argv[1]);
      break;
      }
#endif    
  case 13:
      {
      class GlobalSimulation <std::complex<double>, 2u> h(argv[1]);
      break;
      }
#if COMPILE_MAIN==1  
  case 14:
      {
      class GlobalSimulation <std::complex<double>, 3u> h(argv[1]);
      break;
      }
  case 15:
      {
      class GlobalSimulation <std::complex<long double>, 1u> h(argv[1]);
      break;
      }
  case 16:
      {
      class GlobalSimulation <std::complex<long double>, 2u> h(argv[1]);
      break;
      }
  case 17:
      {
      class GlobalSimulation <std::complex<long double>, 3u> h(argv[1]);
      break;
      }
#endif
  default:
      { 
      std::cout << "Unexpected parameters. Please use valid values for the precision, dimension and 'complex' flag.";
      std::cout << "Check if the code has been compiled with support for complex functions. Exiting.\n";
      exit(0);
      }
  }
  
  debug_message("Program ended with success!\n");
  verbose_message("Done.\n");
  return 0;
}

