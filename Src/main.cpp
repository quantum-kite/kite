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

#define MEMORY   10
#define PATTERNS  4
#define NGHOSTS   2
#ifndef STRIDE
#define STRIDE    128
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



int main(int argc, char *argv[])
{  
  /* Define General characteristics of the data */ 
  int precision = 1, dim, is_complex;

  H5::H5File *file = new H5::H5File(argv[1], H5F_ACC_RDONLY);
  get_hdf5(&is_complex, file, (char *) "/IS_COMPLEX");
  get_hdf5(&precision,  file, (char *) "/PRECISION");
  get_hdf5(&dim,        file, (char *) "/DIM");
  
  file->close();
  if(dim < 1 || dim > 3)
    exit(0);
  
  if(precision < 0 || precision > 2)
    exit(0);
  
  if(is_complex < 0 || is_complex > 1)
    exit(0);

  int index =   dim - 1 + 3 * precision + is_complex * 3 * 3; 

  switch (index ) {
    /*
     * Float Real 
     */
  case 0:
    {
      class GlobalSimulation <float, 1u> h(argv[1]);
      break;
    }
  case 1:
    {
      class GlobalSimulation <float, 2u> h(argv[1]);
      break;
    }
  case 2:
    {
      class GlobalSimulation <float, 3u> h(argv[1]);
      break;
    }
  case 3:
    {
      class GlobalSimulation <double, 1u> h(argv[1]);
      break;
    }
  case 4:
    {
      class GlobalSimulation <double, 2u> h(argv[1]);
      break;
    }
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
  case 13:
    {
      class GlobalSimulation <std::complex<double>, 2u> h(argv[1]);
      break;
    }
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

  default:
    exit(0);
  } 
  
  return 1;
}

