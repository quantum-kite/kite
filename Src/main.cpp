#include "Generic.hpp"

template<typename T, unsigned D>
class Simulation;
#include "Global.hpp"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"
#include "Random.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
#include "Hamiltonian.hpp"
#include "KPM_VectorBasis.hpp"
#include "KPM_Vector.hpp"
#include "queue.hpp"
#include "Simulation.hpp"
#include "SimulationGlobal.hpp"
#include "messages.hpp"

typedef int indextype;



int main(int argc, char *argv[]){  
  //(void) argc;
    
  print_header_message();
  print_info_message();
  print_flags_message();



  verbose_message("\nStarting program...\n\n");
  debug_message("Starting program. The messages in red are debug messages. They may be turned off by setting DEBUG 0 in main.cpp\n");
  if(argc < 2){
    std::cout << "No configuration file found. Exiting.\n";
    exit(1);
  }

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
    std::cout << "Invalid number of dimensions. The code is only valid for 2D or 3D. Exiting.\n";
    exit(1);
  }
  if(precision < 0 || precision > 2){
    std::cout << "Please use a valid value for the numerical precision. Accepted values: 0, 1, 2. Exiting.\n";
    exit(1);
  }
  if(is_complex != 0 && is_complex != 1){
    std::cout << "Bad complex flag. It has to be either 0 or 1. Exiting.\n";
    exit(1);
  }


  // Decide which version of the program should run. This depends on the
  // precision, the dimension and whether or not we want complex functions.
  int index =   dim - 1 + 3 * precision + is_complex * 3 * 3; 
  switch (index ) {
  case 0:
    {
      class GlobalSimulation <float, 1u> h(argv[1]); // float real 1D
      break;
    }
  case 1:
    {
      class GlobalSimulation <float, 2u> h(argv[1]); // float real 2D
      break;
    }
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
  case 4:
      {
      class GlobalSimulation <double, 2u> h(argv[1]); //double real 2D. You get the picture.
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
      { 
      std::cout << "Unexpected parameters. Please use valid values for the precision, dimension and 'complex' flag.";
      std::cout << "Check if the code has been compiled with support for complex functions. Exiting.\n";
      exit(1);
      }
  }
  
  verbose_message("Done.\n");
  return 0;
}

