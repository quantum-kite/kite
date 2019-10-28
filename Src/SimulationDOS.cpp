#include "Generic.hpp"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"
#include "Global.hpp"
#include "Random.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
template <typename T, unsigned D>
class Hamiltonian;
template <typename T, unsigned D>
class KPM_Vector;
#include "queue.hpp"
#include "Simulation.hpp"
#include "Hamiltonian.hpp"
#include "KPM_VectorBasis.hpp"
#include "KPM_Vector.hpp"

template <typename T,unsigned D>
void Simulation<T,D>::store_MU(Eigen::Array<T, -1, -1> *gamma){
    debug_message("Entered store_mu\n");

    long int nMoments   = gamma->rows();
    long int nPositions = gamma->cols();

#pragma omp master
	Global.general_gamma = Eigen::Array<T, -1, -1 > :: Zero(nMoments, nPositions);
#pragma omp barrier
#pragma omp critical
	Global.general_gamma += *gamma;
#pragma omp barrier
    
    
#pragma omp master
{
    H5::H5File * file = new H5::H5File(name, H5F_ACC_RDWR);
    write_hdf5(Global.general_gamma, file, "/Calculation/dos/MU");
    file->close();
    delete file;
}
#pragma omp barrier    
    debug_message("Left store_mu\n");
}


template <typename T,unsigned D>
void Simulation<T,D>::calc_DOS(){
    debug_message("Entered Simulation::calc_DOS\n");

    // Make sure that all the threads are ready before opening any files
    // Some threads could still be inside the Simulation constructor
    // This barrier is essential
#pragma omp barrier

  int NMoments, NRandom, NDisorder;
  bool local_calculate_dos = false;
#pragma omp master
{
  H5::H5File * file = new H5::H5File(name, H5F_ACC_RDONLY);
  Global.calculate_dos = false;
  try{
    int dummy_variable;
    get_hdf5<int>(&dummy_variable,  file, (char *)   "/Calculation/dos/NumMoments");
    Global.calculate_dos = true;
  } catch(H5::Exception& e) {debug_message("DOS: no need to calculate DOS.\n");}
  file->close();
  delete file;
}
#pragma omp barrier
#pragma omp critical
  local_calculate_dos = Global.calculate_dos;

#pragma omp barrier

if(local_calculate_dos){
#pragma omp master
      {
        std::cout << "Calculating DOS.\n";
      }
#pragma omp barrier
#pragma omp critical
{
    H5::H5File * file = new H5::H5File(name, H5F_ACC_RDONLY);

    debug_message("DOS: checking if we need to calculate DOS.\n");
    get_hdf5<int>(&NMoments,  file, (char *)   "/Calculation/dos/NumMoments");
    get_hdf5<int>(&NDisorder, file, (char *)   "/Calculation/dos/NumDisorder");
    get_hdf5<int>(&NRandom,   file, (char *)   "/Calculation/dos/NumRandoms");
    file->close();
    delete file;

    if(NDisorder <= 0){
      std::cout << "Cannot calculate Density of states with nonpositive NDisorder\n";
      exit(0);
    }
    if(NMoments <= 0){
      std::cout << "Cannot calculate Density of states with nonpositive NMoments\n";
      exit(0);
    }
    if(NRandom <= 0){
      std::cout << "Cannot calculate Density of states with nonpositive NRandom\n";
      exit(0);
    }

}
#pragma omp barrier
  DOS(NMoments, NRandom, NDisorder);
  }

}
template <typename T,unsigned D>

void Simulation<T,D>::DOS(int NMoments, int NRandom, int NDisorder){
  debug_message("Entered Simulation::DOS\n");
  std::vector<std::vector<unsigned>> indices = process_string("");
  Gamma1D(NRandom, NDisorder, NMoments, indices, "/Calculation/dos/MU");
  debug_message("Left Simulation::DOS\n");
}



template class Simulation<float ,1u>;
template class Simulation<double ,1u>;
template class Simulation<long double ,1u>;
template class Simulation<std::complex<float> ,1u>;
template class Simulation<std::complex<double> ,1u>;
template class Simulation<std::complex<long double> ,1u>;

template class Simulation<float ,3u>;
template class Simulation<double ,3u>;
template class Simulation<long double ,3u>;
template class Simulation<std::complex<float> ,3u>;
template class Simulation<std::complex<double> ,3u>;
template class Simulation<std::complex<long double> ,3u>;

template class Simulation<float ,2u>;
template class Simulation<double ,2u>;
template class Simulation<long double ,2u>;
template class Simulation<std::complex<float> ,2u>;
template class Simulation<std::complex<double> ,2u>;
template class Simulation<std::complex<long double> ,2u>;
