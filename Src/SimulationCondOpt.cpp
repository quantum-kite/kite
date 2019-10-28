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
void Simulation<T,D>::calc_condopt(){
    debug_message("Entered Simulation::calc_condopt\n");

    // Make sure that all the threads are ready before opening any files
    // Some threads could still be inside the Simulation constructor
    // This barrier is essential
#pragma omp barrier

  int NMoments, NRandom, NDisorder, direction;
  bool local_calculate_condopt = false;
#pragma omp master
{
  H5::H5File * file = new H5::H5File(name, H5F_ACC_RDONLY);
  Global.calculate_condopt = false;
  try{
    int dummy_variable;
    get_hdf5<int>(&dummy_variable,  file, (char *)   "/Calculation/conductivity_optical/NumMoments");
    Global.calculate_condopt = true;
  } catch(H5::Exception& e) {debug_message("CondOpt: no need to calculate CondOpt.\n");}
  file->close();
  delete file;
}
#pragma omp barrier
#pragma omp critical
  local_calculate_condopt = Global.calculate_condopt;

#pragma omp barrier

if(local_calculate_condopt){
#pragma omp master
      {
        std::cout << "Calculating the optical conductivity.\n";
      }
#pragma omp barrier

#pragma omp critical
{
    H5::H5File * file = new H5::H5File(name, H5F_ACC_RDONLY);

    debug_message("Optical conductivity: checking if we need to calculate Condopt.\n");
    get_hdf5<int>(&direction, file, (char *) "/Calculation/conductivity_optical/Direction");
    get_hdf5<int>(&NMoments, file, (char *)  "/Calculation/conductivity_optical/NumMoments");
    get_hdf5<int>(&NRandom, file, (char *)   "/Calculation/conductivity_optical/NumRandoms");
    get_hdf5<int>(&NDisorder, file, (char *)   "/Calculation/conductivity_optical/NumDisorder");

    file->close();
    delete file;

}
  CondOpt(NMoments, NRandom, NDisorder, direction);
  }

}
template <typename T,unsigned D>

void Simulation<T,D>::CondOpt(int NMoments, int NRandom, int NDisorder, int direction){
  std::string dir(num2str2(direction));
  std::string dirc = dir.substr(0,1)+","+dir.substr(1,2);
  Gamma1D(NRandom, NDisorder, NMoments, process_string(dir), "/Calculation/conductivity_optical/Lambda"+dir);
  Gamma2D(NRandom, NDisorder, {NMoments,NMoments}, process_string(dirc), "/Calculation/conductivity_optical/Gamma"+dir);
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
