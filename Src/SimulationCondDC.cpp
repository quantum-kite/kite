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
void Simulation<T,D>::calc_conddc(){
    debug_message("Entered Simulation::calc_conddc\n");

    // Make sure that all the threads are ready before opening any files
    // Some threads could still be inside the Simulation constructor
    // This barrier is essential
#pragma omp barrier

  int NMoments, NRandom, NDisorder, direction;
  bool local_calculate_conddc = false;
#pragma omp master
{
  H5::H5File * file = new H5::H5File(name, H5F_ACC_RDONLY);
  Global.calculate_conddc = false;
  try{
    int dummy_variable;
    get_hdf5<int>(&dummy_variable,  file, (char *)   "/Calculation/conductivity_dc/NumMoments");
    Global.calculate_conddc = true;
  } catch(H5::Exception& e) {debug_message("CondDC: no need to calculate CondDC.\n");}
  file->close();
  delete file;
}
#pragma omp barrier
#pragma omp critical
  local_calculate_conddc = Global.calculate_conddc;

#pragma omp barrier

if(local_calculate_conddc){
#pragma omp master
      {
        std::cout << "Calculating CondDC.\n";
      }
#pragma omp barrier
#pragma omp critical
{
    H5::H5File * file = new H5::H5File(name, H5F_ACC_RDONLY);

    debug_message("DC conductivity: checking if we need to calculate DC conductivity.\n");
    get_hdf5<int>(&direction, file, (char *) "/Calculation/conductivity_dc/Direction");
    get_hdf5<int>(&NMoments, file, (char *)  "/Calculation/conductivity_dc/NumMoments");
    get_hdf5<int>(&NRandom, file, (char *)   "/Calculation/conductivity_dc/NumRandoms");
    get_hdf5<int>(&NDisorder, file, (char *) "/Calculation/conductivity_dc/NumDisorder");

    file->close();
    delete file;

}
  CondDC(NMoments, NRandom, NDisorder, direction);
  }

}
template <typename T,unsigned D>

void Simulation<T,D>::CondDC(int NMoments, int NRandom, int NDisorder, int direction){
  std::string dir(num2str2(direction));
  std::string dirc = dir.substr(0,1)+","+dir.substr(1,2);
  Gamma2D(NRandom, NDisorder, {NMoments,NMoments}, process_string(dirc), "/Calculation/conductivity_dc/Gamma"+dir);
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
