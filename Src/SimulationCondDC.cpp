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


    std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n";
    //std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n";
  //int NMoments, NRandom, NDisorder, direction;
    //std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n" << std::flush;
  //bool local_calculate_conddc = false;
    //std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n" << std::flush;
//#pragma omp master
//{
    //std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n";
  //H5::H5File * file = new H5::H5File(name, H5F_ACC_RDONLY);
    //std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n";
  //Global.calculate_conddc = false;
  //try{
    //int dummy_variable;
    //std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n" << std::flush;
    //get_hdf5<int>(&dummy_variable,  file, (char *)   "/Calculation/conductivity_dc/NumMoments");
    //std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n" << std::flush;
    //Global.calculate_conddc = true;
    //std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n" << std::flush;
  //} catch(H5::Exception& e) {debug_message("CondDC: no need to calculate CondDC.\n");}
  //file->close();
    //std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n" << std::flush;
  //delete file;
    //std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n" << std::flush;
//}
//#pragma omp barrier
//#pragma omp critical
  //local_calculate_conddc = Global.calculate_conddc;
    //std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n" << std::flush;

//#pragma omp barrier

    //std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n" << std::flush;
//if(local_calculate_conddc){
    //std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n" << std::flush;
//#pragma omp critical
//{
    //std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n";
    //H5::H5File * file = new H5::H5File(name, H5F_ACC_RDONLY);

    //debug_message("DC conductivity: checking if we need to calculate DC conductivity.\n");
    //get_hdf5<int>(&direction, file, (char *) "/Calculation/conductivity_dc/Direction");
    //std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n";
    //get_hdf5<int>(&NMoments, file, (char *)  "/Calculation/conductivity_dc/NumMoments");
    //std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n";
    //get_hdf5<int>(&NRandom, file, (char *)   "/Calculation/conductivity_dc/NumRandoms");
    //std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n";
    //get_hdf5<int>(&NDisorder, file, (char *) "/Calculation/conductivity_dc/NumDisorder");
    //std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n";

    //std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n";
    //file->close();
    //delete file;

//}
    //std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n";
  //CondDC(NMoments, NRandom, NDisorder, direction);
    //std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n";
  //}
    std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n";

}
template <typename T,unsigned D>

void Simulation<T,D>::CondDC(int NMoments, int NRandom, int NDisorder, int direction){
  std::string dir(num2str2(direction));
    std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n";
  std::string dirc = dir.substr(0,1)+","+dir.substr(1,2);
    std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n";
  Gamma2D(NRandom, NDisorder, {NMoments,NMoments}, process_string(dirc), "/Calculation/conductivity_dc/Gamma"+dir);
    std::cout << "line: " << __LINE__ << " in file " << __FILE__ << "\n";
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
