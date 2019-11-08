#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
#include "Generic.hpp"
#include "Global.hpp"
#include "ComplexTraits.hpp"
#include "Random.hpp"
template <typename T, unsigned D>
class Hamiltonian;
template <typename T, unsigned D>
class KPM_Vector;
//#include "queue.hpp"
#include "Simulation.hpp"
#include "Hamiltonian.hpp"
#include "KPM_VectorBasis.hpp"
#include "KPM_Vector.hpp"


template<typename T, unsigned D>
KPM_VectorBasis<T,D>::KPM_VectorBasis(int mem,  Simulation<T,D> & sim) : memory(mem), simul(sim) {
  index  = 0;
  v = Eigen::Matrix <T, Eigen::Dynamic,  Eigen::Dynamic >::Zero(simul.r.Sized, memory);
}

template<typename T, unsigned D>
void KPM_VectorBasis<T,D>::set_index(int i) {
  index = i;
}

template<typename T, unsigned D>
void KPM_VectorBasis<T,D>::inc_index() {
  index = (index + 1) % memory;
}

template<typename T, unsigned D>
unsigned KPM_VectorBasis<T,D>::get_index() {
  return index;
}

template<typename T, unsigned D>
bool KPM_VectorBasis<T,D>::aux_test(T & x, T & y ) {
  return (abs(x - y) > std::numeric_limits<double>::epsilon());
}


// Instantiate KPM_VectorBasis
template class KPM_VectorBasis<float ,1u>;
template class KPM_VectorBasis<double ,1u>;
template class KPM_VectorBasis<long double ,1u>;
template class KPM_VectorBasis<std::complex<float> ,1u>;
template class KPM_VectorBasis<std::complex<double> ,1u>;
template class KPM_VectorBasis<std::complex<long double> ,1u>;

template class KPM_VectorBasis<float ,3u>;
template class KPM_VectorBasis<double ,3u>;
template class KPM_VectorBasis<long double ,3u>;
template class KPM_VectorBasis<std::complex<float> ,3u>;
template class KPM_VectorBasis<std::complex<double> ,3u>;
template class KPM_VectorBasis<std::complex<long double> ,3u>;

template class KPM_VectorBasis<float ,2u>;
template class KPM_VectorBasis<double ,2u>;
template class KPM_VectorBasis<long double ,2u>;
template class KPM_VectorBasis<std::complex<float> ,2u>;
template class KPM_VectorBasis<std::complex<double> ,2u>;
template class KPM_VectorBasis<std::complex<long double> ,2u>;
