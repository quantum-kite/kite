template<typename T, unsigned D>
class Simulation;
#include "Generic.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
#include "ComplexTraits.hpp"
#include "Global.hpp"
#include "Random.hpp"
#include "Hamiltonian.hpp"
#include "KPM_VectorBasis.hpp"
#include "KPM_Vector.hpp"
//#include "queue.hpp"
#include "Simulation.hpp"

template <typename T, unsigned D>
KPM_Vector<T,D>::KPM_Vector(int mem, Simulation<T,D> & sim) : KPM_VectorBasis<T,D>(mem, sim), std(x.basis[1]), r(sim.r),h(sim.h),x(sim.r.Ld)  {
}

template <typename T, unsigned D>
KPM_Vector<T,D>::~KPM_Vector(void){}

template <typename T, unsigned D>
void KPM_Vector<T,D>::initiate_vector(){}

template <typename T, unsigned D>
T KPM_Vector<T,D>::get_point(){return v(0,0);}

template <typename T, unsigned D>
void KPM_Vector<T,D>::build_wave_packet(Eigen::Matrix<double,-1,-1> & k, Eigen::Matrix<T,-1,-1> & psi0, double & sigma,
                                        Eigen::Matrix<double, 1, 2> & vb){
  (void) sigma;
  (void) k;
  (void) psi0;
  (void) vb;
}

template <typename T, unsigned D>
void KPM_Vector<T,D>::build_site(unsigned long pos){
  (void) pos;
}

template <typename T, unsigned D>
void KPM_Vector <T, D>::build_planewave(Eigen::Matrix<double,-1,1> & k, Eigen::Matrix<T,-1,1> & weight){
  (void) k;
  (void) weight;
}

template <typename T, unsigned D>
template < unsigned MULT,bool VELOCITY> 
void KPM_Vector<T,D>::build_regular_phases(int i1, unsigned axis){
  (void) i1 ;
  (void) axis;
}

template <typename T, unsigned D>
template < unsigned MULT> 
void KPM_Vector<T,D>::initiate_stride(std::size_t & istr){}

template <typename T, unsigned D>
template < unsigned MULT> 
void inline KPM_Vector<T,D>::mult_local_disorder(const  std::size_t & j0, const  std::size_t & io){
  (void) j0;
  (void) io;
}

template <typename T, unsigned D>
void inline KPM_Vector<T,D>::mult_regular_hoppings(const  std::size_t & j0, const  std::size_t & io){
  (void) j0;
  (void) io;
}

template <typename T, unsigned D>
void KPM_Vector<T,D>::Velocity(T * phi0,T * phiM1, unsigned axis){
  (void) phi0;
  (void) phiM1;
  (void) axis;
}


template <typename T, unsigned D>
void KPM_Vector<T,D>::Velocity(T * phi0,T * phiM1, int axis){
  (void) phi0;
  (void) phiM1;
  (void) axis;
}

template <typename T, unsigned D>
template <unsigned MULT, bool VELOCITY>
void KPM_Vector<T,D>::KPM_MOTOR (T * phi0a, T * phiM1a, T *phiM2a, unsigned axis){
  (void) phi0a;
  (void) phiM2;
  (void) phiM2a;
  (void) axis;
}

template <typename T, unsigned D>
void KPM_Vector<T,D>::measure_wave_packet(T * bra, T * ket, T * results){
  (void) bra;
  (void) ket;
  (void) results;
}

template <typename T, unsigned D>
void KPM_Vector<T,D>::Exchange_Boundaries(){}

template <typename T, unsigned D>
void KPM_Vector<T,D>::test_boundaries_system(){}

template <typename T, unsigned D>
void KPM_Vector<T,D>::empty_ghosts(int mem_index){(void) mem_index;}

template class KPM_Vector<float ,1u>;
template class KPM_Vector<double ,1u>;
template class KPM_Vector<long double ,1u>;
template class KPM_Vector<std::complex<float> ,1u>;
template class KPM_Vector<std::complex<double> ,1u>;
template class KPM_Vector<std::complex<long double> ,1u>;

template class KPM_Vector<float ,3u>;
template class KPM_Vector<double ,3u>;
template class KPM_Vector<long double ,3u>;
template class KPM_Vector<std::complex<float> ,3u>;
template class KPM_Vector<std::complex<double> ,3u>;
template class KPM_Vector<std::complex<long double> ,3u>;


