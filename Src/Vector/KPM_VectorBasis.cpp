/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2021, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/




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
KPM_VectorBasis<T,D>::KPM_VectorBasis(int mem,  Simulation<T,D> & sim) : memory(mem), simul(sim), h(sim.h)
{
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
void KPM_VectorBasis<T,D>::dec_index() {
  index = (index + memory -1) % memory;
}

template<typename T, unsigned D>
unsigned KPM_VectorBasis<T,D>::get_index() {
  return index;
}

template<typename T, unsigned D>
bool KPM_VectorBasis<T,D>::aux_test(T & x, T & y ) {
  return (abs(x - y) > std::numeric_limits<double>::epsilon());
}

template <typename T, unsigned D>
template <unsigned MULT> 
void KPM_VectorBasis<T,D>::Multiply() {
  vverbose_message("Entered Multiply");
  typedef KPM_Vector<T,D> myvector;
  unsigned i = 0;
  /*
    Mosaic Multiplication using a TILE of TILE x TILE
    Right Now We expect that both ld[0] and ld[1]  are multiple of TILE
    MULT = 0 : For the case of the Velocity/Hamiltonian
    MULT = 1 : For the case of the KPM_iteration
  */
  inc_index();
  KPM_Vector<T,D>* child = static_cast<KPM_Vector<T,D>*>(this);
  child->template KPM_MOTOR<MULT, false>(child, i);
}

template<typename T, unsigned D>
void KPM_VectorBasis<T,D>::Velocity(KPM_Vector<T,D> * kpm_final,  std::vector<std::vector<unsigned>> & indices, int pos)
{
  // kpm_final shoud be different from this instance
  // 
  KPM_Vector<T,D>* child = static_cast<KPM_Vector<T,D>*>(this);
  switch(indices.at(pos).size())
    {
    case 0:
      break;
    default:
      {
        inc_index(); 
        child->template KPM_MOTOR<0u, true>(kpm_final, static_cast<unsigned>(pos));
        dec_index();
      }
    }
}

template<typename T, unsigned D>
void KPM_VectorBasis<T,D>::cheb_iteration(unsigned n)
{
  switch(n)
    {
    case 0:
      break;
    case 1:
      this->template Multiply<0>();
      break;
    default:
      this->template Multiply<1>();
    }
}

template <typename T,unsigned D> 
template <unsigned MULT, bool VELOCITY>  
void KPM_VectorBasis<T,D>::multiply_defect(std::size_t istr, T* & phi0, T* & phiM1, unsigned axis)
{
  Coordinates<std::ptrdiff_t, D + 1> local1(simul.r.Ld);
  
  for(auto id = h.hd.begin(); id != h.hd.end(); id++)
    for(std::size_t i = 0; i <  id->position.at(istr).size(); i++)
      {
        
        // Add deffect Hoppings between lattice sites 
        std::size_t ip = id->position.at(istr)[i];
        std::size_t iv = local1.set_coord(ip).coord[D - 1];
        for(unsigned k = 0; k < id->hopping.size(); k++)
          {
            std::size_t k1 = ip + id->node_position[id->element1[k]];
            std::size_t k2 = ip + id->node_position[id->element2[k]];
            
            if( VELOCITY )
              phi0[k1] += value_type(MULT + 1) * id->v.at(axis).at(k) * id-> new_hopping(k, iv) * phiM1[k2] ;
            else
              phi0[k1] += value_type(MULT + 1) * id->new_hopping(k, iv) * phiM1[k2] ;
          }
        // Add deffect lattice local energies  
        if(!VELOCITY)
          for(std::size_t k = 0; k < id->U.size(); k++)
            {
              std::size_t k1 = ip + id->node_position[id->element[k]];
              phi0[k1] += value_type(MULT + 1) * id->U[k] * phiM1[k1];
            }        
      }

}

template <typename T,unsigned D>
void KPM_VectorBasis<T,D>::build_defect_planewave(Eigen::Matrix<double,-1,1> & k , Eigen::Matrix<T,-1,1> & weight )
{
}


// Instantiate KPM_VectorBasis
#define instantiate(type, dim)               template class KPM_VectorBasis <type,dim>; \
  template void KPM_VectorBasis<type,dim>::template multiply_defect<0u,false>(std::size_t , type* & , type * & , unsigned ); \
  template void KPM_VectorBasis<type,dim>::template multiply_defect<1u,false>(std::size_t , type* & , type * & , unsigned ); \
  template void KPM_VectorBasis<type,dim>::template multiply_defect<0u,true>(std::size_t , type* & , type * & , unsigned ); \
  template void KPM_VectorBasis<type,dim>::template Multiply<0u>(); \
  template void KPM_VectorBasis<type,dim>::template Multiply<1u>();
#include "instantiate.hpp"
