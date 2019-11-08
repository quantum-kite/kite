/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

#include <iomanip>
template <typename T,unsigned D>
class KPM_Vector : public KPM_VectorBasis <T, D> {
private:
  std::size_t 	*MemIndBeg[D][2];
  std::size_t 	*MemIndEnd[D][2];
  std::size_t        block[D][2];
  std::size_t          tile[D];
  std::size_t   tile_ghosts[D];
  std::size_t      transf_max[D]; // [d][edged]
  std::size_t transf_bound[D][2]; // [d][edged]
  T          ***mult_t1_ghost_cor;
  T                         *phi0;
  T                        *phiM1;
  T                        *phiM2;
  const std::size_t           std;
public:
  LatticeStructure<D>         & r;
  Hamiltonian<T,D>            & h;
  Coordinates<std::size_t,D+1>  x;

  typedef typename extract_value_type<T>::value_type value_type;
  using KPM_VectorBasis<T,D>::simul;
  using KPM_VectorBasis<T,D>::index;
  using KPM_VectorBasis<T,D>::v;
  using KPM_VectorBasis<T,D>::memory;
  using KPM_VectorBasis<T,D>::aux_wr;
  using KPM_VectorBasis<T,D>::aux_test;
  using KPM_VectorBasis<T,D>::inc_index;
  using KPM_VectorBasis<T,D>::assign_value;
  using KPM_VectorBasis<T,D>::myconj;
  using KPM_VectorBasis<T,D>::multEiphase;
  
  KPM_Vector(int mem, Simulation<T,D> & sim);
  ~KPM_Vector(void);
  void initiate_vector();
  T get_point();

  void build_wave_packet(Eigen::Matrix<double,-1,-1> & k, Eigen::Matrix<T,-1,-1> & psi0, double & sigma,
                         Eigen::Matrix<double, 1, 2> & vb);
  void build_planewave(Eigen::Matrix<double,-1,1> & k, Eigen::Matrix<T,-1,1> & weight);
  void build_site(unsigned long R);

  template < unsigned MULT,bool VELOCITY> 
  void build_regular_phases(int i1, unsigned axis);
  template < unsigned MULT> 
  void initiate_stride(std::size_t & istr);
  template < unsigned MULT> 
  void inline mult_local_disorder(const  std::size_t & j0, const  std::size_t & io);
  void inline mult_regular_hoppings(const  std::size_t & j0, const  std::size_t & io);
  template <unsigned MULT> 
  void Multiply(){}

  void Velocity(T * phi0,T * phiM1, unsigned axis);
  void Velocity(T * phi0,T * phiM1, int axis);
  template <unsigned MULT, bool VELOCITY>
  void KPM_MOTOR(T * phi0a, T * phiM1a, T *phiM2a, unsigned axis);
  void measure_wave_packet(T * bra, T * ket, T * results);  
  void Exchange_Boundaries();
  void test_boundaries_system();
  void empty_ghosts(int mem_index);

};


#include "KPM_Vector2D.hpp"
#include "KPM_Vector3D.hpp"



