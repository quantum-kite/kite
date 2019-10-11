
/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

#include <iomanip>
template <typename T>
class KPM_Vector <T, 3> : public KPM_VectorBasis <T,3> {
private:
  static const unsigned             D=3u;
  std::size_t  	        *MemIndBeg[3][2];   // Corner from the postion read the boundary for each Orbital
  std::size_t  	        *MemIndEnd[3][2];   // Corner from the position to write the boundary for each Orbital
  std::size_t                block[3][2];
  std::size_t                  tile[3];
  std::size_t           tile_ghosts[3];
  std::size_t           transf_max[3][3]; // [d][Maximum lateral bondary lengh]
  std::size_t      transf_bound[3][2][3]; // [d][edged][Lateral boundary lengh]
  T                 ***mult_t1_ghost_cor;
  T                                *phi0;
  T                               *phiM1;
  T                               *phiM2;
public:
  LatticeStructure<3u>               & r;
  Hamiltonian<T,3u>                  & h;
  //  Coordinates<std::size_t,4>           x;
  typedef typename extract_value_type<T>::value_type value_type;
  using KPM_VectorBasis<T,3>::simul;
  using KPM_VectorBasis<T,3>::index;
  using KPM_VectorBasis<T,3>::v;
  using KPM_VectorBasis<T,3>::memory;
  using KPM_VectorBasis<T,3>::aux_wr;
  using KPM_VectorBasis<T,3>::aux_test;
  using KPM_VectorBasis<T,3>::inc_index;
  using KPM_VectorBasis<T,3>::assign_value;
  using KPM_VectorBasis<T,3>::myconj;
  using KPM_VectorBasis<T,3>::multEiphase;
  
  KPM_Vector(int mem, Simulation<T,3> & sim);
  ~KPM_Vector(void);
  void initiate_vector();
  T get_point();
  void build_wave_packet(Eigen::Matrix<double,-1,-1> & k, Eigen::Matrix<T,-1,-1> & psi0, double & sigma,
                         Eigen::Matrix<double,1,2> & vb);
  template < unsigned MULT,bool VELOCITY> 
  void build_regular_phases(int i1, unsigned axis);
  template < unsigned MULT> 
  void initiate_stride(std::size_t & istr);
  template < unsigned MULT> 
  void inline mult_local_disorder(const  std::size_t & j0, const  std::size_t & io);
  void inline mult_regular_hoppings(const  std::size_t & j0, const  std::size_t & io);
  template <unsigned MULT> 
  void Multiply();
  void Velocity(T * phi0,T * phiM1, unsigned axis);
  void Velocity(T * phi0,T * phiM1, int);
  template <unsigned MULT, bool VELOCITY>
  void KPM_MOTOR(T * phi0a, T * phiM1a, T *phiM2a, unsigned axis);
  void measure_wave_packet(T * bra, T * ket, T * results);  
  void Exchange_Boundaries();
  void test_boundaries_system();
  void empty_ghosts(int mem_index);
  void build_site(unsigned long R);
  void build_planewave(Eigen::Matrix<double,-1,1> & k, Eigen::Matrix<T,-1,1> & weight);
};
      
