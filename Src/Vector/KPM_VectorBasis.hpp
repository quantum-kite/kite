/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2021, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/

template <typename T, unsigned D>
class KPM_Vector;

template <typename T, unsigned D>
class KPM_VectorBasis: public ComplexTraits<T> {
protected:
  int index;
  const int memory;
  Simulation<T,D> & simul;
  Hamiltonian<T,D>           & h;
public:
  typedef typename extract_value_type<T>::value_type value_type;
  using ComplexTraits<T>::assign_value;
  using ComplexTraits<T>::myconj;
  using ComplexTraits<T>::multEiphase;
  using ComplexTraits<T>::aux_wr;
  Eigen::Matrix <T, Eigen::Dynamic,  Eigen::Dynamic > v;  
  KPM_VectorBasis(int mem,  Simulation<T,D> & sim);
  void     set_index(int i);
  void     inc_index();
  void     dec_index();
  unsigned get_index();
  bool     aux_test(T & x, T & y );
  template <unsigned MULT>
  void     Multiply();
  void     Velocity(KPM_Vector<T,D> * kpm_final, std::vector<std::vector<unsigned>> & indices, int axis);
  void     cheb_iteration(unsigned );
  
  template <unsigned MULT, bool VELOCITY>  
  void     multiply_defect(std::size_t istr, T* & phi0, T* & phiM1, unsigned axis);
  void     build_defect_planewave(Eigen::Matrix<double,-1,1> & k , Eigen::Matrix<T,-1,1> & weight );
};
