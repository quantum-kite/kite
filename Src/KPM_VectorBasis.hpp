/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

template <typename T, unsigned D>
class KPM_VectorBasis: public ComplexTraits<T> {
protected:
  int index;
  const int memory;
  Simulation<T,D> & simul;  
public:
  using ComplexTraits<T>::assign_value;
  using ComplexTraits<T>::myconj;
  using ComplexTraits<T>::multEiphase;
  using ComplexTraits<T>::aux_wr;
  Eigen::Matrix <T, Eigen::Dynamic,  Eigen::Dynamic > v;  
  KPM_VectorBasis(int mem,  Simulation<T,D> & sim);
  void set_index(int i);
  void inc_index();  
  unsigned get_index();
  bool aux_test(T & x, T & y );  
};
