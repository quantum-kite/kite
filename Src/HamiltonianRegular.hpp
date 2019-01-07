/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

#ifndef DEBUG
#define DEBUG 0
#endif

template <typename T, unsigned D>
struct Periodic_Operator : public ComplexTraits<T> {
  using ComplexTraits<T>::multEiphase;
  typedef typename extract_value_type<T>::value_type         value_type;
  LatticeStructure <D> & r;
  // Non-diagonal component of the operator
  Eigen::Array<unsigned,    Eigen::Dynamic, 1 >            NHoppings;         // Number of elements different from Zero from each orbital
  Eigen::Array<std::ptrdiff_t, Eigen::Dynamic, Eigen::Dynamic> distance;      // Distance in the basis 
  Eigen::Array<   T, Eigen::Dynamic, Eigen::Dynamic> hopping;                 // Hopping
  std::vector<Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>>        v;
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> dist; 
  Periodic_Operator<T,D>(char *, LatticeStructure <D> & );
  void Convert_Build (  LatticeStructure <D> &  );
  void build_velocity(std::vector<unsigned> & components, unsigned n);  
};
