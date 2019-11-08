/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

#include <vector>
#include <Eigen/Dense>
#include "Global.hpp"

template <typename T>
GLOBAL_VARIABLES<T>::GLOBAL_VARIABLES() { }

template <typename T>
void GLOBAL_VARIABLES<T>::addbond( std::size_t  ele1, std::ptrdiff_t ele2, T hop ) {
  element1.push_back(ele1);
  element2_diff.push_back(ele2); 
  hopping.push_back(hop);
}

template <typename T>
void GLOBAL_VARIABLES<T>::addlocal( std::size_t  ele,  T u ) {
  element.push_back(ele);
  U.push_back(u);
}

template struct GLOBAL_VARIABLES<float>;
template struct GLOBAL_VARIABLES<double>;
template struct GLOBAL_VARIABLES<long double>;
template struct GLOBAL_VARIABLES<std::complex<float>>;
template struct GLOBAL_VARIABLES<std::complex<double>>;
template struct GLOBAL_VARIABLES<std::complex<long double>>;
