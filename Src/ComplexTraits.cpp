#include "ComplexTraits.hpp"

/****** assign *************/
template <typename T>
template <typename U>
typename std::enable_if<!is_tt<std::complex, U>::value, U>::type ComplexTraits <T>::assign_valueB(double x, double y) {
  return T(x);
};

template <typename T>
template <typename U>
typename std::enable_if<is_tt<std::complex, U>::value, U>::type ComplexTraits<T>::assign_valueB(double x, double y) {
  typedef typename extract_value_type<T>::value_type value_type;
  return T(value_type(x),value_type(y));
};

/****** myconj *************/
template <typename T>
template <typename U>
typename std::enable_if<!is_tt<std::complex, U>::value, U>::type ComplexTraits<T>::myconjB(T & x) {
  return x;
};

template <typename T>
template <typename U>
typename std::enable_if<is_tt<std::complex, U>::value, U>::type ComplexTraits<T>::myconjB(T & x) {
  return std::conj(x);
};

template <typename T>
T ComplexTraits<T>::myconj(T & x) {
  return ComplexTraits<T>::myconjB<T>(x);
};

template <typename T>
T ComplexTraits<T>::assign_value(double x, double y) {
  return ComplexTraits<T>::assign_valueB<T>(x,y);
};


template struct ComplexTraits<double>;
template struct ComplexTraits<float>;
template struct ComplexTraits<long double>;
template struct ComplexTraits<std::complex<double>>;
template struct ComplexTraits<std::complex<float>>;
template struct ComplexTraits<std::complex<long double>>;

