#include "ComplexTraits.hpp"

/****** assign *************/
template <typename T>
template <typename U>
typename std::enable_if<!is_tt<std::complex, U>::value, U>::type ComplexTraits <T>::assign_valueB(double x, double y) {
  return T(x +0*y);
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
template <typename U>
typename std::enable_if<is_tt<std::complex, U>::value, U>::type ComplexTraits<T>::multEiphaseB(double phase) {
  return static_cast<U>(exp(std::complex<double>(0,phase)));
};

template <typename T>
template <typename U>
typename std::enable_if<!is_tt<std::complex, U>::value, U>::type ComplexTraits<T>::multEiphaseB(double phase) {
  return T(1.0 + 0*phase);
};

// Define aux_wr for complex T
template <typename T>
template <typename U>
typename std::enable_if<is_tt<std::complex, U>::value, U>::type ComplexTraits<T>::aux_wrB(std::size_t x ) {
  typedef typename extract_value_type<U>::value_type value_type;
  return U(value_type(x), value_type(2*x));
};

// Define aux_wr for non complex T
template <typename T>
template <typename U>
typename std::enable_if<!is_tt<std::complex, U>::value, U>::type ComplexTraits<T>::aux_wrB(std::size_t x) {
  return U(x);
};

template <typename T>
T ComplexTraits<T>::myconj(T & x) {
  return ComplexTraits<T>::myconjB<T>(x);
};

template <typename T>
T ComplexTraits<T>::assign_value(double x, double y) {
  return ComplexTraits<T>::assign_valueB<T>(x,y);
};

template <typename T>
T ComplexTraits<T>::multEiphase(double x) {
  return ComplexTraits<T>::multEiphaseB<T>(x);
};

template <typename T>
T ComplexTraits<T>::aux_wr(std::size_t x) {
  return ComplexTraits<T>::aux_wrB<T>(x);
};






template struct ComplexTraits<double>;
template struct ComplexTraits<float>;
template struct ComplexTraits<long double>;
template struct ComplexTraits<std::complex<double>>;
template struct ComplexTraits<std::complex<float>>;
template struct ComplexTraits<std::complex<long double>>;

