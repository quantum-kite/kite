/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

#include <type_traits>
/*
  Auxiliar code to define specialized methods in templated classes depending on the argument T is complex:

  template <typename U = T>
  typename std::enable_if<is_tt<std::complex, U>::value, U>::type aux_wr(unsigned long x ) {
  CODE
  }; 

  or non-complex:
   template <typename U = T>
   typename std::enable_if<!is_tt<std::complex, U>::value, U>::type aux_wr(unsigned long x ) {
    CODE
  };
  
  Get the template argument of a complex:
  typedef typename extract_value_type<T>::value_type value_type;
 */

template <template <class...> class TT, class... Args>
std::true_type is_tt_impl(TT<Args...>);
template <template <class...> class TT>
std::false_type is_tt_impl(...);

template <template <class...> class TT, class T>
using is_tt = decltype(is_tt_impl<TT>(std::declval<typename std::decay<T>::type>()));

template<typename T>
struct extract_value_type
{
  typedef T value_type;
};

template<template<typename, typename ...> class X, typename T, typename ...Args>
struct extract_value_type<X<T, Args...>>   //specialization                                                                                                                    
{
  typedef T value_type;
};

