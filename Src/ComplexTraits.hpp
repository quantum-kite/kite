#include <type_traits>

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

