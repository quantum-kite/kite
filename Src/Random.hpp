/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

template <typename T> 
class KPMRandom {
  std::mt19937 rng;
  std::uniform_real_distribution <double> dist;
  std::normal_distribution<> gauss;
public:
  
  typedef typename extract_value_type<T>::value_type value_type;
  
  KPMRandom();
  void init_random();
  double get();  
  double uniform(double  mean, double  width);
  double gaussian(double  mean, double  width);
  
  template <typename U = T>
  typename std::enable_if<is_tt<std::complex, U>::value, U>::type initA();
  
  template <typename U = T>
  typename std::enable_if<!is_tt<std::complex, U>::value, U>::type initA();
  
  T init();
  
};

