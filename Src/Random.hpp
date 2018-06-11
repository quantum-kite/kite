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
  
  KPMRandom() {
    init_random();
  };

  void init_random()
  {
    std::random_device r;
    std::array<int, 624> seed_data;
    std::generate(seed_data.begin(), seed_data.end(), std::ref(r));
    std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
    rng.seed(seq); 
  };

  double get() {
    return dist(rng);
  };
  
  double uniform(double  mean, double  width) {
    // mean  : mean value
    // width : root mean square deviation
    return mean + sqrt(3.) * width * (2 * dist(rng)  - 1);
  };
  
  double gaussian(double  mean, double  width) {
    // mean  : mean value
    // width : root mean square deviation
    return mean + width*gauss(rng);
  };
  
  
  template <typename U = T>
  typename std::enable_if<is_tt<std::complex, U>::value, U>::type init() {
    return exp(T(0., value_type(2*M_PI*dist(rng)) ));
  };
  
  template <typename U = T>
  typename std::enable_if<!is_tt<std::complex, U>::value, U>::type init() {
    return (2*dist(rng) - 1.)*sqrt(3);
  };
};

