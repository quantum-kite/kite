template <typename T> 
class KPMRandom {
  std::mt19937 rng;
  std::uniform_real_distribution <double> dist;
public:
  
  typedef typename extract_value_type<T>::value_type value_type;
  
  KPMRandom() {
    std::random_device r;
    std::array<int, 624> seed_data;
    std::generate(seed_data.begin(), seed_data.end(), std::ref(r));
    std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
    rng.seed(seq);
  };
  
  double uniform(double  mean, double  width) {
    // mean  : mean value
    // width : root mean square deviation
    return mean + sqrt(3.) * width * (2 * dist(rng)  - 1);
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

