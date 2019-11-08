#include "Generic.hpp"
#include "ComplexTraits.hpp"
#include "Random.hpp"

template <typename T>
KPMRandom<T>::KPMRandom() {
  init_random();
}

template <typename T>
void KPMRandom<T>::init_random()
{

  char *env;
  env = getenv("SEED");
    if(env==NULL){
      // Didn't find the seed
      std::random_device r;
      std::array<int, 624> seed_data;
      std::generate(seed_data.begin(), seed_data.end(), std::ref(r));
      std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
      rng.seed(seq); 
    }
    else {
      // Found the seed
      rng.seed(atoi(env)); 
    }
}

template <typename T>
double  KPMRandom<T>::get() {
  return dist(rng);
}
template <typename T>
double KPMRandom<T>::uniform(double  mean, double  width) {
    // mean  : mean value
    // width : root mean square deviation
  return mean + sqrt(3.) * width * (2 * dist(rng)  - 1);
}
template <typename T>
double KPMRandom<T>::gaussian(double  mean, double  width) {
  // mean  : mean value
  // width : root mean square deviation
  return mean + width*gauss(rng);
}

template <typename T>
template <typename U>
typename std::enable_if<is_tt<std::complex, U>::value, U>::type KPMRandom<T>::initA() {
  return exp(T(0., value_type(2*M_PI*dist(rng)) ));
}

template <typename T>
template <typename U>
typename std::enable_if<!is_tt<std::complex, U>::value, U>::type KPMRandom<T>::initA() {
  return (2*dist(rng) - 1.)*sqrt(3);
}

template <typename T>
T KPMRandom<T>::init(){
  return initA<T>();
}
  
template class KPMRandom<float>;
template class KPMRandom<double>;
template class KPMRandom<long double>;

template class KPMRandom<std::complex<float>>;
template class KPMRandom<std::complex<double>>;
template class KPMRandom<std::complex<long double>>;
