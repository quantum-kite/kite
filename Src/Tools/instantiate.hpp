#define instantiateTYPE(type)                   \
  instantiate(type,1u)                          \
  instantiate(type,2u)                          \
  instantiate(type,3u)

instantiateTYPE(float)
instantiateTYPE(double)
instantiateTYPE(long double)
instantiateTYPE(std::complex<float>)
instantiateTYPE(std::complex<double>)
instantiateTYPE(std::complex<long double>)
