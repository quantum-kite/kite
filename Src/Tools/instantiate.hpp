/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2022, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/

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
