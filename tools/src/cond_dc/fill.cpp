//Eigen::Matrix<std::complex<T>, -1, -1, Eigen::RowMajor> 
//fill_dgreenR(int NEnergies, int NMoments, double scat){

//}


#include <Eigen/Dense>
#include <complex>
#include "H5Cpp.h"
#include "../parse_input.hpp"
#include "../systemInfo.hpp"
#include "conductivity_dc.hpp"
#include "../functions.hpp"
#include <fstream>


template <typename U, unsigned DIM>
Eigen::Matrix<std::complex<U>, -1, -1, Eigen::ColMajor> conductivity_dc<U, DIM>::fill_delta(){

  Eigen::Matrix<std::complex<U>, -1, -1, Eigen::ColMajor> greenR;
  greenR = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(Moments_G, NEnergies);

  // Imaginary part of the Green's function: Dirac delta
  U factor;
  std::complex<U> complexEnergy;
  for(int i = 0; i < NEnergies; i++){
    complexEnergy = std::complex<U>(energies(i), deltascat);
    for(int m = 0; m < Moments_G; m++){
      factor = -1.0/(1.0 + U(m==0))/M_PI;
      greenR(m, i) = green(m, 1, complexEnergy).imag()*factor;
    }
  }
  return greenR;
}

template <typename U, unsigned DIM>
Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> conductivity_dc<U, DIM>::fill_dgreenR(){
  std::complex<U> complexEnergyP;

  Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> dgreenR;
  dgreenR = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(NEnergies, Moments_D);

  U factor;
  for(int i = 0; i < NEnergies; i++){
    complexEnergyP = std::complex<U>(energies(i), scat);
    for(int m = 0; m < Moments_D; m++){
      factor = 1.0/(1.0 + U(m==0));
      dgreenR(i, m) = dgreen<U>(m,  1, complexEnergyP)*factor;
    }
  }
  return dgreenR;
}


template <typename U, unsigned DIM>
Eigen::Matrix<std::complex<U>, -1, -1> conductivity_dc<U, DIM>::triple_product(
  Eigen::Matrix<std::complex<U>, -1, -1, Eigen::ColMajor> greenR,
  Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> dgreenR){

  // GammaE has NE elements
  Eigen::Array<std::complex<U>, -1, -1> GammaE;
  GammaE = Eigen::Array<std::complex<U>, -1, -1>::Zero(NEnergies, 1);

  omp_set_num_threads(NumThreads);
#pragma omp parallel firstprivate(greenR)
  {

#pragma omp for schedule(static, 1) nowait
    for(int thread = 0; thread < NumThreads; thread++){
      int localMoments = Moments_D/NumThreads;

      Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> LocalGamma;
      Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> local_dgreenR;
      Eigen::Matrix<std::complex<U>, -1, -1, Eigen::ColMajor> GammaEN;
      Eigen::Array<std::complex<U>, -1, -1> LocalGammaE;

      LocalGamma = Gamma_Padded.matrix().block(thread*localMoments, 0, localMoments, Moments_G);
      local_dgreenR = dgreenR.matrix().block(0, thread*localMoments, NEnergies, localMoments);

      // GammaEN has NE * ND/T elements
      GammaEN = LocalGamma*greenR;

      // LocalGammaE has NE elements
      LocalGammaE = Eigen::Array<std::complex<U>, -1, -1>::Zero(NEnergies, 1);
      for(int i = 0; i < NEnergies; i++)
        LocalGammaE(i) = std::complex<U>(local_dgreenR.row(i)*GammaEN.col(i));

#pragma omp critical
      {
      GammaE += 2*LocalGammaE.imag();
      }
    }
#pragma omp barrier
  }
  return GammaE;

}



template <typename U, unsigned DIM>
Eigen::Matrix<std::complex<U>, -1, 1> conductivity_dc<U, DIM>::calc_cond(Eigen::Matrix<std::complex<U>, -1, -1> GammaE){

  Eigen::Matrix<std::complex<U>, -1, 1> condDC;
  Eigen::Matrix<std::complex<U>, -1, 1> integrand;


  condDC = Eigen::Matrix<std::complex<U>, -1, 1>::Zero(NFermiEnergies, 1);
  integrand = Eigen::Matrix<std::complex<U>, -1, 1>::Zero(NEnergies, 1);

  U fermi;
  for(int i = 0; i < NFermiEnergies; i++){
    fermi = fermiEnergies(i);
    for(int j = 0; j < NEnergies; j++){
      integrand(j) = GammaE(j)*fermi_function(energies(j), fermi, beta);
    }
    condDC(i) = integrate(energies, integrand);
  }

  return condDC;
}



template <typename U, unsigned DIM>
void conductivity_dc<U, DIM>::save_to_file(Eigen::Matrix<std::complex<U>, -1, -1> condDC){

  std::complex<U> cond;
  U energy;
  std::ofstream myfile;
  myfile.open(filename);
  for(int i=0; i < NFermiEnergies; i++){
    energy = fermiEnergies(i)*systemInfo.energy_scale + systemInfo.energy_shift;
    cond = condDC(i);
    myfile  << energy << " " << cond.real() << " " << cond.imag() << "\n";
  }
  
  myfile.close();


}

template Eigen::Matrix<std::complex<float>, -1, -1, Eigen::ColMajor> conductivity_dc<float, 1u>::fill_delta();
template Eigen::Matrix<std::complex<float>, -1, -1, Eigen::ColMajor> conductivity_dc<float, 2u>::fill_delta();
template Eigen::Matrix<std::complex<float>, -1, -1, Eigen::ColMajor> conductivity_dc<float, 3u>::fill_delta();

template Eigen::Matrix<std::complex<double>, -1, -1, Eigen::ColMajor> conductivity_dc<double, 1u>::fill_delta();
template Eigen::Matrix<std::complex<double>, -1, -1, Eigen::ColMajor> conductivity_dc<double, 2u>::fill_delta();
template Eigen::Matrix<std::complex<double>, -1, -1, Eigen::ColMajor> conductivity_dc<double, 3u>::fill_delta();

template Eigen::Matrix<std::complex<long double>, -1, -1, Eigen::ColMajor> conductivity_dc<long double, 1u>::fill_delta();
template Eigen::Matrix<std::complex<long double>, -1, -1, Eigen::ColMajor> conductivity_dc<long double, 2u>::fill_delta();
template Eigen::Matrix<std::complex<long double>, -1, -1, Eigen::ColMajor> conductivity_dc<long double, 3u>::fill_delta();






template Eigen::Matrix<std::complex<float>, -1, -1, Eigen::RowMajor> conductivity_dc<float, 1u>::fill_dgreenR();
template Eigen::Matrix<std::complex<float>, -1, -1, Eigen::RowMajor> conductivity_dc<float, 2u>::fill_dgreenR();
template Eigen::Matrix<std::complex<float>, -1, -1, Eigen::RowMajor> conductivity_dc<float, 3u>::fill_dgreenR();

template Eigen::Matrix<std::complex<double>, -1, -1, Eigen::RowMajor> conductivity_dc<double, 1u>::fill_dgreenR();
template Eigen::Matrix<std::complex<double>, -1, -1, Eigen::RowMajor> conductivity_dc<double, 2u>::fill_dgreenR();
template Eigen::Matrix<std::complex<double>, -1, -1, Eigen::RowMajor> conductivity_dc<double, 3u>::fill_dgreenR();

template Eigen::Matrix<std::complex<long double>, -1, -1, Eigen::RowMajor> conductivity_dc<long double, 1u>::fill_dgreenR();
template Eigen::Matrix<std::complex<long double>, -1, -1, Eigen::RowMajor> conductivity_dc<long double, 2u>::fill_dgreenR();
template Eigen::Matrix<std::complex<long double>, -1, -1, Eigen::RowMajor> conductivity_dc<long double, 3u>::fill_dgreenR();




template Eigen::Matrix<std::complex<float>, -1, -1> conductivity_dc<float, 1u>::triple_product(Eigen::Matrix<std::complex<float>, -1, -1, Eigen::ColMajor>, Eigen::Matrix<std::complex<float>, -1, -1, Eigen::RowMajor>);
template Eigen::Matrix<std::complex<float>, -1, -1> conductivity_dc<float, 2u>::triple_product(Eigen::Matrix<std::complex<float>, -1, -1, Eigen::ColMajor>, Eigen::Matrix<std::complex<float>, -1, -1, Eigen::RowMajor>);
template Eigen::Matrix<std::complex<float>, -1, -1> conductivity_dc<float, 3u>::triple_product(Eigen::Matrix<std::complex<float>, -1, -1, Eigen::ColMajor>, Eigen::Matrix<std::complex<float>, -1, -1, Eigen::RowMajor>);

template Eigen::Matrix<std::complex<double>, -1, -1> conductivity_dc<double, 1u>::triple_product(Eigen::Matrix<std::complex<double>, -1, -1, Eigen::ColMajor>, Eigen::Matrix<std::complex<double>, -1, -1, Eigen::RowMajor>);
template Eigen::Matrix<std::complex<double>, -1, -1> conductivity_dc<double, 2u>::triple_product(Eigen::Matrix<std::complex<double>, -1, -1, Eigen::ColMajor>, Eigen::Matrix<std::complex<double>, -1, -1, Eigen::RowMajor>);
template Eigen::Matrix<std::complex<double>, -1, -1> conductivity_dc<double, 3u>::triple_product(Eigen::Matrix<std::complex<double>, -1, -1, Eigen::ColMajor>, Eigen::Matrix<std::complex<double>, -1, -1, Eigen::RowMajor>);

template Eigen::Matrix<std::complex<long double>, -1, -1> conductivity_dc<long double, 1u>::triple_product(Eigen::Matrix<std::complex<long double>, -1, -1, Eigen::ColMajor>, Eigen::Matrix<std::complex<long double>, -1, -1, Eigen::RowMajor>);
template Eigen::Matrix<std::complex<long double>, -1, -1> conductivity_dc<long double, 2u>::triple_product(Eigen::Matrix<std::complex<long double>, -1, -1, Eigen::ColMajor>, Eigen::Matrix<std::complex<long double>, -1, -1, Eigen::RowMajor>);
template Eigen::Matrix<std::complex<long double>, -1, -1> conductivity_dc<long double, 3u>::triple_product(Eigen::Matrix<std::complex<long double>, -1, -1, Eigen::ColMajor>, Eigen::Matrix<std::complex<long double>, -1, -1, Eigen::RowMajor>);




template Eigen::Matrix<std::complex<float>, -1, 1> conductivity_dc<float, 1u>::calc_cond(Eigen::Matrix<std::complex<float>, -1, -1>);
template Eigen::Matrix<std::complex<float>, -1, 1> conductivity_dc<float, 2u>::calc_cond(Eigen::Matrix<std::complex<float>, -1, -1>);
template Eigen::Matrix<std::complex<float>, -1, 1> conductivity_dc<float, 3u>::calc_cond(Eigen::Matrix<std::complex<float>, -1, -1>);

template Eigen::Matrix<std::complex<double>, -1, 1> conductivity_dc<double, 1u>::calc_cond(Eigen::Matrix<std::complex<double>, -1, -1>);
template Eigen::Matrix<std::complex<double>, -1, 1> conductivity_dc<double, 2u>::calc_cond(Eigen::Matrix<std::complex<double>, -1, -1>);
template Eigen::Matrix<std::complex<double>, -1, 1> conductivity_dc<double, 3u>::calc_cond(Eigen::Matrix<std::complex<double>, -1, -1>);

template Eigen::Matrix<std::complex<long double>, -1, 1> conductivity_dc<long double, 1u>::calc_cond(Eigen::Matrix<std::complex<long double>, -1, -1>);
template Eigen::Matrix<std::complex<long double>, -1, 1> conductivity_dc<long double, 2u>::calc_cond(Eigen::Matrix<std::complex<long double>, -1, -1>);
template Eigen::Matrix<std::complex<long double>, -1, 1> conductivity_dc<long double, 3u>::calc_cond(Eigen::Matrix<std::complex<long double>, -1, -1>);



template void conductivity_dc<float, 1u>::save_to_file(Eigen::Matrix<std::complex<float>, -1, -1>);
template void conductivity_dc<float, 2u>::save_to_file(Eigen::Matrix<std::complex<float>, -1, -1>);
template void conductivity_dc<float, 3u>::save_to_file(Eigen::Matrix<std::complex<float>, -1, -1>);

template void conductivity_dc<double, 1u>::save_to_file(Eigen::Matrix<std::complex<double>, -1, -1>);
template void conductivity_dc<double, 2u>::save_to_file(Eigen::Matrix<std::complex<double>, -1, -1>);
template void conductivity_dc<double, 3u>::save_to_file(Eigen::Matrix<std::complex<double>, -1, -1>);

template void conductivity_dc<long double, 1u>::save_to_file(Eigen::Matrix<std::complex<long double>, -1, -1>);
template void conductivity_dc<long double, 2u>::save_to_file(Eigen::Matrix<std::complex<long double>, -1, -1>);
template void conductivity_dc<long double, 3u>::save_to_file(Eigen::Matrix<std::complex<long double>, -1, -1>);
