/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

#include <Eigen/Dense>
#include <string>
#include <vector>
#include <iostream>
#include "H5Cpp.h"
#include "../parse_input.hpp"
#include "../systemInfo.hpp"
#include "../functions.hpp"
#include "conductivity_2order.hpp"
#include <omp.h>


template <typename U, unsigned DIM>
Eigen::Matrix<std::complex<U>, -1, -1> conductivity_nonlinear<U, DIM>::Gamma3shgContract_RA(){
  // Calculate the first of the three three-velocity terms.
  // this is the term with two Green's functions that depend on the 
  // frequency, making it more difficult to calculate.
  

  // Number of moments in each direction
  int N0 = NumMoments;
  int N1 = NumMoments;
  int N2 = NumMoments;

  // number of threads, thread number and number of moments alocated
  // to each thread. These have to be computed inside a threaded block
  int N_threads;
  int thread_num;
  int local_NumMoments;

  // Functions that are going to be used by the contractor
  int NumMoments1 = NumMoments; U beta1 = beta; U e_fermi1 = e_fermi;
  std::function<U(int, U)> deltaF = [beta1, e_fermi1, NumMoments1](int n, U energy)->U{
    return delta(n, energy)*U(1.0/(1.0 + U(n==0)))*fermi_function(energy, e_fermi1, beta1)*kernel_jackson<U>(n, NumMoments1);
  };
  
  Eigen::Matrix<std::complex<U>, -1, -1> global_omega_energies;
  Eigen::Matrix<std::complex<U>, -1, -1> omega_energies;
  global_omega_energies = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  omega_energies = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  
  omp_set_num_threads(systemInfo.NumThreads);
  // Start the parallelization. It is done in the direction p
#pragma omp parallel shared(N_threads, global_omega_energies) firstprivate(omega_energies)
{
  
#pragma omp master
{
  N_threads = omp_get_num_threads();
  // check if each thread will get the same number of moments
  if(N2%N_threads != 0){
    std::cout << "The number of Chebyshev moments in the nonlinear optical conductivity must"
      "be a multiple of the number of threads\n" << std::flush;
    exit(1);
  }
}
#pragma omp barrier

#pragma omp for schedule(static, 1) nowait
  for(int i = 0; i < N_threads; i++){
    local_NumMoments = N2/N_threads;
    
    Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> Gamma3Aligned;
    Gamma3Aligned = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>::Zero(N0*local_NumMoments, N1);
    for(long int p = i*local_NumMoments; p < local_NumMoments*(i+1); p++)
      for(long int m = 0; m < N1; m++)
        for(long int n = 0; n < N0; n++)
          Gamma3Aligned((p-i*local_NumMoments)*N0 + n, m) = Gamma3(N0*N1*p + N0*m + n);
      
    // Delta matrix of chebyshev moments and energies
    Eigen::Matrix<std::complex<U>,-1, -1> DeltaMatrix;
    
    DeltaMatrix = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::ColMajor>::Zero(N1, N_energies);
    for(int n = 0; n < NumMoments; n++)
      for(int e = 0; e < N_energies; e++)
        DeltaMatrix(n,e) = deltaF(n, energies(e)); 

    Eigen::Matrix<std::complex<U>, -1, -1> Gamma3NNE;
    Gamma3NNE = Gamma3Aligned*DeltaMatrix;
    
    
    // Matrix of Green's functions
    Eigen::Matrix<std::complex<U>, -1, -1> GreenR, GreenA;
    GreenR  = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>::Zero(N_energies, N0);
    GreenA  = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>::Zero(N_energies, local_NumMoments);

    
    
    U w1, w2;
    for(int w = 0; w < N_omegas; w++){
      w1 = frequencies2(w,0);
      w2 = frequencies2(w,1);
      
      // The scat term is the same in both cases because greenRscat and greenAscat already
      // take into account that the sign of scat is different in those cases
      for(int n = 0; n < N0; n++)
        for(int e = 0; e < N_energies; e++)
          GreenR(e, n) = greenRscat<U>(scat)(n, energies(e) + w1); 
      
      for(int p = 0; p < local_NumMoments; p++)
        for(int e = 0; e < N_energies; e++)
          GreenA(e, p) = greenAscat<U>(scat)(i*local_NumMoments + p, energies(e) - w2); 
      
      Eigen::Matrix<std::complex<U>, -1, -1> temp;
      temp = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(1,1);
      for(int col = 0; col < N_energies; col++){
        for(int p = 0; p < local_NumMoments; p++){
          temp = GreenR.row(col)*Gamma3NNE.block(p*N0, col, N0, 1);
          omega_energies(col, w) += temp(0,0)*GreenA(col, p);
        }
      }
    }
  }
#pragma omp critical
      global_omega_energies += omega_energies;
#pragma omp barrier
}

  return global_omega_energies;

}


template <typename U, unsigned DIM>
Eigen::Matrix<std::complex<U>, -1, -1> conductivity_nonlinear<U, DIM>::Gamma3shgContract_RR(){
  // Calculate the first of the three three-velocity terms.
  // this is the term with two Green's functions that depend on the 
  // frequency, making it more difficult to calculate.
  

  // Number of moments in each direction
  int N0 = NumMoments;
  int N1 = NumMoments;
  int N2 = NumMoments;

  // number of threads, thread number and number of moments alocated
  // to each thread. These have to be computed inside a threaded block
  int N_threads;
  int thread_num;
  int local_NumMoments;

  // Functions that are going to be used by the contractor
  int NumMoments1 = NumMoments; U beta1 = beta; U e_fermi1 = e_fermi;
  std::function<U(int, U)> deltaF = [beta1, e_fermi1, NumMoments1](int n, U energy)->U{
    return delta(n, energy)*U(1.0/(1.0 + U(n==0)))*fermi_function(energy, e_fermi1, beta1)*kernel_jackson<U>(n, NumMoments1);
  };
  
  Eigen::Matrix<std::complex<U>, -1, -1> global_omega_energies;
  Eigen::Matrix<std::complex<U>, -1, -1> omega_energies;
  global_omega_energies = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  omega_energies = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  
  omp_set_num_threads(systemInfo.NumThreads);
  // Start the parallelization. It is done in the direction p
#pragma omp parallel shared(N_threads, global_omega_energies) firstprivate(omega_energies)
{
  
#pragma omp master
{
  N_threads = omp_get_num_threads();
  // check if each thread will get the same number of moments
  if(N0%N_threads != 0){
    std::cout << "The number of Chebyshev moments in the nonlinear optical conductivity must"
      "be a multiple of the number of threads\n" << std::flush;
    exit(1);
  }
}
#pragma omp barrier

#pragma omp for schedule(static, 1) nowait
  for(int i = 0; i < N_threads; i++){
    local_NumMoments = N0/N_threads;
    
    Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> Gamma3Aligned;
    Gamma3Aligned = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>::Zero(N1*local_NumMoments, N2);
    for(long int n = i*local_NumMoments; n < local_NumMoments*(i+1); n++)
      for(long int p = 0; p < N2; p++)
        for(long int m = 0; m < N1; m++)
          Gamma3Aligned((n-i*local_NumMoments)*N1 + m, p) = Gamma3(N0*N1*p + N0*m + n);
      
    // Delta matrix of chebyshev moments and energies
    Eigen::Matrix<std::complex<U>,-1, -1> DeltaMatrix;
    
    DeltaMatrix = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::ColMajor>::Zero(N2, N_energies);
    for(int n = 0; n < NumMoments; n++)
      for(int e = 0; e < N_energies; e++)
        DeltaMatrix(n,e) = deltaF(n, energies(e)); 

    Eigen::Matrix<std::complex<U>, -1, -1> Gamma3NNE;
    Gamma3NNE = Gamma3Aligned*DeltaMatrix;
    
    
    // Matrix of Green's functions
    Eigen::Matrix<std::complex<U>, -1, -1> Green2R, GreenR;
    Green2R = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>::Zero(N_energies, local_NumMoments);
    GreenR  = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>::Zero(N_energies, N1);

    
    
    U w1, w2;
    for(int w = 0; w < N_omegas; w++){
      w1 = frequencies2(w,0);
      w2 = frequencies2(w,1);
      //std::cout << "w1: " << w1 << " w2: " << w2 << "\n";
      
      
      // The scat term is the same in both cases because greenRscat and greenAscat already
      // take into account that the sign of scat is different in those cases
      for(int m = 0; m < N1; m++)
        for(int e = 0; e < N_energies; e++)
          GreenR(e, m) = greenRscat<U>(scat)(m, energies(e) + w2); 
      
      for(int n = 0; n < local_NumMoments; n++)
        for(int e = 0; e < N_energies; e++)
          Green2R(e, n) = greenRscat<U>(2.0*scat)(i*local_NumMoments + n, energies(e) + w1 + w2); 
      
      Eigen::Matrix<std::complex<U>, -1, -1> temp;
      temp = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(1,1);
      for(int col = 0; col < N_energies; col++){
        for(int n = 0; n < local_NumMoments; n++){
          temp = GreenR.row(col)*Gamma3NNE.block(n*N1, col, N1, 1);
          omega_energies(col, w) += temp(0,0)*Green2R(col, n);
        }
      }
    }
  }
#pragma omp critical
      global_omega_energies += omega_energies;
#pragma omp barrier
}

  return global_omega_energies;

}


template <typename U, unsigned DIM>
Eigen::Matrix<std::complex<U>, -1, -1> conductivity_nonlinear<U, DIM>::Gamma3shgContract_AA(){
  // Calculate the first of the three three-velocity terms.
  // this is the term with two Green's functions that depend on the 
  // frequency, making it more difficult to calculate.
  

  // Number of moments in each direction
  int N0 = NumMoments;
  int N1 = NumMoments;
  int N2 = NumMoments;

  // number of threads, thread number and number of moments alocated
  // to each thread. These have to be computed inside a threaded block
  int N_threads;
  int thread_num;
  int local_NumMoments;

  // Functions that are going to be used by the contractor
  int NumMoments1 = NumMoments; U beta1 = beta; U e_fermi1 = e_fermi;
  std::function<U(int, U)> deltaF = [beta1, e_fermi1, NumMoments1](int n, U energy)->U{
    return delta(n, energy)*U(1.0/(1.0 + U(n==0)))*fermi_function(energy, e_fermi1, beta1)*kernel_jackson<U>(n, NumMoments1);
  };
  
  Eigen::Matrix<std::complex<U>, -1, -1> global_omega_energies;
  Eigen::Matrix<std::complex<U>, -1, -1> omega_energies;
  global_omega_energies = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  omega_energies = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  
  omp_set_num_threads(systemInfo.NumThreads);
  // Start the parallelization. It is done in the direction p
#pragma omp parallel shared(N_threads, global_omega_energies) firstprivate(omega_energies)
{
  
#pragma omp master
{
  N_threads = omp_get_num_threads();
  // check if each thread will get the same number of moments
  if(N1%N_threads != 0){
    std::cout << "The number of Chebyshev moments in the nonlinear optical conductivity must"
      "be a multiple of the number of threads\n" << std::flush;
    exit(1);
  }
}
#pragma omp barrier

#pragma omp for schedule(static, 1) nowait
  for(int i = 0; i < N_threads; i++){
    local_NumMoments = N1/N_threads;
    
    Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> Gamma3Aligned;
    Gamma3Aligned = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>::Zero(N2*local_NumMoments, N0);
    for(long int m = i*local_NumMoments; m < local_NumMoments*(i+1); m++)
      for(long int n = 0; n < N0; n++)
        for(long int p = 0; p < N2; p++)
          Gamma3Aligned((m-i*local_NumMoments)*N2 + p, n) = Gamma3(N0*N1*p + N0*m + n);
      
    // Delta matrix of chebyshev moments and energies
    Eigen::Matrix<std::complex<U>,-1, -1> DeltaMatrix;
    
    DeltaMatrix = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::ColMajor>::Zero(N0, N_energies);
    for(int n = 0; n < NumMoments; n++)
      for(int e = 0; e < N_energies; e++)
        DeltaMatrix(n,e) = deltaF(n, energies(e)); 

    Eigen::Matrix<std::complex<U>, -1, -1> Gamma3NNE;
    Gamma3NNE = Gamma3Aligned*DeltaMatrix;
    
    
    // Matrix of Green's functions
    Eigen::Matrix<std::complex<U>, -1, -1> GreenA, Green2A;
    Green2A = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>::Zero(N_energies, N2);
    GreenA  = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>::Zero(N_energies, local_NumMoments);

    
    
    U w1, w2;
    for(int w = 0; w < N_omegas; w++){
      w1 = frequencies2(w,0);
      w2 = frequencies2(w,1);
      
      // The scat term is the same in both cases because greenRscat and greenAscat already
      // take into account that the sign of scat is different in those cases
      for(int p = 0; p < N2; p++)
        for(int e = 0; e < N_energies; e++)
          Green2A(e, p) = greenAscat<U>(2*scat)(p, energies(e) - w1 - w2); 
      
      for(int m = 0; m < local_NumMoments; m++)
        for(int e = 0; e < N_energies; e++)
          GreenA(e, m) = greenAscat<U>(scat)(i*local_NumMoments + m, energies(e) - w1); 

      
      Eigen::Matrix<std::complex<U>, -1, -1> temp;
      temp = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(1,1);
      for(int col = 0; col < N_energies; col++){
        for(int m = 0; m < local_NumMoments; m++){
          temp = Green2A.row(col)*Gamma3NNE.block(m*N2, col, N2, 1);
          omega_energies(col, w) += temp(0,0)*GreenA(col, m);
        }
      }
    }
  }
#pragma omp critical
      global_omega_energies += omega_energies;
#pragma omp barrier
}
  return global_omega_energies;
}




template Eigen::Matrix<std::complex<float>, -1, -1> conductivity_nonlinear<float, 1u>::Gamma3shgContract_RA();
template Eigen::Matrix<std::complex<float>, -1, -1> conductivity_nonlinear<float, 2u>::Gamma3shgContract_RA();
template Eigen::Matrix<std::complex<float>, -1, -1> conductivity_nonlinear<float, 3u>::Gamma3shgContract_RA();

template Eigen::Matrix<std::complex<double>, -1, -1> conductivity_nonlinear<double, 1u>::Gamma3shgContract_RA();
template Eigen::Matrix<std::complex<double>, -1, -1> conductivity_nonlinear<double, 2u>::Gamma3shgContract_RA();
template Eigen::Matrix<std::complex<double>, -1, -1> conductivity_nonlinear<double, 3u>::Gamma3shgContract_RA();

template Eigen::Matrix<std::complex<long double>, -1, -1> conductivity_nonlinear<long double, 1u>::Gamma3shgContract_RA();
template Eigen::Matrix<std::complex<long double>, -1, -1> conductivity_nonlinear<long double, 2u>::Gamma3shgContract_RA();
template Eigen::Matrix<std::complex<long double>, -1, -1> conductivity_nonlinear<long double, 3u>::Gamma3shgContract_RA();



template Eigen::Matrix<std::complex<float>, -1, -1> conductivity_nonlinear<float, 1u>::Gamma3shgContract_RR();
template Eigen::Matrix<std::complex<float>, -1, -1> conductivity_nonlinear<float, 2u>::Gamma3shgContract_RR();
template Eigen::Matrix<std::complex<float>, -1, -1> conductivity_nonlinear<float, 3u>::Gamma3shgContract_RR();

template Eigen::Matrix<std::complex<double>, -1, -1> conductivity_nonlinear<double, 1u>::Gamma3shgContract_RR();
template Eigen::Matrix<std::complex<double>, -1, -1> conductivity_nonlinear<double, 2u>::Gamma3shgContract_RR();
template Eigen::Matrix<std::complex<double>, -1, -1> conductivity_nonlinear<double, 3u>::Gamma3shgContract_RR();

template Eigen::Matrix<std::complex<long double>, -1, -1> conductivity_nonlinear<long double, 1u>::Gamma3shgContract_RR();
template Eigen::Matrix<std::complex<long double>, -1, -1> conductivity_nonlinear<long double, 2u>::Gamma3shgContract_RR();
template Eigen::Matrix<std::complex<long double>, -1, -1> conductivity_nonlinear<long double, 3u>::Gamma3shgContract_RR();


template Eigen::Matrix<std::complex<float>, -1, -1> conductivity_nonlinear<float, 1u>::Gamma3shgContract_AA();
template Eigen::Matrix<std::complex<float>, -1, -1> conductivity_nonlinear<float, 2u>::Gamma3shgContract_AA();
template Eigen::Matrix<std::complex<float>, -1, -1> conductivity_nonlinear<float, 3u>::Gamma3shgContract_AA();

template Eigen::Matrix<std::complex<double>, -1, -1> conductivity_nonlinear<double, 1u>::Gamma3shgContract_AA();
template Eigen::Matrix<std::complex<double>, -1, -1> conductivity_nonlinear<double, 2u>::Gamma3shgContract_AA();
template Eigen::Matrix<std::complex<double>, -1, -1> conductivity_nonlinear<double, 3u>::Gamma3shgContract_AA();

template Eigen::Matrix<std::complex<long double>, -1, -1> conductivity_nonlinear<long double, 1u>::Gamma3shgContract_AA();
template Eigen::Matrix<std::complex<long double>, -1, -1> conductivity_nonlinear<long double, 2u>::Gamma3shgContract_AA();
template Eigen::Matrix<std::complex<long double>, -1, -1> conductivity_nonlinear<long double, 3u>::Gamma3shgContract_AA();
