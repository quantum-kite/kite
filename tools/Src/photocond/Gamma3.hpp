/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/


template <typename U, unsigned DIM>
Eigen::Matrix<std::complex<U>, -1, -1> conductivity_nonlinear<U, DIM>::Gamma3Contract_RRandAAblocks(){
  // Calculation of the second and third three-velocity terms. These have the 
  // product of two functions which do not depend on the frequency
  //
  int N0 = NumMoments;
  int N1 = NumMoments;
  int N2 = NumMoments;
  
  // number of threads, thread number and number of moments alocated
  // to each thread. These have to be computed inside a threaded block
  int N_threads;
  int local_NumMoments;
  

  Eigen::Matrix<std::complex<U>, -1, -1> global_omega_energies;
  Eigen::Matrix<std::complex<U>, -1, -1> omega_energies;
  global_omega_energies = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  omega_energies = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);

  // Functions that are going to be used by the contractor
  int NumMoments1 = NumMoments; U beta1 = beta; U e_fermi1 = e_fermi;
  std::function<U(int, U)> deltaF = [beta1, e_fermi1, NumMoments1](int n, U energy)->U{
    return delta(n, energy)*U(1.0/(1.0 + U(n==0)))*fermi_function(energy, e_fermi1, beta1)*kernel_jackson<U>(n, NumMoments1);
  };

  // Delta matrix of chebyshev moments and energies
  // This only needs to be calculated once. Then, it is copied to all threads. These are
  // huge matrices that are going to be copied to each thread. It has to be done in blocks
  // otherwise there will not be enough memory.

  int N_blocks = 1;
  if(const char* env_p = std::getenv("OMP_NUM_THREADS"))
    N_blocks = atoi(std::getenv("OMP_NUM_THREADS"));


  for(int block = 0; block < N_blocks; block ++){
    Eigen::Matrix<std::complex<U>,-1, -1, Eigen::RowMajor> DeltaGreenRMatrix;
    Eigen::Matrix<std::complex<U>,-1, -1, Eigen::RowMajor> DeltaGreenAMatrix;
    DeltaGreenRMatrix = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>::Zero(N_energies, N0*N2/N_blocks);
    DeltaGreenAMatrix = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>::Zero(N_energies, N0*N2/N_blocks);
    for(int n = 0; n < N0; n++)
      for(int p = 0; p < N2/N_blocks; p++)
        for(int e = 0; e < N_energies; e++){
          DeltaGreenRMatrix(e, p*N0 + n) = deltaF(p + block*N2/N_blocks, energies(e))*greenRscat<U>(2*scat)(n, energies(e)); 
          DeltaGreenAMatrix(e, p*N0 + n) = deltaF(n, energies(e))*greenAscat<U>(2*scat)(p + block*N2/N_blocks, energies(e)); 
        }

    omp_set_num_threads(systemInfo.NumThreads);
#pragma omp parallel shared(N_threads, global_omega_energies) firstprivate(DeltaGreenAMatrix, DeltaGreenRMatrix, omega_energies)
  {

    // block to determine the number of threads and check if it's a valid number
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


      // Parallelization is done in the m direction this time
#pragma omp for schedule(static, 1) nowait
    for(int i = 0; i < N_threads; i++){
      local_NumMoments = N1/N_threads;

      // The operations that will be done next on the Gamma3 matrix are very time-consuming
      // and are of order NumMoments^3 * N_energies, so we want to reduce the computation time
      // as much as possible. For that reason, this block aligns the matrix entries in thememory so
      // that the matrix multiplication that follows is as fast as possible. This operation is
      // of order NumMoments^3
      Eigen::Matrix<std::complex<U>, -1, -1, Eigen::ColMajor> Gamma3Aligned;
      Gamma3Aligned = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::ColMajor>::Zero(N0*N2/N_blocks, local_NumMoments);
      for(long int p = 0; p < N2/N_blocks; p++)
        for(long int m = i*local_NumMoments; m < (i+1)*local_NumMoments; m++)
          for(long int n = 0; n < N0; n++)
            Gamma3Aligned(p*N0 + n, m - i*local_NumMoments) = Gamma3(N0*N1*(p + block*N2/N_blocks) + N0*m + n);
        
      
      // Perform the matrix product of Gamma3Aligned with the matrix of 
      // Green functions and Dirac Deltas
      Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> Gamma3NER;
      Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> Gamma3NEA;

      Gamma3NER = DeltaGreenRMatrix*Gamma3Aligned;
      Gamma3NEA = DeltaGreenAMatrix*Gamma3Aligned;
      
      // Matrix of Green's functions
      Eigen::Matrix<std::complex<U>, -1, -1> GreenR, GreenA;
      GreenR = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>::Zero(local_NumMoments, N_energies);
      GreenA = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>::Zero(local_NumMoments, N_energies);
      
      for(int w = 0; w < N_omegas; w++){
        // The scat term is the same in both cases because greenRscat and greenAscat already
        // take into account that the sign of scat is different in those cases
        for(int m = 0; m < local_NumMoments; m++)
          for(int e = 0; e < N_energies; e++){
            GreenR(m, e) = greenRscat<U>(scat)(m + i*local_NumMoments, energies(e) - frequencies(w)); 
            GreenA(m, e) = greenAscat<U>(scat)(m + i*local_NumMoments, energies(e) - frequencies(w)); 
          }
        
        Eigen::Matrix<std::complex<U>, -1, -1> temp;
        for(int e = 0; e < N_energies; e++){
          temp  = Gamma3NER.row(e)*GreenR.col(e); 
          temp += Gamma3NEA.row(e)*GreenA.col(e); 
          omega_energies(e, w) = temp(0,0);
        }
      }
    }
#pragma omp critical
        global_omega_energies += omega_energies;
#pragma omp barrier
  }
  }
  return global_omega_energies;
}

template <typename U, unsigned DIM>
Eigen::Matrix<std::complex<U>, -1, -1> conductivity_nonlinear<U, DIM>::Gamma3Contract_RA(){
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
  
  // Start the parallelization. It is done in the direction p
  omp_set_num_threads(systemInfo.NumThreads);
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
    GreenR = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>::Zero(N_energies, N0);
    GreenA = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>::Zero(N_energies, local_NumMoments);



    
    for(int w = 0; w < N_omegas; w++){
      
      // The scat term is the same in both cases because greenRscat and greenAscat already
      // take into account that the sign of scat is different in those cases
      
      for(int n = 0; n < N0; n++)
        for(int e = 0; e < N_energies; e++)
          GreenR(e, n) = greenRscat<U>(scat)(n, energies(e) + frequencies(w)); 

      for(int p = 0; p < local_NumMoments; p++)
        for(int e = 0; e < N_energies; e++)
          GreenA(e, p) = greenAscat<U>(scat)(i*local_NumMoments + p, energies(e) + frequencies(w)); 

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
{
      
      global_omega_energies += omega_energies;
}
#pragma omp barrier
}
  return global_omega_energies;
}
