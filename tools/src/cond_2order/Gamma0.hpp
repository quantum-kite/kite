/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/


template <typename U, unsigned DIM>
Eigen::Matrix<std::complex<U>, -1, -1> conductivity_nonlinear<U, DIM>::Gamma0contract(){

  Eigen::Matrix<std::complex<U>, -1, -1> omega_energies;
  Eigen::Matrix<std::complex<U>, -1, -1> temp;

  omega_energies = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  temp = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(1,1);

  return omega_energies;
}
