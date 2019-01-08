#include <Eigen/Dense>
#include <iostream>
#include <complex>
#include <string>
#include "macros.hpp"

template <typename T>	
std::complex<T> integrate(Eigen::Matrix<T, -1, 1> energies, Eigen::Matrix<std::complex<T>, -1, 1> integrand);


template <typename T>
T fermi_function(T energy, T mu, T beta);


std::string num2str3(int dir_num);
std::string num2str2(int dir_num);
template <typename T>
std::complex<T> contract1(
    std::function<T(int, T)> f0, int N_moments, 
    const Eigen::Array<std::complex<T>, -1, -1>& Gamma, 
    const Eigen::Matrix<T, -1, 1>& energies);

int int_pow(int base, int exp);

template <typename T>
T kernel_jackson(int n, int M);



double lorentz(double lambda, int n, int M);
std::function<double(int, int)> kernel_lorentz(double lambda);
double kernel_dirichlet(int n, int M);

template <typename TC>
TC green(int n, int sigma, TC energy);

template <typename T>
std::complex<T> dgreen(int n, int sigma, std::complex<T> energy);

template <typename T>
std::complex<T> greenR(int n, T energy, T scat);

template <typename T>
std::complex<T> greenA(int n, T energy, T scat);

template <typename T>
std::function<std::complex<T>(int, T)> greenRscat(T scat);

template <typename T>
std::function<std::complex<T>(int, T)> greenAscat(T scat);

template <typename TR>
TR delta(int n, TR energy);
