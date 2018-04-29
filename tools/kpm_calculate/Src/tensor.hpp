#include "headers.hpp"
#include "value.hpp"

int int_pow(int base, int exp){
    int result = 1;
    while (exp)
    {
        if (exp & 1)
           result *= base;
        exp /= 2;
        base *= base;
    }
    return result;
}

using namespace std::placeholders;  // for _1, _2, _3...

double kernel_jackson(int n, int M){
	double f = M_PI/(M+1);
	return ((M+1-n)*cos(n*f) + sin(n*f)/tan(f))/(M+1);
}

double lorentz(double lambda, int n, int M){
	return sinh(lambda*(1.0 - (double)n/M))/sinh(lambda);
}
std::function<double(int, int)> kernel_lorentz(double lambda){
	return std::bind(lorentz, lambda, _1, _2);
}

double kernel_dirichlet(int n, int M){
	return (double)(n < M);
}

template <typename TC>
TC green(int n, int sigma, TC energy){
	const TC i(0.0,1.0); 
	TC sq = sqrt(1.0 - energy*energy);
	return 2.0*sigma/sq*i*exp(-sigma*n*1.0*acos(energy)*i);
}

template <typename T>
std::complex<T> greenR(int n, T energy, T scat){
  return green(n,  1, std::complex<T>(energy,  scat))/(1.0 + T(n==0));
}

template <typename T>
std::complex<T> greenA(int n, T energy, T scat){
  return green(n, -1, std::complex<T>(energy, -scat))/(1.0 + T(n==0));
}

template <typename T>
std::function<std::complex<T>(int, T)> greenRscat(T scat){
  return std::bind(greenR<T>, _1, _2, scat);
}

template <typename T>
std::function<std::complex<T>(int, T)> greenAscat(T scat){
  return std::bind(greenA<T>, _1, _2, scat);
}


template <typename TR>
TR delta(int n, TR energy){
	TR sq = sqrt(1.0 - energy*energy);
	if(energy < 1 and energy > -1)
		return 2.0/M_PI/sq*cos(n*acos(energy));
	else
		return 0;
}



