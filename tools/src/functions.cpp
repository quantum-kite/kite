#include "functions.hpp"

template <typename T>	
std::complex<T> integrate(Eigen::Matrix<T, -1, 1> energies, Eigen::Matrix<std::complex<T>, -1, 1> integrand){
	if(energies.rows() != integrand.rows() or energies.cols() != integrand.cols()){
		std::cout << "x and y arrays in the integrator must have the same number of elements. Exiting.\n";
		exit(1);
	}
	
	int N = energies.cols()*energies.rows();
	std::complex<T> sum(0,0);
	
	for(int i = 0; i < N - 1; i++){
		sum += (energies(i) - energies(i+1))*(integrand(i) + integrand(i+1))/T(2.0);
	}
	
	return sum;
}

// Instantiations
template std::complex<float> integrate<float>(Eigen::Matrix<float, -1, 1>, Eigen::Matrix<std::complex<float>, -1, 1>);
template std::complex<double> integrate<double>(Eigen::Matrix<double, -1, 1>, Eigen::Matrix<std::complex<double>, -1, 1>);
template std::complex<long double> integrate<long double>(Eigen::Matrix<long double, -1, 1>, Eigen::Matrix<std::complex<long double>, -1, 1>);

template <typename T>
T fermi_function(T energy, T mu, T beta){
	return 1.0/(1.0 + exp(beta*(energy - mu)));
}
template float fermi_function(float, float, float);
template double fermi_function(double, double, double);
template long double fermi_function(long double, long double, long double);

std::string num2str3(int dir_num){
  std::string dir;
 
  switch(dir_num){
    case 0:
      dir = "xxx"; break;
    case 1:
      dir = "xxy"; break;
    case 2:
      dir = "xxz"; break;
    case 3:
      dir = "xyx"; break;
    case 4:
      dir = "xyy"; break;
    case 5:
      dir = "xyz"; break;
    case 6:
      dir = "xzx"; break;
    case 7:
      dir = "xzy"; break;
    case 8:
      dir = "xzz"; break;
    case 9:
      dir = "yxx"; break;
    case 10:
      dir = "yxy"; break;
    case 11:
      dir = "yxz"; break;
    case 12:
      dir = "yyx"; break;
    case 13:
      dir = "yyy"; break;
    case 14:
      dir = "yyz"; break;
    case 15:
      dir = "yzx"; break;
    case 16:
      dir = "yzy"; break;
    case 17:
      dir = "yzz"; break;
    case 18:
      dir = "zxx"; break;
    case 19:
      dir = "zxy"; break;
    case 20:
      dir = "zxz"; break;
    case 21:
      dir = "zyx"; break;
    case 22:
      dir = "zyy"; break;
    case 23:
      dir = "zyz"; break;
    case 24:
      dir = "zzx"; break;
    case 25:
      dir = "zzy"; break;
    case 26:
      dir = "zzz"; break;
    default:
      std::cout << "Invalid direction in num2str_dir3.\n"; exit(1);
  }
  return dir;
}

std::string num2str2(int dir_num){
  std::string dir;
 
  switch(dir_num){
    case 0:
      dir = "xx"; break;
    case 1:
      dir = "yy"; break;
    case 2:
      dir = "zz"; break;
    case 3:
      dir = "xy"; break;
    case 4:
      dir = "xz"; break;
    case 5:
      dir = "yx"; break;
    case 6:
      dir = "yz"; break;
    case 7:
      dir = "zx"; break;
    case 8:
      dir = "zy"; break;
    default:
      std::cout << "Invalid direction for the optical conductivity.\n"; exit(1);
  }
  return dir;
}


template <typename T>
std::complex<T> contract1(
    std::function<T(int, T)> f0, int N_moments, 
    const Eigen::Array<std::complex<T>, -1, -1>& Gamma, 
    const Eigen::Matrix<T, -1, 1>& energies){
    debug_message("Entered contract1\n");

    int N_energies = energies.rows();

    T energy;
    Eigen::Matrix<std::complex<T>, -1, 1> GammaE;
    GammaE = Eigen::Matrix<std::complex<T>, -1, 1>::Zero(N_energies, 1);
    
    for(int e = 0; e < N_energies; e++){
      energy = energies(e);
      for(int m = 0; m < N_moments; m++){
        GammaE(e) += Gamma(m)*f0(m, energy);
      }
    }
    debug_message("Left contract1");
    return integrate(energies, GammaE);
}

template std::complex<float> contract1( std::function<float(int, float)>, int, const Eigen::Array<std::complex<float>, -1, -1>&, const Eigen::Matrix<float, -1, 1>&);
template std::complex<double> contract1( std::function<double(int, double)>, int, const Eigen::Array<std::complex<double>, -1, -1>&, const Eigen::Matrix<double, -1, 1>&);
template std::complex<long double> contract1( std::function<long double(int, long double)>, int, const Eigen::Array<std::complex<long double>, -1, -1>&, const Eigen::Matrix<long double, -1, 1>&);


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

template <typename T>
T kernel_jackson(int n, int M){
	T f = M_PI/(M+1);
	return ((M+1-n)*cos(n*f) + sin(n*f)/tan(f))/(M+1);
}

// Instantiations
template float kernel_jackson(int, int);
template double kernel_jackson(int, int);
template long double kernel_jackson(int, int);







using namespace std::placeholders;  // for _1, _2, _3...

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
	TC sq = sqrt(TC(1.0) - energy*energy);
	return -TC(2.0*sigma)/sq*i*exp(-TC(sigma*n)*acos(energy)*i);
}

template std::complex<float> green(int, int, std::complex<float>);
template std::complex<double> green(int, int, std::complex<double>);
template std::complex<long double> green(int, int, std::complex<long double>);




template <typename T>
std::complex<T> dgreen(int n, int sigma, std::complex<T> energy){
	const std::complex<T> i(0.0,1.0); 
  
  std::complex<T> den = T(1.0) - energy*energy;
  std::complex<T>  sq = sqrt(den);
	return -T(2.0*sigma)/den*i*exp(-T(sigma*n)*acos(energy)*i)*(T(n*sigma)*i + energy/sq);
}
template std::complex<float> dgreen(int, int, std::complex<float>);
template std::complex<double> dgreen(int, int, std::complex<double>);
template std::complex<long double> dgreen(int, int, std::complex<long double>);



template <typename T>
std::complex<T> greenR(int n, T energy, T scat){
  return green(n,  1, std::complex<T>(energy,  scat))*T(1.0/(1.0 + T(n==0)));
}

template <typename T>
std::complex<T> greenA(int n, T energy, T scat){
  return green(n, -1, std::complex<T>(energy, -scat))*T(1.0/(1.0 + T(n==0)));
}

template std::complex<float> greenR(int, float, float);
template std::complex<double> greenR(int, double, double);
template std::complex<long double> greenR(int, long double, long double);
template std::complex<float> greenA(int, float, float);
template std::complex<double> greenA(int, double, double);
template std::complex<long double> greenA(int, long double, long double);




template <typename T>
std::function<std::complex<T>(int, T)> greenRscat(T scat){
  return std::bind(greenR<T>, _1, _2, scat);
}

template <typename T>
std::function<std::complex<T>(int, T)> greenAscat(T scat){
  return std::bind(greenA<T>, _1, _2, scat);
}

template std::function<std::complex<float>(int, float)> greenRscat(float);
template std::function<std::complex<double>(int, double)> greenRscat(double);
template std::function<std::complex<long double>(int, long double)> greenRscat(long double);
template std::function<std::complex<float>(int, float)> greenAscat(float);
template std::function<std::complex<double>(int, double)> greenAscat(double);
template std::function<std::complex<long double>(int, long double)> greenAscat(long double);


template <typename TR>
TR delta(int n, TR energy){
	TR sq = sqrt(1.0 - energy*energy);
	if(energy < 1 and energy > -1)
		return 2.0/M_PI/sq*cos(n*acos(energy));
	else
		return 0;
}

template float delta(int, float);
template double delta(int, double);
template long double delta(int, long double);



