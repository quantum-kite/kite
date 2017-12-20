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
	return pow(-1,n)*2.0*sigma/sq*i*exp(-sigma*n*1.0*acos(energy)*i);
}

template <typename TR>
TR delta(int n, TR energy){
	TR sq = sqrt(1.0 - energy*energy);
	if(energy < 1 and energy > -1)
		return 2.0/M_PI/sq*cos(n*acos(energy));
	else
		return 0;
}



template <int ENER, int CHEB, typename T>
class tensor {
	
	private:
		int initialized;
		Eigen::Array<T, -1, 1> energies;
	
	public:
	
		Eigen::Array<std::complex<T>, -1, -1> data;
		
		int N_energies;
		int N_moments;
		
		tensor(Eigen::Array<T, 1, -1>, int);
		void set_delta(std::function<double(int, int)> kernel);
		void set_gamma(Eigen::Array<std::complex<T>,-1,-1>);
		
		void set_green(std::function<double(int, int)>, std::string);
		void set_green(double, std::string);
		void set_green_imag(double, std::string);
		
		
		
		tensor<ENER+1, CHEB-1, T> operator* (tensor<1, 1, T>);
		
		
		// Overloading the access operator
		template <int CHEB1 = CHEB, int ENER1 = ENER>
		typename std::enable_if<CHEB1 == 0 and ENER1 == 1, std::complex<T>>::type operator()(int E1){
			return data(E1);
		}
		
		template <int CHEB1 = CHEB, int ENER1 = ENER>
		typename std::enable_if<CHEB1 == 0 and ENER1 == 2, std::complex<T>>::type operator()(int E1, int E2){
			return data(E1*N_energies + E2);
		}
		
		template <int CHEB1 = CHEB, int ENER1 = ENER>
		typename std::enable_if<CHEB1 == 0 and ENER1 == 3, std::complex<T>>::type operator()(int E1, int E2, int E3){
			return data((E1*N_energies + E2)*N_energies + E3);
		}
		
		template <int CHEB1 = CHEB, int ENER1 = ENER>
		typename std::enable_if<CHEB1 == 1 and ENER1 == 1, std::complex<T>>::type operator()(int n1, int n2){
			return data(n1, n2);
		}

		// Overloading the access operator for "value" types
		template <int CHEB1 = CHEB, int ENER1 = ENER>
		typename std::enable_if<CHEB1 == 0 and ENER1 == 1, std::complex<T>>::type operator()(value<T> E1){
			return data(E1.ind);
		}
		
		template <int CHEB1 = CHEB, int ENER1 = ENER>
		typename std::enable_if<CHEB1 == 0 and ENER1 == 2, std::complex<T>>::type operator()(value<T> E1, value<T> E2){
			return data(E1.ind*N_energies + E2.ind);
		}
		
		template <int CHEB1 = CHEB, int ENER1 = ENER>
		typename std::enable_if<CHEB1 == 0 and ENER1 == 3, std::complex<T>>::type operator()(value<T> E1, value<T> E2, value<T> E3){
			return data((E1.ind*N_energies + E2.ind)*N_energies + E3.ind);
		}
		
};









template <int ENER, int CHEB, typename T>
tensor<ENER, CHEB, T>::tensor(Eigen::Array<T, 1, -1> input_energies, int moments) {
	N_energies = input_energies.cols();
	N_moments  = moments;
	energies   = input_energies;
}


template <int ENER, int CHEB, typename T>
void tensor<ENER, CHEB, T>::set_green(std::function<double(int, int)> kernel, std::string adv_ret){
	/* The 'data' variable is still unknown. In this function, we decide what it'll be: a Green's function.
	 * It therefore has only two entries: (1) the Chebyshev moment and (2) the Energy value.
	 * This function declares the 'data' variable to be a matrix of 'N_energies' rows and 'N_moments' columns. Then the 
	 * function fills the matrix with the values calculated from the Green's function.
	 * This function is overloaded, so it can either take a kernel parameter or, instead, use
	 * a finite GAMMA as a broadening parameter
	 * */
	 
	data = Eigen::Array<std::complex<T>, -1, -1>::Zero(N_energies, N_moments);				// declare the array
	initialized	= 1; 																																			// set the 'initialized' flag to true
	
	
	// Verify whether it's an advanced or retarded Green's function
	int sigma;
	if(adv_ret == "retarded")
		sigma = -1;
	else{
		if(adv_ret == "advanced")
			sigma = 1;
		else{
			std::cout << "Unrecognized parameter for the type of Green's function. Use 'advanced' or 'retarded'. Exiting program.\n";
			exit(0);
		}
	}
		
	
	double energy;	
	for(int i = 0; i < N_energies; i++){
		energy = energies(i);
		for(int n = 0; n < N_moments; n++)
			data(i,n) = kernel(n, N_moments)/(1.0 + (double)(n==0))*green(n, sigma , energy);
	}
}

template <int ENER, int CHEB, typename T>
void tensor<ENER, CHEB, T>::set_green(double GAMMA, std::string adv_ret){
	/* The 'data' variable is still unknown. In this function, we decide what it'll be: a Green's function.
	 * It therefore has only two entries: (1) the Chebyshev moment and (2) the Energy value.
	 * This function declares the 'data' variable to be a matrix of 'N_energies' rows and 'N_moments' columns. Then the 
	 * function fills the matrix with the values calculated from the Green's function.
	 * This function is overloaded, so it can either take a kernel parameter or, instead, use
	 * a finite GAMMA as a broadening parameter
	 * */
	
	data = Eigen::Array<std::complex<T>, -1, -1>::Zero(N_energies, N_moments);				// declare the array
	initialized	= 1; 																																			// set the 'initialized' flag to true
	
	
	// Verify whether it's an advanced or retarded Green's function
	int sigma;
	if(adv_ret == "retarded")
		sigma = -1;
	else{
		if(adv_ret == "advanced")
			sigma = 1;
		else{
			std::cout << "Unrecognized parameter for the type of Green's function. Use 'advanced' or 'retarded'. Exiting program.\n";
			exit(0);
		}
	}
		

	std::complex<double> complex_energy;	
	for(int i = 0; i < N_energies; i++){
		complex_energy = std::complex<double> (energies(i), GAMMA);
		for(int n = 0; n < N_moments; n++)
			data(i,n) = 1.0/(1.0 + (double)(n==0))*green(n, sigma, complex_energy);
	}
}

template <int ENER, int CHEB, typename T>
void tensor<ENER, CHEB, T>::set_green_imag(double GAMMA, std::string adv_ret){
	/* The 'data' variable is still unknown. In this function, we decide what it'll be: a Green's function.
	 * It therefore has only two entries: (1) the Chebyshev moment and (2) the Energy value.
	 * This function declares the 'data' variable to be a matrix of 'N_energies' rows and 'N_moments' columns. Then the 
	 * function fills the matrix with the values calculated from the Green's function.
	 * This function is overloaded, so it can either take a kernel parameter or, instead, use
	 * a finite GAMMA as a broadening parameter
	 * */
	data = Eigen::Array<std::complex<T>, -1, -1>::Zero(N_energies, N_moments);				// declare the array
	initialized	= 1; 																																			// set the 'initialized' flag to true
	
	
	// Verify whether it's an advanced or retarded Green's function
	int sigma;
	if(adv_ret == "retarded")
		sigma = -1;
	else{
		if(adv_ret == "advanced")
			sigma = 1;
		else{
			std::cout << "Unrecognized parameter for the type of Green's function. Use 'advanced' or 'retarded'. Exiting program.\n";
			exit(0);
		}
	}
		

	std::complex<double> complex_energy;	
	for(int i = 0; i < N_energies; i++){
		complex_energy = std::complex<double> (energies(i), GAMMA);
		for(int n = 0; n < N_moments; n++)
			data(i,n) = 1.0/(1.0 + (double)(n==0))*green(n, sigma, complex_energy).imag();
	}
}


template <int ENER, int CHEB, typename T>
void tensor<ENER, CHEB, T>::set_delta(std::function<double(int, int)> kernel){
	/* The 'data' variable is still unknown. In this function, we decide what it'll be: a Delta function.
	 * It therefore has only two entries: (1) the Chebyshev moment and (2) the Energy value.
	 * This function declares the 'data' variable to be a matrix of 'N_energies' rows and 'N_moments' columns. Then the 
	 * function fills the matrix with the values calculated from the Delta function.
	 * */
	
	data = Eigen::Array<std::complex<T>, -1, -1>::Zero(N_energies, N_moments);		// declare the array
	initialized	= 1; 																// set the 'initialized' flag to true

	double energy;	
	for(int i = 0; i < N_energies; i++){
		energy = energies(i);
		if(energy > -0.9999999 and energy < 0.9999999){
			for(int n = 0; n < N_moments; n++){
				data(i,n) = kernel(n, N_moments)/(1.0 + (double)(n==0))*delta(n, energy);
			}
		}
	}
	
	//std::cout << "Successfully exited delta\n";fflush(stdout);
}


template <int ENER, int CHEB, typename T>
void tensor<ENER, CHEB, T>::set_gamma(Eigen::Array<std::complex<T>,-1,-1> input){
	data = input;
};







template <int ENER, int CHEB, typename T>
tensor<ENER+1, CHEB-1, T> tensor<ENER, CHEB, T>::operator* (tensor<1, 1, T> other_matrix){
	tensor<ENER+1, CHEB-1, T> temp(energies, N_moments);
	int energy = pow(N_energies, ENER);
	int moments  = pow(N_moments , CHEB);
	
	temp.data = Eigen::Array<std::complex<T>,-1,-1>::Zero(energy*N_energies, moments/N_moments);	
	
	for(int i = 0; i < energy; i++){
		for(int j = 0; j < N_energies; j++){
			for(int n = 0; n < moments/N_moments; n++){
				for(int m = 0; m < N_moments; m++){
					temp.data(i*N_energies + j, n) += data(i, n*N_moments + m)*other_matrix(j, m);
				}
			}
		}
	}
	return temp;
}
