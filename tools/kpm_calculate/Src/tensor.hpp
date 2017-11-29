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

std::complex<double> green(int n, int sigma, std::complex<double> energy){
	const std::complex<double> i(0.0,1.0); 
	std::complex<double> sq = sqrt(1.0 - energy*energy);
	return pow(-1,n)*2.0*sigma/sq*i*exp(-sigma*n*1.0*acos(energy)*i);
}

double delta(int n, double energy){
	double sq = sqrt(1.0 - energy*energy);
	if(energy < 1 and energy > -1)
		return 2.0/M_PI/sq*cos(n*acos(energy));
	else
		return 0;
}
/*
std::complex<double> delta(int n, std::complex<double> energy){
	std::complex<double> sq = sqrt(1.0 - energy*energy);
	return 2.0/M_PI/sq*cos(n*1.0*acos(energy));
}
*/






template <int ENER, int CHEB>
class tensor {
	
	
	int initialized;
	Eigen::Array<double, 1, -1> energies;
	
	public:
		Eigen::Matrix<std::complex<double>, -1, -1> data;
		int N_energies;
		int N_moments;
	
		tensor(Eigen::Array<double, 1, -1>, int);
		void set_green(std::function<double(int, int)>, std::string);
		void set_green(double, std::string);
		void set_green_imag(double, std::string);
		void set_delta(std::function<double(int, int)> kernel);
		void set_gamma(Eigen::Array<std::complex<double>,-1,-1>);
		tensor<ENER+1, CHEB-1>  operator* (const tensor<1, 1> );
		
		std::complex<double> operator()(int);
		std::complex<double> operator()(int, int);
		std::complex<double> operator()(int, int, int);
		
		
		std::complex<double> operator()(value<double>);
		std::complex<double> operator()(value<double>, value<double>);
		std::complex<double> operator()(value<double>, value<double>, value<double>);
};








template <int ENER, int CHEB>
tensor<ENER, CHEB>::tensor(Eigen::Array<double, 1, -1> input_energies, int moments) {
	N_energies = input_energies.cols();
	N_moments  = moments;
	energies   = input_energies;
}

template <int ENER, int CHEB>
void tensor<ENER, CHEB>::set_green(std::function<double(int, int)> kernel, std::string adv_ret){
	/* The 'data' variable is still unknown. In this function, we decide what it'll be: a Green's function.
	 * It therefore has only two entries: (1) the Chebyshev moment and (2) the Energy value.
	 * This function declares the 'data' variable to be a matrix of 'N_energies' rows and 'N_moments' columns. Then the 
	 * function fills the matrix with the values calculated from the Green's function.
	 * This function is overloaded, so it can either take a kernel parameter or, instead, use
	 * a finite GAMMA as a broadening parameter
	 * */
	
	data = Eigen::Array<std::complex<double>, -1, -1>::Zero(N_energies, N_moments);				// declare the array
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

template <int ENER, int CHEB>
void tensor<ENER, CHEB>::set_green(double GAMMA, std::string adv_ret){
	/* The 'data' variable is still unknown. In this function, we decide what it'll be: a Green's function.
	 * It therefore has only two entries: (1) the Chebyshev moment and (2) the Energy value.
	 * This function declares the 'data' variable to be a matrix of 'N_energies' rows and 'N_moments' columns. Then the 
	 * function fills the matrix with the values calculated from the Green's function.
	 * This function is overloaded, so it can either take a kernel parameter or, instead, use
	 * a finite GAMMA as a broadening parameter
	 * */
	
	
	data = Eigen::Array<std::complex<double>, -1, -1>::Zero(N_energies, N_moments);				// declare the array
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

template <int ENER, int CHEB>
void tensor<ENER, CHEB>::set_green_imag(double GAMMA, std::string adv_ret){
	/* The 'data' variable is still unknown. In this function, we decide what it'll be: a Green's function.
	 * It therefore has only two entries: (1) the Chebyshev moment and (2) the Energy value.
	 * This function declares the 'data' variable to be a matrix of 'N_energies' rows and 'N_moments' columns. Then the 
	 * function fills the matrix with the values calculated from the Green's function.
	 * This function is overloaded, so it can either take a kernel parameter or, instead, use
	 * a finite GAMMA as a broadening parameter
	 * */
	
	
	data = Eigen::Array<std::complex<double>, -1, -1>::Zero(N_energies, N_moments);				// declare the array
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

template <int ENER, int CHEB>
void tensor<ENER, CHEB>::set_delta(std::function<double(int, int)> kernel){
	/* The 'data' variable is still unknown. In this function, we decide what it'll be: a Delta function.
	 * It therefore has only two entries: (1) the Chebyshev moment and (2) the Energy value.
	 * This function declares the 'data' variable to be a matrix of 'N_energies' rows and 'N_moments' columns. Then the 
	 * function fills the matrix with the values calculated from the Delta function.
	 * */
	
	//std::cout << "entered delta\n";fflush(stdout);
	data = Eigen::Array<std::complex<double>, -1, -1>::Zero(N_energies, N_moments);				// declare the array
	initialized	= 1; 																																			// set the 'initialized' flag to true

	//std::cout << "energies inside the delta function \n";fflush(stdout);
	double energy;	
	for(int i = 0; i < N_energies; i++){
		energy = energies(i);
		//std::cout << "i: " << energy << "\n\n"; fflush(stdout);
		if(energy > -0.9999999 && energy < 0.9999999){
			for(int n = 0; n < N_moments; n++){
				//std::cout << "n: " << n << " " << kernel(n, N_moments)/(1.0 + (double)(n==0))*delta(n, energy) << " ";
				data(i,n) = kernel(n, N_moments)/(1.0 + (double)(n==0))*delta(n, energy);
				//std::cout << "n: " << n << " " << data(i,n) << " ";
			}
		}
	}
	
	//std::cout << "Successfully exited delta\n";fflush(stdout);
}


template <int ENER, int CHEB>
void tensor<ENER, CHEB>::set_gamma(Eigen::Array<std::complex<double>,-1,-1> input){
	data = input;
	//std::cout << "GAMMA\n" << input << "\n";fflush(stdout);
};

/* Overloading the access operator for integers
 * */

template <>
std::complex<double> tensor<1, 0>::operator()(int E1){
	return data(E1);
}

template <>
std::complex<double> tensor<2, 0>::operator()(int E1, int E2){
	return data(E1*N_energies + E2);
}

template <>
std::complex<double> tensor<3,0>::operator()(int E1, int E2, int E3){
	
	return data((E1*N_energies + E2)*N_energies + E3);
}

template <>
std::complex<double> tensor<1, 1>::operator()(int n1, int n2){
	return data(n1, n2);
}
/* Overloading the access operator for "value" types
 * */


template <>
std::complex<double> tensor<1, 0>::operator()(value<double> E1){
	return data(E1.ind);
}

template <>
std::complex<double> tensor<2, 0>::operator()(value<double> E1, value<double> E2){
	return data(E1.ind*N_energies + E2.ind);
}

template <>
std::complex<double> tensor<3,0>::operator()(value<double> E1, value<double> E2, value<double> E3){
	return data((E1.ind*N_energies + E2.ind)*N_energies + E3.ind);
}


/*
template <>
std::complex<double> tensor<1, 0>::operator()(double E1){
	value<double> e(energies);
	return data((int)E1);
}*/
/*
template <> 
std::complex<double> tensor<2, 0>::operator()(int n1, int n2){
	return data(n1, n2);
}*/



template <int ENER, int CHEB>
tensor<ENER+1, CHEB-1> tensor<ENER, CHEB>::operator* (tensor<1, 1> other_matrix){
	tensor<ENER+1, CHEB-1> temp(energies, N_moments);
	int energy = pow(N_energies, ENER);
	int moments  = pow(N_moments , CHEB);
	
	temp.data = Eigen::Array<std::complex<double>,-1,-1>::Zero(energy*N_energies, moments/N_moments);	
	/*
	std::cout << "ENER: " << ENER << " CHEB: " << CHEB << "\n";fflush(stdout);
	std::cout << "energy: " << energy << " moments: " << moments << "\n";fflush(stdout);
	
	
	std::cout << "division: " <<  moments/N_moments << "\n";fflush(stdout);
	
	std::cout << "N_energies: " << N_energies << " N_moments: " << N_moments << "\n";
	
	std::cout << data.rows() << "x" << data.cols() << "\n";
	std::cout << temp.data.rows() << "x" << temp.data.cols() << "\n";
	std::cout << other_matrix.data.rows() << "x" << other_matrix.data.cols() << "\n";
	
	std::cout << "#####\n";fflush(stdout);
	
	
	std::cout << "temp.data\n";fflush(stdout);
	
	int NNN1 = 0;
	int NNN2 = 10;
	std::cout << "temp2: \n";fflush(stdout);
	for(int i = NNN1; i < NNN2; i++){
		for(int j = NNN1; j < NNN2; j++){
			std::cout << temp.data(i, j) << "\t";
		}
		std::cout << "\n";
	}
	
	
	std::cout << "does this exist for mu? ";fflush(stdout);
	std::cout << data(0,0) << "\n";fflush(stdout);
	*/
	for(int i = 0; i < energy; i++){
		//std::cout << "i" << i << " ";fflush(stdout);
		for(int j = 0; j < N_energies; j++){
			//std::cout << "j" << j << " ";fflush(stdout);
			for(int n = 0; n < moments/N_moments; n++){
				//std::cout << "n" << n << " ";fflush(stdout);
				//std::cout << "(i*N_energies + j)" << i*N_energies + j << " ";fflush(stdout);
				temp.data(i*N_energies + j, n) = 0;
				for(int m = 0; m < N_moments; m++){
					//std::cout << "m" << m << " ";fflush(stdout);
					//other_matrix(j, m);
					//std::cout << "...";fflush(stdout)
					//std::cout << "why tho: " <<  data(i, n*N_moments + m)*other_matrix(j, m) << " ";fflush(stdout);
					temp.data(i*N_energies + j, n) += data(i, n*N_moments + m)*other_matrix(j, m);
				}
			}
		}
	}
	/*
	std::cout << "temp2: \n";fflush(stdout);
	for(int i = NNN1; i < NNN2; i++){
		for(int j = NNN1; j < NNN2; j++){
			std::cout << temp.data(i, j) << "\t";
		}
		std::cout << "\n";
	}*/
	
	return temp;
}
