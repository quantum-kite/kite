#include "headers.hpp"


template <typename T>
T fermi_function(T energy, T mu, T beta){
	return 1.0/(1.0 + exp(beta*(energy - mu)));
}

template <typename T>	
T fermi_function(value<T> E1, T mu,T beta){
	return 1.0/(1.0 + exp(beta*(E1.val - mu)));
}

template <typename T>	
std::complex<T> integrate(Eigen::Array<T, -1, -1> energies, Eigen::Array<std::complex<T>, -1, -1> integrand){
	if(energies.rows() != integrand.rows() or energies.cols() != integrand.cols()){
		std::cout << "x and y arrays in the integrator must have the same number of elements. Exiting.\n";
		exit(0);
	}
	
	int N = energies.cols()*energies.rows();
	std::complex<T> sum(0,0);
	
	for(int i = 0; i < N - 1; i++){
		sum += (energies(i) - energies(i+1))*(integrand(i) + integrand(i+1))/T(2.0);
	}
	
	return sum;
}

template <typename T, unsigned DIM>
void calc_dos(info<T,DIM> *config, double lim, int N_energies){
	/* Calculates the density of states for N_energies energies in the range ]-lim, lim[. 
	 * These energies are in the KPM scale. 
	 */
	 
	// Array of energies in which the density of states will be calculated
	Eigen::Matrix<T, Eigen::Dynamic, 1> energies;
	energies = Eigen::Matrix<T, Eigen::Dynamic, 1>::LinSpaced(N_energies, -lim, lim); 
	int n_moments = config->num_moments(0);
	
	
	tensor<1,1,T> D(energies, n_moments);		// Dirac delta
	tensor<0,1,T> mu(energies, n_moments);		// Matrix calculated from the kpm_transport program. Only depends on the Chebyshev moments
	tensor<1,0,T> DOS(energies, n_moments);		// Density of states. Only depends on the energy
	
	// Set the values and perform the calculation
	D.set_delta(kernel_jackson);
	mu.set_gamma(config->MU);
	DOS = mu*D;	
	
	// Output to a file
	std::ofstream myfile;
	myfile.open ("DOS.dat");
	for(int i=0; i < energies.rows()*energies.cols(); i++){
		myfile  << energies(i)*config->energy_scale << " " << DOS(i).real() << " " << DOS(i).imag()<< "\n";
	}
	
	myfile.close();
}

template <typename U, unsigned DIM>
void calc_optical_cond(info<U,DIM> *config, int sigma, int R, int T, U beta, U e_fermi, std::string cond_axis){
	/* Calculates the optical conductivity for a set of frequencies in the range [-sigma, sigma].
	 * These frequencies are in the KPM scale, that is, the scale where the energy is in the range ]-1,1[.
	 * T is the number of frequency intervals in the range [0,sigma].
	 * R is the number of energy intervals inside each frequency interval.
	 * beta and e_fermi are not in the KPM scale, so they must be converted.
	 */
	
	// Convert beta and e_fermi to the KPM scale
	beta *= config->energy_scale;
	e_fermi /= config->energy_scale;
	
	// Calculate the number of frequencies and energies needed to perform the calculation.
	int N_omegas = 2*T*sigma + 1;
	int N_energies = 2*R*T*(sigma + 1) + 1;
	std::complex<U> imaginary(0.0, 1.0);
	
	Eigen::Matrix<U, -1, 1> energies;
	energies = Eigen::Matrix<U, -1, 1>::LinSpaced(N_energies, -sigma-1, sigma+1); 
	
	
	int n_moments = -1; // Will crash if used improperly
	if(cond_axis == "xx")
		n_moments = config->num_moments(config->CondXX);
	
	if(cond_axis == "xy")
		n_moments = config->num_moments(config->CondXY);
		
	
	
	// Dirac Delta and Green's functions
	tensor<1,1,U> Gr(energies, n_moments);
	tensor<1,1,U> Ga(energies, n_moments);
	tensor<1,1,U> D(energies, n_moments);
	
	// Matrices calculated from the kpm_transport program. They only depend on the Chebyshev moments
	tensor<0,2,U> Gamma(energies, n_moments);
	tensor<0,1,U> Lambda(energies, n_moments);
	
	// Objects calculated from those matrices that only depend on energy, not Chebyshev moments
	tensor<2,0,U> temp11(energies, n_moments);
	tensor<2,0,U> temp12(energies, n_moments);
	tensor<1,0,U> temp2(energies, n_moments);
	
	// The Green's functions are not calculated with a kernel, but with a finite imaginary part
	U finite_gamma = 0.01;
	Gr.set_green(-finite_gamma, "retarded");
	Ga.set_green( finite_gamma, "advanced");
	D.set_delta(kernel_jackson);
	
	// Choose which conductivity to calculate. XX or XY.
	if(cond_axis == "xx"){
		Gamma.set_gamma(config->GammaXX);
		Lambda.set_gamma(config->LambdaXX);
	}
	
	if(cond_axis == "xy"){
		Gamma.set_gamma(config->GammaXY);
		Lambda.set_gamma(config->LambdaXY);
	}
	
	// Perform the calculations to find the objects 
	temp11 = Gamma*Gr*D;
	temp12 = Gamma*D*Ga;
	temp2  = Lambda*D;
	
	// This is an auxiliary class to aid evaluating arrays
	value<U> e(energies);
	value<U> w(energies);
	
	// Array that's going to be integrated for each frequency
	Eigen::Array<std::complex<U>, -1, -1> integrand;
	Eigen::Array<U, -1, -1> integrand_energies;
	integrand = Eigen::Array<std::complex<U>, -1, -1>::Zero(1, 2*R*T + 1);
	integrand_energies = Eigen::Array<U, -1, -1>::Zero(1, 2*R*T + 1);
	
	// Frequency and conductivity arrays
	Eigen::Array<std::complex<U>, -1, -1> cond;
	cond = Eigen::Array<std::complex<U>, -1, -1>::Zero(1, N_omegas);
	Eigen::Array<U, -1, -1> freqs_array;
	freqs_array = Eigen::Array<U, -1, -1>::Zero(1, N_omegas);
	
	
	// Calculation of the optical conductivity
	for(int om = 0; om <= 2*T*sigma; om++){
		w.update_index(R*T + R*om);
		for(int i = 0; i <= 2*R*T; i++){			
			e.update_index(i + (N_energies-1)/2 - R*T);
			integrand_energies(i) = e.val;
			integrand(i) = fermi_function(e, e_fermi, beta)*(temp11(-e-w,e) + temp12(e,-e+w) + temp2(e))/(w.val + imaginary*finite_gamma);			
		}
		freqs_array(om) = w.val;
		cond(om) = imaginary*integrate(integrand_energies, integrand)*U(4.0*config->num_orbitals*config->spin_degeneracy/config->unit_cell_area);
	}
	
	// Output to a file
	std::ofstream myfile;
	myfile.open ("optical_cond.dat");
	for(int i=0; i < N_omegas; i++)
		myfile  << freqs_array(i)*config->energy_scale << " " << cond.real()(i) << " " << cond.imag()(i) << "\n";
	
	myfile.close();
}

template <typename T, unsigned DIM>
void single_shot(info<T,DIM> *config, Eigen::Array<T, -1, 1> energies){
	/* Calculates the DC XX conductivity at zero temperature.
	 * Evaluated at the values in the energies array.
	 */
	 
	int n_moments = config->num_moments(config->CondXX);
	
	tensor<1,1,T> imGa(energies, n_moments);
	tensor<0,2,T> GammaXX(energies, n_moments);
	tensor<2,0,T> cond_DC(energies, n_moments);
	
	T finite_gamma = 0.001;
	imGa.set_green_imag(finite_gamma, "advanced");
	GammaXX.set_gamma(config->GammaXX);
	cond_DC = GammaXX*imGa*imGa;
	
	// This is the constant that multiplies the conductivity. The formula is:
	// sigma = a*a*g_s*g_0/(A_w*pi)*(e*e/hbar)
	// a is the lattice constant that cancels with the one coming from A_w = unit_cell_area
	// g_s = spin_degeneracy accounts for the spin degeneracy, that was not contemplated by the Hamiltonian
	// g_0 is the number of orbitals per unit cell excluding spin. This appears due to our definition 
	//    of the random vectors. They are normalized by sqrt(N*g0), where N is the number of unit cells
	// The factor of 4 comes from the definition of the graphene universal conductivity
	// The minus sign comes from considering the imaginary part.
	
	double factor = -4*config->num_orbitals*config->spin_degeneracy/config->unit_cell_area/M_PI;
	
	std::ofstream myfile;
	myfile.open ("single_shot_matrix.dat");
	for(int i=0; i < energies.rows()*energies.cols(); i++)
		myfile  << energies(i) << " " << cond_DC(i,i).real()*factor << "\n";
	
	myfile.close();
	
	std::cout << "Left single_shot\n";fflush(stdout);
}


template <typename U, unsigned DIM>
void calculate(char *name){
	
	info<U,DIM> config("test_f.h5");
	config.read();
	
	if(config.DOS >= 0){
		std::cout << "Calculating the density of states... " << std::flush;
		double lim;
		int N_energies;
		calc_dos<U,DIM>(&config, lim = 1.0, N_energies = 1001);
		std::cout << "Done.\n" << std::flush;
	}
	
	if(config.CondXX >= 0){
		std::cout << "Calculating the XX optical conductivity... " << std::flush;
		int max_freq;
		int R;
		int T;
		U beta;
		U e_fermi;
		calc_optical_cond<U, DIM>(&config, max_freq = 2, R = 2, T=100, beta=200.0, e_fermi = 0.5, "xx");
		std::cout << "Done.\n" << std::flush;
	}
	
	if(config.CondXY >= 0){
		std::cout << "Calculating the XY optical conductivity... " << std::flush;
		int max_freq;
		int R;
		int T;
		U beta;
		U e_fermi;
		calc_optical_cond<U, DIM>(&config, max_freq = 2, R = 2, T=100, beta=200.0, e_fermi = 0.5, "xy");
		std::cout << "Done.\n" << std::flush;
	}
	
	
	// This is mostly used for testing, and as such is left commented out.
	/*
	int N_energies2 = 100;
	Eigen::Array<U, -1, 1> energies;
	energies = Eigen::Array<U, -1, 1>::Zero(N_energies2,1);
	for(int i=0; i < N_energies2; i++)
		energies(i) = -1 + 2.0*i/N_energies2;
	single_shot<U,DIM>(&config, energies);
	*/
}

