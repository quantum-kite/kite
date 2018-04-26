#include "headers.hpp"

/*
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

template <typename TR>
TR delta(int n, TR energy){
	TR sq = sqrt(1.0 - energy*energy);
	if(energy < 1 and energy > -1)
		return 2.0/M_PI/sq*cos(n*acos(energy));
	else
		return 0;
}
*/



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

template <typename T>
std::complex<T> contract(std::complex<T>(*f)(int, std::complex<T>), T omega, int N_moments,
                         Eigen::Matrix<std::complex<T>, -1, -1> gammaEN, Eigen::Matrix<T, -1, 1> energies) {

  // This function performs the multiplication of the gamma matrix with the functions
  // specified in the arguments. To increase efficiency, one very important aspect
  // must be considered. One of the functions multiplying by the original Gamma matrix
  // is always a delta function (it originates from the Fermi function). The important
  // fact is that it does not depend on the frequencies at all. Therefore, we can
  // perform this multiplication beforehand and then provide THAT resulting product 
  // to THIS function. That product is called gammaEN to reflect the fact that it
  // has one energy argument and one Chebychev argument. 
  //
  // It is assumed that gammaEN already includes the Fermi function
  

  int N_energies = energies.rows();
  double energy;



  // Second and last step, contract the final Chebyshev index, the left-most one
  // This yields a matrix with dimensions (N_energies, 1)
  Eigen::Matrix<std::complex<T>, -1, -1> gammaE;
  gammaE = Eigen::Matrix<std::complex<T>, -1, -1>::Zero(N_energies, 1);

  for(int e = 0; e < N_energies; e++){
    energy = energies(e);
    for(int n = 0; n < N_moments; n++){
      gammaE(e) += f(n, energy + omega)*gammaEN(e, n)/(1.0 + double(n==0));
    }
  }

  

  // All the indices have now been contracted. We may proceed to the integration
  return integrate(energies, gammaE);
}



template <typename T>
std::complex<T> contract(std::complex<T>(*f1)(int, std::complex<T>), T omega1, int N_moments1,
                         std::complex<T>(*f2)(int, std::complex<T>), T omega2, int N_moments2,
                         Eigen::Matrix<std::complex<T>, -1, -1> gammaENN, Eigen::Matrix<T, -1, 1> energies) {

  // This function performs the multiplication of the gamma matrix with the functions
  // specified in the arguments. To increase efficiency, one very important aspect
  // must be considered. One of the functions multiplying by the original Gamma matrix
  // is always a delta function (it originates from the Fermi function). The important
  // fact is that it does not depend on the frequencies at all. Therefore, we can
  // perform this multiplication beforehand and then provide THAT resulting product 
  // to THIS function. That product is called gammaENN to reflect the fact that it
  // has one energy argument and two Chebychev arguments. 
  //
  // It is assumed that gammaENN already includes the Fermi function
  

  int N_energies = energies.rows();
  double energy;


  // First step, contract one of the Chebyshev indices, the right-most one
  // This yields a matrix with dimensions (N_energies, N_moments1)
  Eigen::Matrix<std::complex<T>, -1, -1> gammaEN;
  gammaEN = Eigen::Matrix<std::complex<T>, -1, -1>::Zero(N_energies, N_moments1);

  for(int e = 0; e < N_energies; e++){
    energy = energies(e);
    for(int n = 0; n < N_moments1; n++){
      for(int m = 0; m < N_moments2; m++){
        gammaEN(e,n) += f2(m, energy + omega2)*gammaENN(e,n*N_moments1 + m)/(1.0 + double(m==0));
      }
    }
  }





  // Second and last step, contract the final Chebyshev index, the left-most one
  // This yields a matrix with dimensions (N_energies, 1)
  Eigen::Matrix<std::complex<T>, -1, -1> gammaE;
  gammaE = Eigen::Matrix<std::complex<T>, -1, -1>::Zero(N_energies, 1);

  for(int e = 0; e < N_energies; e++){
    energy = energies(e);
    for(int n = 0; n < N_moments1; n++){
      gammaE(e) += f1(n, energy + omega1)*gammaEN(e, n)/(1.0 + double(n==0));
    }
  }

  

  // All the indices have now been contracted. We may proceed to the integration
  return integrate(energies, gammaE);
}



template <typename T>
T fermi_function(T energy, T mu, T beta){
	return 1.0/(1.0 + exp(beta*(energy - mu)));
}


template<typename T, unsigned DIM>
void calc_dos_simple(info<T,DIM> *config, double lim, int N_energies){
  // Calculates the density of states
  debug_message("Entered calc_dos_simple.\n");




  int N_moments = config->num_moments(0);

  Eigen::Matrix<T, -1, 1> energies;
  Eigen::Matrix<std::complex<T>, -1, 1> DOS;
	energies = Eigen::Matrix<T, -1, 1>::LinSpaced(N_energies, -lim, lim); 
	DOS = Eigen::Matrix<std::complex<T>, -1, 1>::Zero(N_energies); 
   


  for(int i = 0; i < N_moments; i++){
    for(int j = 0; j < N_energies; j++){
      DOS(j) += config->MU(i)*delta(i, energies(j))*kernel_jackson(i,N_moments)/(1.0 + (double)(i==0));
    }
  }


  
	std::ofstream myfile;
	myfile.open ("DOS.dat");
	for(int i=0; i < energies.rows()*energies.cols(); i++){
		myfile  << energies(i)*config->energy_scale << " " << DOS(i).real() << " " << DOS(i).imag()<< "\n";
	}
	
	myfile.close();
	debug_message("Left calc_dos.\n");
};


template <typename U, unsigned DIM>
void calc_optical_cond_simple_explicit(info<U,DIM> *config, std::string cond_axis){
	debug_message("Entered calc_optical_cond.\n");
	/* Calculates the optical conductivity for a set of frequencies in the range [-sigma, sigma].
	 * These frequencies are in the KPM scale, that is, the scale where the energy is in the range ]-1,1[.
	 * T is the number of frequency intervals in the range [0,sigma].
	 * R is the number of energy intervals inside each frequency interval.
	 * beta and e_fermi are not in the KPM scale, so they must be converted.
	 */
	double beta = 200;
  double e_fermi = 0.0;

	// Convert beta and e_fermi to the KPM scale
	beta *= config->energy_scale;
	e_fermi /= config->energy_scale;
	
	// Calculate the number of frequencies and energies needed to perform the calculation.
	int N_omegas = 500;
	int N_energies = 1001;
  double lim = 0.999;
	std::complex<U> imaginary(0.0, 1.0);

  Eigen::Matrix<U, -1, 1> energies;
  int N_moments;

  energies  = Eigen::Matrix<U, -1, 1>::LinSpaced(N_energies, -lim, lim);
  

  /*Eigen::Array<std::complex<U>, -1, -1> Gamma;
  Eigen::Array<std::complex<U>, -1, -1> Lambda;
  Gamma  = Eigen::Array<std::complex<U>, -1, -1>::Zero(1, N_moments*N_moments);
  Lambda = Eigen::Array<std::complex<U>, -1, -1>::Zero(1, N_moments);*/
  
  std::complex<U> *pointer_gamma;
  std::complex<U> *pointer_lambda;
	if(cond_axis == "xx"){
		//N_moments.at(0) = config->num_moments(config->CondXX);
		//N_moments.at(1) = config->num_moments(config->CondXX);
		N_moments = config->num_moments(config->CondXX);
    pointer_gamma = &config->GammaXX(0);
    pointer_lambda = &config->LambdaXX(0);
  }
	
	if(cond_axis == "xy"){
		//N_moments.at(0) = config->num_moments(config->CondXY);
		//N_moments.at(1) = config->num_moments(config->CondXY);
		N_moments = config->num_moments(config->CondXY);
    pointer_gamma = &config->GammaXY(0);
    pointer_lambda = &config->LambdaXY(0);
  }
  
  else{
    std::cout << "Unknown direction in calc_optical_cond_simple.\n";
  }
  
  Eigen::Map<Eigen::Matrix<std::complex<U>, 1, -1>> Gamma(pointer_gamma, N_moments*N_moments);
  Eigen::Map<Eigen::Matrix<std::complex<U>, 1, -1>> Lambda(pointer_lambda, N_moments*N_moments);
  
  /*
  for(int i = 0; i < N_moments; i++){
    for(int j = 0; j < N_moments; j++){
      std::cout << Gamma(i*N_moments + j) << " ";
    }
    std::cout << "\n";
  }*/

  




  // Perform the first part of the multiplication. This is the multiplication of
  // one of the indices of the Gamma matrix by the Delta function. This does not
  // depend on the frequencies, so it may be done separately to increase speed
  
  Eigen::Matrix<std::complex<U>, -1, -1> DG;
  Eigen::Matrix<std::complex<U>, -1, -1> GD;
  DG = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_moments);
  GD = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_moments);
  double factor;

  for(int n = 0; n < N_moments; n++){
    for(int e = 0; e < N_energies; e++){
      for(int m = 0; m < N_moments; m++){
        factor = delta(m, energies(e))*kernel_jackson(m, N_moments)/(1.0 + (double)(m==0));
        DG(e,n) += factor*Gamma(m + n*N_moments);
        GD(e,n) += factor*Gamma(n + n*N_moments);
      }
    }
  }





  // Now perform the part of the calculation that does depend on the frequency.

  Eigen::Matrix<U, -1, 1> frequencies;
  frequencies = Eigen::Matrix<U, -1, 1>::LinSpaced(N_omegas, -lim, lim);
  Eigen::Matrix<std::complex<U>, -1, 1> integrand;
  Eigen::Matrix<std::complex<U>, -1, 1> cond;

  integrand = Eigen::Matrix<std::complex<U>, -1, 1>::Zero(N_energies, 1);
  double freq;

  double scat = 0.001;
  std::complex<U> complex_energy_p;
  std::complex<U> complex_energy_m;

  for(int w = 0; w < N_omegas; w++){
    freq = frequencies[w];
    for(int e = 0; e < N_energies; e++){
      complex_energy_p = std::complex<U> (energies(e) + freq,  scat);
      complex_energy_m = std::complex<U> (energies(e) - freq, -scat);
      integrand(e) = 0;
      for(int m = 0; m < N_moments; m++){
        integrand(e) += fermi_function(energies(e), e_fermi, beta)/(1.0 + (U)(m==0))*
          (GD(e,m)*green(m,  1, complex_energy_p) + DG(e,m)*green(m, -1, complex_energy_m));
      }
    }


    cond(w) = imaginary*integrate(energies, integrand)/(freq + imaginary*scat)*U(4.0*
        config->num_orbitals*config->spin_degeneracy/config->unit_cell_area);

  }




  


  // Second part of the multiplication, diamagnetic term
  double energy;
  std::complex<U> diam(0,0);
  for(int n = 0; n < N_moments; n++){
    for(int e = 0; e < N_energies; e++){
      energy = energies(e);
      integrand(e) = fermi_function(energy, e_fermi, beta)*delta(n, energy)*
                      kernel_jackson(n, N_moments)*Lambda(n);
    }

    diam += integrate(energies, integrand);
  }

  for(int w = 0; w < N_omegas; w++){

    freq = frequencies[w];
    cond(w) += imaginary*diam/(freq + imaginary*scat)*
      U(4.0*config->num_orbitals*config->spin_degeneracy/config->unit_cell_area);
  }





	// Output to a file
	std::ofstream myfile;
	myfile.open ("optical_cond.dat");
	for(int i=0; i < N_omegas; i++)
		myfile  << frequencies[i]*config->energy_scale << " " << cond.real()(i) << " " << cond.imag()(i) << "\n";
	
	myfile.close();
	debug_message("Left calc_optical_cond.\n");




      
    

  


  


  

  
};



template <typename U, unsigned DIM>
void calc_optical_cond_simple(info<U,DIM> *config, std::string cond_axis){
	debug_message("Entered calc_optical_cond.\n");
	/* Calculates the optical conductivity for a set of frequencies in the range [-sigma, sigma].
	 * These frequencies are in the KPM scale, that is, the scale where the energy is in the range ]-1,1[.
	 * T is the number of frequency intervals in the range [0,sigma].
	 * R is the number of energy intervals inside each frequency interval.
	 * beta and e_fermi are not in the KPM scale, so they must be converted.
	 */
	double beta = 200;
  double e_fermi = 0.0;

	// Convert beta and e_fermi to the KPM scale
	beta *= config->energy_scale;
	e_fermi /= config->energy_scale;
	
	// Calculate the number of frequencies and energies needed to perform the calculation.
	int N_omegas = 500;
	int N_energies = 1001;
  double lim = 0.999;
	std::complex<U> imaginary(0.0, 1.0);

  Eigen::Matrix<U, -1, 1> energies;
  int N_moments;

  energies  = Eigen::Matrix<U, -1, 1>::LinSpaced(N_energies, -lim, lim);
  

  
  std::complex<U> *pointer_gamma;
  std::complex<U> *pointer_lambda;
	if(cond_axis == "xx"){
		//N_moments.at(0) = config->num_moments(config->CondXX);
		//N_moments.at(1) = config->num_moments(config->CondXX);
		N_moments = config->num_moments(config->CondXX);
    pointer_gamma = &config->GammaXX(0);
    pointer_lambda = &config->LambdaXX(0);
  }
	
	if(cond_axis == "xy"){
		//N_moments.at(0) = config->num_moments(config->CondXY);
		//N_moments.at(1) = config->num_moments(config->CondXY);
		N_moments = config->num_moments(config->CondXY);
    pointer_gamma = &config->GammaXY(0);
    pointer_lambda = &config->LambdaXY(0);
  }
  
  else{
    std::cout << "Unknown direction in calc_optical_cond_simple.\n";
  }
  
  Eigen::Map<Eigen::Matrix<std::complex<U>, 1, -1>> Gamma(pointer_gamma, N_moments*N_moments);
  Eigen::Map<Eigen::Matrix<std::complex<U>, 1, -1>> Lambda(pointer_lambda, N_moments*N_moments);
  
  /*
  for(int i = 0; i < N_moments; i++){
    for(int j = 0; j < N_moments; j++){
      std::cout << Gamma(i*N_moments + j) << " ";
    }
    std::cout << "\n";
  }*/

  
  





	// Output to a file
	std::ofstream myfile;
	myfile.open ("optical_cond.dat");
	for(int i=0; i < N_omegas; i++)
		myfile  << frequencies[i]*config->energy_scale << " " << cond.real()(i) << " " << cond.imag()(i) << "\n";
	
	myfile.close();
	debug_message("Left calc_optical_cond.\n");




      
    

  


  


  

  
};



template <typename U, unsigned DIM>
void calculate_simple(char *name){
	debug_message("Entered calculate_simple.\n");

	info<U,DIM> config(name);
	config.read();
	
	if(config.DOS >= 0){
		double lim;
		int N_energies = 1001;
		verbose_message("Calculating the density of states with "); verbose_message(N_energies); verbose_message(" energies... ");
		calc_dos_simple<U,DIM>(&config, lim = 1.0, N_energies);
		verbose_message("Done.\n");
	}
	
  if(config.CondXX >= 0){
    verbose_message("calculating the optical conductivity");
    calc_optical_cond_simple(&config, "xx");
  }
	debug_message("Left calculate_simple.\n");
}



