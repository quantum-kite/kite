#include "headers.hpp"

#include "H5Cpp.h"
#include "tensor.hpp"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"
#include "info.hpp"

void calc_dos(double, int);
void calc_optical_cond(int, int, int, double, double);
//void single_shot(Eigen::Array<double, -1, 1> energies);
void single_shot_matrix(Eigen::Array<double, -1, 1> energies);

//#################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

//TO DO:
//Make sure program crashes when the dataset read from h5 does not exist.

//https://support.hdfgroup.org/HDF5/doc/cpplus_RM/readdata_8cpp-example.html

double fermi_function(double energy, double mu, double beta){
	return 1.0/(1.0 + exp(beta*(energy - mu)));
}
	
double fermi_function(value<double> E1, double mu, double beta){
	return 1.0/(1.0 + exp(beta*(E1.val - mu)));
}
	
std::complex<double> integrate(Eigen::Array<double, -1, -1> energies, Eigen::Array<std::complex<double>, -1, -1> integrand){
	if(energies.rows() != integrand.rows() or energies.cols() != integrand.cols()){
		std::cout << "x and y arrays in the integrator must have the same number of elements. Exiting.\n";
		exit(0);
	}
	
	int N = energies.cols()*energies.rows();
	std::complex<double> sum(0,0);
	
	for(int i = 0; i < N - 1; i++){
		sum += (energies(i) - energies(i+1))*(integrand(i) + integrand(i+1))/2.0;
	}
	
	return sum;
}

int main(){
	double lim;
	int N_energies;
	calc_dos(lim = 1.0, N_energies = 101);
	
	Eigen::Array<double, -1, -1> energies;
	energies = Eigen::Array<double, -1, 1>::LinSpaced(100,-1,1);
	single_shot_matrix(energies);
	/*
	double sigma;
	int R;
	int T;
	double beta;
	double e_fermi;
	calc_optical_cond(sigma = 2.0, R = 2, T=100, beta=500, e_fermi = 0.0);
	*/
	
	return 0;
}


void calc_dos(double lim, int N_energies){
	Eigen::Matrix<double, Eigen::Dynamic, 1> energies;
	energies = Eigen::Matrix<double, Eigen::Dynamic, 1>::LinSpaced(N_energies, -lim, lim); 
	
	info config("test_f.h5");
	config.read();
	int n_moments = config.num_moments(0);
	
	tensor<1,1> D(energies, n_moments);
	tensor<0,1> mu(energies, n_moments);
	tensor<1,0> DOS(energies, n_moments);
	
	D.set_delta(kernel_jackson);
	mu.set_gamma(config.MU);
	DOS = mu*D;	
	
	std::ofstream myfile;
  myfile.open ("DOS.dat");
	for(int i=0; i < energies.rows()*energies.cols(); i++)
		myfile  << energies(i) << " " << DOS(i).real() << " " << DOS(i).imag()<< "\n";
	
	myfile.close();
}

void calc_optical_cond(int sigma, int R, int T, double beta, double e_fermi){
	int N_omegas = 2*T*sigma + 1;
	int N_energies = 2*R*T*(sigma + 1) + 1;
	std::complex<double> imaginary(0.0, 1.0);
	
	Eigen::Matrix<double, Eigen::Dynamic, 1> energies;
	energies = Eigen::Matrix<double, Eigen::Dynamic, 1>::LinSpaced(N_energies, -sigma-1, sigma+1); 
	
	// Read the data from the HDF5 file
	info config("data.h5");
	config.read();
	int n_moments = config.num_moments(0);
	
	// Dirac Delta and Green's functions
	tensor<1,1> Gr(energies, n_moments);
	tensor<1,1> Ga(energies, n_moments);
	tensor<1,1> D(energies, n_moments);
	
	// Matrices calculated from the kpm_transport program
	tensor<0,2> GammaXX(energies, n_moments);
	tensor<0,1> LambdaXX(energies, n_moments);
	
	tensor<2,0> temp11(energies, n_moments);
	tensor<2,0> temp12(energies, n_moments);
	tensor<1,0> temp2(energies, n_moments);
	
	double finite_gamma = 0.01;
	Gr.set_green(-finite_gamma, "retarded");
	Ga.set_green( finite_gamma, "advanced");
	D.set_delta(kernel_jackson);
	
	GammaXX.set_gamma(config.GammaXX);
	LambdaXX.set_gamma(config.LambdaXX);
	
	temp11 = GammaXX*Gr*D;
	temp12 = GammaXX*D*Ga;
	temp2  = LambdaXX*D;
	
	value<double> e(energies);
	value<double> w(energies);
	
	Eigen::Array<std::complex<double>, -1, -1> integrand;
	Eigen::Array<double, -1, -1> integrand_energies;
	integrand = Eigen::Array<std::complex<double>, -1, -1>::Zero(1, 2*R*T + 1);
	integrand_energies = Eigen::Array<double, -1, -1>::Zero(1, 2*R*T + 1);
	
	Eigen::Array<std::complex<double>, -1, -1> cond;
	cond = Eigen::Array<std::complex<double>, -1, -1>::Zero(1, N_omegas);
	Eigen::Array<double, -1, -1> freqs_array;
	freqs_array = Eigen::Array<double, -1, -1>::Zero(1, N_omegas);
	
	
	for(int om = 0; om <= 2*T*sigma; om++){
		w.update_index(R*T + R*om);
		for(int i = 0; i <= 2*R*T; i++){			
			e.update_index(i + (N_energies-1)/2 - R*T);
			integrand_energies(i) = e.val;
			integrand(i) = fermi_function(e, e_fermi, beta)*(temp11(-e-w,e) + temp12(e,-e+w) + temp2(e))/(w.val + imaginary*finite_gamma);			
		}
		freqs_array(om) = w.val;
		cond(om) = imaginary*integrate(integrand_energies, integrand);
	}
	
	std::ofstream myfile;
  myfile.open ("optical_cond.dat");
	for(int i=0; i < cond.rows()*cond.cols(); i++)
		myfile  << freqs_array(i) << " " << cond.real()(i) << " " << cond.imag()(i) << "\n";
	
	myfile.close();
}

void single_shot_matrix(Eigen::Array<double, -1, 1> energies){
	std::cout << "Calculating single-shot conductivity. This function requires the GammaXX matrix from the hdf5 file. The normalization factor requires knowledge of:\n";
	std::cout << " -spin degeneracy: g_s\n";
	std::cout << " -number of orbitals apart from spin: g_0\n";
	std::cout << " -area of the unit cell: unit_cell_area\n";
	std::cout << "All these should be obtainable from the hdf file, except the spin degeneracy. Please make sure that you need the spin degeneracy\n";
	
	// Read the data from the HDF5 file
	info config("test_f.h5");
	config.read();
	int n_moments = config.num_moments(0);
	
	tensor<1,1> imGa(energies, n_moments);
	tensor<0,2> GammaXX(energies, n_moments);
	tensor<2,0> cond_DC(energies, n_moments);
	
	double finite_gamma = 0.001;
	imGa.set_green_imag(finite_gamma, "advanced");
	GammaXX.set_gamma(config.GammaXX);
	cond_DC = GammaXX*imGa*imGa;
	
	// This is the constant that multiplies the conductivity. The formula is:
	// sigma = a*a*g_s*g_0/(A_w*pi)*(e*e/hbar)
	// a is the lattice constant that cancels with the one coming from A_w = unit_cell_area
	// g_s = spin_degeneracy accounts for the spin degeneracy, that was not contemplated by the Hamiltonian
	// g_0 is the number of orbitals per unit cell excluding spin. This appears due to our definition 
	//    of the random vectors. They are normalized by sqrt(N*g0), where N is the number of unit cells
	// The factor of 4 comes from the definition of the graphene universal conductivity
	// The minus sign comes from considering the imaginary part.
	
	double factor = -4*config.num_orbitals*config.spin_degeneracy/config.unit_cell_area/M_PI;
	
	std::ofstream myfile;
  myfile.open ("single_shot_matrix.dat");
	for(int i=0; i < energies.rows()*energies.cols(); i++)
		myfile  << energies(i) << " " << cond_DC(i,i).real()*factor << "\n";
	
	myfile.close();
	
	std::cout << "Left single_shot\n";fflush(stdout);
}
