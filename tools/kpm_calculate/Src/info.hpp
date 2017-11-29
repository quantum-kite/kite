#include "headers.hpp"

using std::cout;
using namespace H5;
using std::endl;

class info{
	H5::H5File file;
	public:
		int dim;
		int num_orbitals;
		Eigen::Array<int,1,-1> size;
		Eigen::Array<double,-1,-1> vectors;
		Eigen::Array<double, -1, -1> orbital_positions;
		Eigen::Array<int,1,-1> num_moments;
			
		double unit_cell_area;
		double spin_degeneracy;
		
		//Eigen::Array<std::complex<double>,-1,-1> GammaXoXX;
		//Eigen::Array<std::complex<double>,-1,-1> GammaXXoX;
		Eigen::Array<std::complex<double>,-1,-1> LambdaXX;
		Eigen::Array<std::complex<double>,-1,-1> GammaXX;
		Eigen::Array<std::complex<double>,-1,-1> MU;
		Eigen::Array<double,-1,-1> MUreal;
		
		info(std::string);
		void read();
	
};

info::info(std::string name){
	file = H5::H5File(name, H5F_ACC_RDONLY);
}
	


void info::read(){
	
	get_hdf5(&dim, &file, (char*)"DIM");
	
	//std::cout << "DIM\n";fflush(stdout);
	
	size = Eigen::Array<int,1,-1>::Zero(1,dim);
	get_hdf5(size.data(), &file, (char*)"L");
	
	//std::cout << "L\n";fflush(stdout);
	
	vectors = Eigen::Array<double,-1,-1>::Zero(dim,dim);
	get_hdf5(vectors.data(), &file, (char*)"LattVectors");
	unit_cell_area = fabs(vectors.matrix().determinant());
	std::cout << "cell_area: " << unit_cell_area;
	
	//std::cout << "LattVectors\n";fflush(stdout);
	
	get_hdf5(&num_orbitals, &file, (char*)"NOrbitals");
	//std::cout << "NOrbitals\n";fflush(stdout);
	
	orbital_positions = Eigen::Array<double,-1,-1>::Zero(num_orbitals, num_orbitals);
	get_hdf5(orbital_positions.data(), &file, (char*)"OrbPositions");
	//std::cout << "OrbPositions\n";fflush(stdout);
	
	// Time for the Gamma matrices
	
	num_moments = Eigen::Array<int, 1, -1>::Zero(1,3);
	get_hdf5(num_moments.data(), &file, (char*)"Calculation/NumMoments");
	//std::cout << "NumMoments: " << num_moments << "\n";fflush(stdout);
	
	
	//std::cout << "GammaXX\n";fflush(stdout);
	GammaXX = Eigen::Array<std::complex<double>,-1,-1>::Zero(1,num_moments(0)*num_moments(0));
	get_hdf5(GammaXX.data(), &file, (char*)"GammaXX");
	/*
	
	//std::cout << "LambdaXX\n";fflush(stdout);
	LambdaXX = Eigen::Array<std::complex<double>,-1,-1>::Zero(1,num_moments(0)*num_moments(0));
	get_hdf5(LambdaXX.data(), &file, (char*)"LambdaXX");
	*/
	//std::cout << "MU\n";fflush(stdout);
	MU = Eigen::Array<std::complex<double>,-1,-1>::Zero(1,num_moments(0));
	//MUreal = Eigen::Array<double,-1,-1>::Zero(1,num_moments(0));
	get_hdf5(MU.data(), &file, (char*)"MU");
	/*
	for(int i = 0; i < num_moments(0); i++){
		MU(i) = MUreal(i);
		std::cout << MU(i);
	}*/
	
	//std::cout << "finished MU\n";fflush(stdout);
	
	spin_degeneracy = 2;
	std::cout << "DIM: " << dim << "\n";
	std::cout << "L: " << size << "\n";
	std::cout << "vectors: " << vectors << "\n";
	std::cout << "num_orbitals: " << num_orbitals << "\n";
	std::cout << "positions of orbitals: " << orbital_positions << "\n";
	
}

