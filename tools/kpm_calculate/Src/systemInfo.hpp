#include "headers.hpp"

using std::cout;
using namespace H5;
using std::endl;

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
T fermi_function(T energy, T mu, T beta){
	return 1.0/(1.0 + exp(beta*(energy - mu)));
}

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
      dir = "xy"; break;
    case 2:
      dir = "xz"; break;
    case 3:
      dir = "yx"; break;
    case 4:
      dir = "yy"; break;
    case 5:
      dir = "yz"; break;
    case 6:
      dir = "zx"; break;
    case 7:
      dir = "zy"; break;
    case 8:
      dir = "zz"; break;
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


    int N_energies = energies.rows();
    std::cout << "N_energies:" << N_energies << "\n";

    T energy;
    Eigen::Matrix<std::complex<T>, -1, 1> GammaE;
    GammaE = Eigen::Matrix<std::complex<T>, -1, 1>::Zero(N_energies, 1);
    
    for(int e = 0; e < N_energies; e++){
      energy = energies(e);
      for(int m = 0; m < N_moments; m++){
        GammaE(e) += Gamma(m)*f0(m, energy);
      }
    }
    return integrate(energies, GammaE);
}

template <typename T>
Eigen::Matrix<std::complex<T>, -1, 1> contract2(
    std::function<T(int, T)> f0, int delta_position, 
    std::function<std::complex<T>(int, T)> f1, int N_moments, 
    const Eigen::Array<std::complex<T>, -1, -1>& Gamma, 
    const Eigen::Matrix<T, -1, 1>& energies, 
    const Eigen::Matrix<T, -1, 1>& frequencies){

    // First of all, contract the index that will not depend on the frequency.
    // That is the index of the delta function, that's why we need to know its
    // position.
    verbose_message("Entered contract 2d.\n");

    int N_energies = energies.rows();
    int N_freqs = frequencies.rows();

    T energy;
    Eigen::Matrix<std::complex<T>, -1, -1> GammaEN;
    GammaEN = Eigen::Matrix<std::complex<T>, -1, -1>::Zero(N_energies, N_moments);
   
    verbose_message("First part of the calculation.\n");

    if(delta_position == 0){
      for(int e = 0; e < N_energies; e++){
        energy = energies(e);
        for(int n = 0; n < N_moments; n++){
          for(int m = 0; m < N_moments; m++){
            GammaEN(e, n) += Gamma(n,m)*f0(m, energy);
            //std::cout << "e:" << energy;
            //std::cout << "n:" << n;
            //std::cout << "m:" << m;
            //std::cout << "GammaNN:" << Gamma(n.m);
            //std::cout << "f0:" << f0(m;
          }
        }
      }
    }

    
    else if(delta_position == 1){
      for(int e = 0; e < N_energies; e++){
        energy = energies(e);
        for(int n = 0; n < N_moments; n++){
          for(int m = 0; m < N_moments; m++){
            GammaEN(e, n) += Gamma(m,n)*f0(m, energy);
          }
        }
      }
    }
    else{
      std::cout << "delta_position must be specified. Exiting.\n";
      exit(1);
    }

    

    verbose_message("Second part of the calculation.\n");
    // Now contract the remaining index. This index does depend on the
    // frequency argument
    
    T freq;
    Eigen::Matrix<std::complex<T>, -1, 1> GammaE;
    GammaE = Eigen::Matrix<std::complex<T>, -1, 1>::Zero(N_energies, 1);
    
    Eigen::Matrix<std::complex<T>, -1, 1> cond;
    cond = Eigen::Matrix<std::complex<T>, -1, 1>::Zero(N_freqs, 1);

    for(int w1 = 0; w1 < N_freqs; w1++){
      freq = frequencies(w1);
      for(int e = 0; e < N_energies; e++){
        energy = energies(e) + freq;
        GammaE(e) = 0;
        for(int n = 0; n < N_moments; n++){
          GammaE(e) += GammaEN(e, n)*f1(n, energy);
        }
      }
      cond(w1) = integrate(energies, GammaE);
    }
    verbose_message("Left contract.\n");
    return cond;
}

template <typename T>
Eigen::Matrix<std::complex<T>, -1, 1> contract3(
    std::function<T(int, T)> f0, int delta_position, 
    std::function<std::complex<T>(int, T)> f1, int N_moments, 
    std::function<std::complex<T>(int, T)> f2, int N_moments2, 
    const Eigen::Array<std::complex<T>, -1, -1>& Gamma, 
    const Eigen::Matrix<T, -1, 1>& energies, 
    const Eigen::Matrix<T, -1, 1>& frequencies,
    const Eigen::Matrix<T, -1, 1>& frequencies2){

    // First of all, contract the index that will not depend on the frequency.
    // That is the index of the delta function, that's why we need to know its
    // position.
    verbose_message("Entered contract 2d.\n");

    int N_energies = energies.rows();
    int N_freqs1 = frequencies.rows();
    int N_freqs2 = frequencies2.rows();

    T energy;
    Eigen::Matrix<std::complex<T>, -1, -1> GammaENN;
    GammaENN = Eigen::Matrix<std::complex<T>, -1, -1>::Zero(N_energies, N_moments*N_moments);
   
    verbose_message("First part of the calculation.\n");

    if(delta_position == 0){
      for(int e = 0; e < N_energies; e++){
        energy = energies(e);
        for(int n = 0; n < N_moments; n++){
          for(int m = 0; m < N_moments; m++){
            for(int p = 0; p < N_moments; p++){
             GammaENN(e, n*N_moments + m) += Gamma(n*N_moments*N_moments + m*N_moments + p)*f0(p, energy);
            }
          }
        }
      }
    }
      

    if(delta_position == 1){
      for(int e = 0; e < N_energies; e++){
        energy = energies(e);
        for(int n = 0; n < N_moments; n++){
          for(int m = 0; m < N_moments; m++){
            for(int p = 0; p < N_moments; p++){
             GammaENN(e, n*N_moments + m) += Gamma(n*N_moments*N_moments + p*N_moments + m)*f0(p, energy);
            }
          }
        }
      }
    }
      

    if(delta_position == 2){
      for(int e = 0; e < N_energies; e++){
        energy = energies(e);
        for(int n = 0; n < N_moments; n++){
          for(int m = 0; m < N_moments; m++){
            for(int p = 0; p < N_moments; p++){
             GammaENN(e, n*N_moments + m) += Gamma(p*N_moments*N_moments + m*N_moments + n)*f0(p, energy);
            }
          }
        }
      }
    }
    else{
      std::cout << "delta_position must be specified. Exiting.\n";
      exit(1);
    }

    

    verbose_message("Second part of the calculation.\n");
    // Now contract the remaining index. This index does depend on the
    // frequency argument
    
    Eigen::Matrix<std::complex<T>, -1, 1> GammaEN;
    GammaEN = Eigen::Matrix<std::complex<T>, -1, 1>::Zero(N_energies, N_moments);
    Eigen::Matrix<std::complex<T>, -1, 1> GammaE;
    GammaE = Eigen::Matrix<std::complex<T>, -1, 1>::Zero(N_energies, 1);
    
    T freq1, freq2;

    Eigen::Matrix<std::complex<T>, -1, -1> cond;
    cond = Eigen::Matrix<std::complex<T>, -1, -1>::Zero(N_freqs1, N_freqs2);

    for(int w2 = 0; w2 < N_freqs2; w2++){
      freq2 = frequencies2(w2);
      for(int w1 = 0; w1 < N_freqs1; w1++){
        freq1 = frequencies(w1);


        for(int m = 0; m < N_moments; m++){
          for(int e = 0; e < N_energies; e++){
            energy = energies(e) + freq1;
            GammaEN(e, m) = 0;
            for(int n = 0; n < N_moments; n++){
              GammaEN(e, m) += GammaENN(e, m*N_moments + n)*f1(n, energy);
            }
          }
        }

        for(int e = 0; e < N_energies; e++){
          Gamma(e) = 0;
          energy = energies(e) + freq2;
          for(int m = 0; m < N_moments; m++){
            GammaE(e) += GammaEN(e,m)*f2(m, energy);
          }
        }
        cond(w1, w2) = integrate(energies, GammaE);
    
      }
    }
    verbose_message("Left contract.\n");
    return cond;
}

template <typename T, unsigned DIM>
class system_info{
	H5::H5File file;
	public:
		int dim;
		int num_orbitals;
    int isComplex;
		Eigen::Array<int,1,-1> size;
		Eigen::Array<double,-1,-1> vectors;
		Eigen::Array<double, -1, -1> orbital_positions;
			
		double unit_cell_area;
		double spin_degeneracy;
		double energy_scale;
    std::string filename;
	
    // These two quantities may only be known after the DOS is calculated
    bool EnergyLimitsKnown = false;
    T minEnergy, maxEnergy;

		system_info(std::string);
		system_info();
		void read();
	
};

template <typename T, unsigned DIM>
system_info<T, DIM>::system_info(){};

template <typename T, unsigned DIM>
system_info<T, DIM>::system_info(std::string name){
  filename = name;
	file = H5::H5File(name, H5F_ACC_RDONLY);
}
	




template <typename T, unsigned DIM>
void system_info<T, DIM>::read(){
	debug_message("Entered info::read.\n");
	// This function reads from the h5 file all the data that pertrains to
  // the Hamiltonian. Nothing about whether or not we need to calculate the
  // conductivity, or density of states, etc.
	
	// Basic information about the lattice 
	verbose_message("Reading basic information about the lattice: Dimension DIM, Length L and primitive lattice vectors LattVectors\n");
	dim = DIM;										// two-dimensional or three-dimensional
	size = Eigen::Array<int,1,-1>::Zero(1,dim);		// size of the sample
	get_hdf5(size.data(), &file, (char*)"L");
	get_hdf5(&isComplex, &file, (char*)"IS_COMPLEX"); // is the Hamiltonian a complex matrix?
	
	vectors = Eigen::Array<double,-1,-1>::Zero(dim,dim);	// Basis of primitive vectors that generate the lattice
	get_hdf5(vectors.data(), &file, (char*)"LattVectors");
	unit_cell_area = fabs(vectors.matrix().determinant());	// Use the basis vectors to determine the area of the unit cell
	
	verbose_message("Reading the energy scale, number of orbitals NOrbitals and their positions OrbPositions\n");
	get_hdf5(&energy_scale, &file, (char*)"EnergyScale");								// energy scale
	get_hdf5(&num_orbitals, &file, (char*)"NOrbitals");									// number of orbitals in each unit cell	
	orbital_positions = Eigen::Array<double,-1,-1>::Zero(num_orbitals, num_orbitals);	// position of each of those orbitals
	get_hdf5(orbital_positions.data(), &file, (char*)"OrbPositions");
	
	spin_degeneracy = 2; // put by hand?
	
	
	// Information about the data types
	verbose_message("Reading data type and checking whether it is complex.\n");
	int precision = 1, complex;
	get_hdf5(&complex, &file, (char *) "/IS_COMPLEX");
	get_hdf5(&precision,  &file, (char *) "/PRECISION");
	
	file.close();
	debug_message("Left info::read.\n");
}

