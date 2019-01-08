/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

#include "headers.hpp"

using std::cout;
using namespace H5;
using std::endl;


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

int custom_find(int *arr, int arr_size, int value_to_search){
	/* This function searches array 'arr' for any occurence of number 'value to search'.
	 * If it exists, it exits and returns the index where it occurs.
	 * Otherwise, returns -1
	 */
	for(int i = 0; i < arr_size; i++){
		if(*(arr+i) == value_to_search)
			return i;
	}
	
	return -1;	
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
		
		system_info(std::string);
		system_info();
		void read();
	
};

template <typename T, unsigned DIM>
system_info<T, DIM>::system_info(){};

template <typename T, unsigned DIM>
system_info<T, DIM>::system_info(std::string name){
	file = H5::H5File(name, H5F_ACC_RDONLY);
}
	




template <typename T, unsigned DIM>
void system_info<T, DIM>::read(){
	debug_message("Entered info::read.\n");
	// This function reads from the h5 file all the data that pertrains to
  // the Hamiltonian. Nothing about whether or not we need to calculate the
  // conductivity, or density of states, etc.
	
	// Basic information about the lattice 
	debug_message("Reading basic information about the lattice: Dimension DIM, Length L and primitive lattice vectors LattVectors\n");
	dim = DIM;										// two-dimensional or three-dimensional
	size = Eigen::Array<int,1,-1>::Zero(1,dim);		// size of the sample
	get_hdf5(size.data(), &file, (char*)"L");
	get_hdf5(&isComplex, &file, (char*)"IS_COMPLEX"); // is the Hamiltonian a complex matrix?
	
	vectors = Eigen::Array<double,-1,-1>::Zero(dim,dim);	// Basis of primitive vectors that generate the lattice
	get_hdf5(vectors.data(), &file, (char*)"LattVectors");
	unit_cell_area = fabs(vectors.matrix().determinant());	// Use the basis vectors to determine the area of the unit cell
	
	debug_message("Reading the energy scale, number of orbitals NOrbitals and their positions OrbPositions\n");
	get_hdf5(&energy_scale, &file, (char*)"EnergyScale");								// energy scale
	get_hdf5(&num_orbitals, &file, (char*)"NOrbitals");									// number of orbitals in each unit cell	
	orbital_positions = Eigen::Array<double,-1,-1>::Zero(num_orbitals, num_orbitals);	// position of each of those orbitals
	get_hdf5(orbital_positions.data(), &file, (char*)"OrbPositions");
	
	spin_degeneracy = 1; // put by hand?
	
	
	// Information about the data types
	debug_message("Reading data type and checking whether it is complex.\n");
	int precision = 1, complex;
	get_hdf5(&complex, &file, (char *) "/IS_COMPLEX");
	get_hdf5(&precision,  &file, (char *) "/PRECISION");
	
	file.close();
	debug_message("Left info::read.\n");
}



