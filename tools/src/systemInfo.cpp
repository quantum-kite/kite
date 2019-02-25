/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

#include <Eigen/Dense>
#include <iostream>

#include "ComplexTraits.hpp"
#include "H5Cpp.h"
#include "myHDF5.hpp"

#include <complex>
#include <string>
#include "systemInfo.hpp"
#include "macros.hpp"

template <typename T, unsigned DIM>
system_info<T, DIM>::system_info(){
    EnergyLimitsKnown = false;
}

template <typename T, unsigned DIM>
system_info<T, DIM>::system_info(std::string name){
  filename = name;
    EnergyLimitsKnown = false;
}

template <typename T, unsigned DIM>
void system_info<T, DIM>::print_info(){
  if(VERBOSE == 1){
  std::cout << "Information about the lattice obtained from the configuration file:\n"
    "Dimension:           " << dim              << "\n"
    "Linear dimension:    " << size             << "\n"
    "Complex:             " << isComplex        << "\n"
    //"Lattice vectors:     " << LattVectors      << "\n"
    "Unit cell volume:    " << unit_cell_area   << "\n"
    "Energy scale:        " << energy_scale     << "\n"
    "Energy shift:        " << energy_shift     << "\n"
    //"Number of orbitals:  " << num_orbitals     << "\n"
    //"Orbital positions:   " << OrbPositions     << "\n"
    //"Divisions:           " << divisions        << "\n"
    "Number of threads:   " << NumThreads       << "\n";
    //"Precision:           " << precision        << "\n";
  }
}

template <typename T, unsigned DIM>
void system_info<T, DIM>::read(){
	debug_message("Entered info::read.\n");
	/* This function reads from the h5 file all the data that pertrains to
   * the Hamiltonian. Nothing about whether or not we need to calculate the
   * conductivity, or density of states, etc. */
	
	// Basic information about the lattice 
	debug_message("Reading basic information about the lattice: Dimension DIM, Length L and primitive lattice vectors LattVectors\n");


	file = H5::H5File(filename, H5F_ACC_RDONLY);
	dim = DIM;										// two-dimensional or three-dimensional
	size = Eigen::Array<int,1,-1>::Zero(1,dim);		// size of the sample
	get_hdf5(size.data(), &file, (char*)"L");
	get_hdf5(&isComplex, &file, (char*)"IS_COMPLEX"); // is the Hamiltonian a complex matrix?
	
	vectors = Eigen::Array<double,-1,-1>::Zero(dim,dim);	// Basis of primitive vectors that generate the lattice
	get_hdf5(vectors.data(), &file, (char*)"LattVectors");
	unit_cell_area = fabs(vectors.matrix().determinant());	// Use the basis vectors to determine the area of the unit cell
	
  debug_message("Reading the energy scale, number of orbitals NOrbitals and their positions OrbPositions\n");
  get_hdf5(&energy_scale, &file, (char*)"EnergyScale");								// energy scale
  get_hdf5(&energy_shift, &file, (char*)"EnergyShift");								// energy shift
  get_hdf5(&num_orbitals, &file, (char*)"NOrbitals");									// number of orbitals in each unit cell	
  orbital_positions = Eigen::Array<double,-1,-1>::Zero(num_orbitals, dim);	// position of each of those orbitals
  get_hdf5(orbital_positions.data(), &file, (char*)"OrbPositions");

	
  Eigen::Array<int, -1, -1> divisions;
  divisions = Eigen::Array<int, -1, -1>::Zero(dim, 1);

  get_hdf5(divisions.data(), &file, (char*)"Divisions");
  NumThreads = 1;
  for(int ntr = 0; ntr < dim; ntr++)
    NumThreads *= divisions(ntr);
  debug_message("NumThreads: "); debug_message(NumThreads); debug_message("\n");

  spin_degeneracy = 1; // put by hand?
  
   //Information about the data types
  //debug_message("Reading data type and checking whether it is complex.\n");
  //int precision = 1, complex;
  //get_hdf5(&complex, &file, (char *) "/IS_COMPLEX");
  //get_hdf5(&precision,  &file, (char *) "/PRECISION");
	
	file.close();
	debug_message("Left info::read.\n");
}

// Instantiations

template class system_info<float, 1u>;
template class system_info<float, 2u>;
template class system_info<float, 3u>;

template class system_info<double, 1u>;
template class system_info<double, 2u>;
template class system_info<double, 3u>;

template class system_info<long double, 1u>;
template class system_info<long double, 2u>;
template class system_info<long double, 3u>;
