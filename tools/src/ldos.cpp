/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <string>
#include <omp.h>

#include "H5Cpp.h"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"

#include "parse_input.hpp"
#include "systemInfo.hpp"
#include "ldos.hpp"

#include "functions.hpp"
#include "macros.hpp"

template <typename T, unsigned DIM>
ldos<T, DIM>::ldos(system_info<T, DIM>& sysinfo, shell_input & vari){
    // Class constructor
    
    systemInfo  = &sysinfo;              // retrieve the information about the Hamiltonian
    variables   = vari;                   // retrieve the shell input
    dirName     = "/Calculation/ldos/";     // location of the information about the conductivity
    
    isRequired = is_required() && variables.lDOS_is_required;         // was the local density of states requested?
    isPossible = false;                 // do we have all we need to calculate the density of states?
    if(isRequired){
        set_default_parameters();
        isPossible = fetch_parameters();
        override_parameters();

      if(isPossible){
          printLDOS();                  // Print all the parameters used
          calculate();
      } else {
        std::cout << "ERROR. The LDOS was requested but the data "
            "needed for its computation was not found in the input .h5 file. "
            "Make sure KITEx has processed the file first. Exiting.";
        exit(1);
      }
    }
}

template <typename T, unsigned DIM>
void ldos<T, DIM>::printLDOS(){
  double scale = systemInfo->energy_scale;
    std::cout << "The local density of states will be calculated with the following parameters:\n"
        "   Number of energies: " << NumEnergies << "\n"
        "   Number of positions: " << NumPositions << "\n"
        "   Filename: " << filename  << "X.dat" << ((default_filename)?" (default)":"") << "\n"
        "   Kernel: "               << kernel           << ((default_kernel)?           " (default)":"") << "\n";
    if(kernel == "green"){
        std::cout << "   Kernel parameter: "     << kernel_parameter*scale << ((default_kernel_parameter)? " (default)":"") << "\n";
    }
}

template <typename T, unsigned DIM>
bool ldos<T, DIM>::is_required(){
    // check whether the local density of states was asked for
    // if this quantity exists, so should all the others.

    name = systemInfo->filename;
	H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);
    bool result = false;
    try{
        H5::Exception::dontPrint();
        get_hdf5(&NumMoments, &file, (char*)(dirName+"NumMoments").c_str());									
        result = true;
    } catch(H5::Exception& e){}
  

    file.close();

    return result;
}
	
template <typename T, unsigned DIM>
void ldos<T, DIM>::override_parameters(){
    if(variables.lDOS_Name != ""){
        filename         = variables.lDOS_Name;
        default_filename = false;
    }

    if(variables.lDOS_NumMoments != -1){
        NumMoments         = variables.lDOS_NumMoments;
        default_NumMoments = false;
        if(variables.lDOS_NumMoments > MaxMoments){
          std::cout << "lDOS: The number of Chebyshev moments specified"
            " cannot be larger than the number of moments calculated by KITEx."
            " Please specify a smaller number. Exiting.\n";
          exit(1);
        }
    }

    //std::cout << "variables kernel:" << variables.lDOS_kernel << "\n";
    if(variables.lDOS_kernel != ""){
      //std::cout << "entered if \n";
        kernel         = variables.lDOS_kernel;
        default_kernel = false;
    }
    //std::cout << "kernel: " << kernel << "\n";

    if(kernel == "green"){
      if(variables.lDOS_kernel_parameter != -8888.8){
        kernel_parameter = variables.lDOS_kernel_parameter/systemInfo->energy_scale;
        default_kernel_parameter = false;
      }
    }
}

template <typename T, unsigned DIM>
void ldos<T, DIM>::set_default_parameters(){
    filename = "ldos";
    default_filename = true;
    MaxMoments = -1;

    // kernel options
    kernel = "jackson";
    default_kernel = true;
    default_kernel_parameter = true;

}


template <typename T, unsigned DIM>
bool ldos<T, DIM>::fetch_parameters(){
	debug_message("Entered ldos::fetch_parameters.\n");
	//This function reads all the data from the hdf5 file that's needed to 
    //calculate the LDoS
	 
    // Check if the data for the ldos exists
    if(!isRequired){
        std::cout << "Data for LDoS does not exist. Exiting.\n";
        exit(1);
    }
  
    H5::DataSet * dataset;
    H5::DataSpace * dataspace;
    hsize_t dim[2];
	H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

    dataset            = new H5::DataSet(file.openDataSet("/Calculation/ldos/Orbitals")  );
    dataspace          = new H5::DataSpace(dataset->getSpace());
    dataspace -> getSimpleExtentDims(dim, NULL);
    dataspace->close(); delete dataspace;
    dataset->close();   delete dataset;
    NumPositions = dim[0];

    dataset            = new H5::DataSet(file.openDataSet("/Calculation/ldos/Energy")  );
    dataspace          = new H5::DataSpace(dataset->getSpace());
    dataspace -> getSimpleExtentDims(dim, NULL);
    dataspace->close(); delete dataspace;
    dataset->close();   delete dataset;
    NumEnergies = dim[0];

    ldos_Orbitals = Eigen::Matrix<unsigned long, -1, -1>::Zero(NumPositions,1);
    ldos_Positions = Eigen::Matrix<unsigned long, -1, -1>::Zero(NumPositions,1);
    energies = Eigen::Matrix<float, -1, -1>::Zero(NumEnergies,1);

     //Fetch the relevant parameters from the hdf file
    get_hdf5(&MaxMoments, &file, (char*)(dirName+"NumMoments").c_str());	
    get_hdf5(ldos_Orbitals.data(), &file, (char*)"/Calculation/ldos/Orbitals");
    get_hdf5(ldos_Positions.data(), &file, (char*)"/Calculation/ldos/FixPosition");
    get_hdf5(energies.data(), &file, (char*)"/Calculation/ldos/Energy");


    global_positions = Eigen::Matrix<unsigned long, -1, -1>::Zero(NumPositions,3);
    for(long i = 0; i < NumPositions; i++){
      int L = systemInfo->size[0];
      global_positions(i,0) = ldos_Positions(i)%L;
      global_positions(i,1) = ldos_Positions(i)/L;
      global_positions(i,2) = ldos_Orbitals(i);
    }



  // Check whether the matrices we're going to retrieve are complex or not
  int complex = systemInfo->isComplex;

  bool result = false;
  // Retrieve the lmu Matrix
  std::string MatrixName = dirName + "lMU";
  try{
		debug_message("Filling the lMU matrix.\n");
		lMU = Eigen::Matrix<std::complex<T>,-1,-1>::Zero(MaxMoments, NumPositions);
		
		if(complex)
			get_hdf5(lMU.data(), &file, (char*)MatrixName.c_str());
		if(!complex){
			Eigen::Matrix<T,-1,-1> lMUReal; 
			lMUReal = Eigen::Matrix<T,-1,-1>::Zero(MaxMoments, NumPositions); 
			get_hdf5(lMUReal.data(), &file, (char*)MatrixName.c_str()); 

      lMU = lMUReal.template cast<std::complex<T>>();
		}				

    result = true;
  } catch(H5::Exception& e) {debug_message("lDOS: There is no lMU matrix.\n");}
	

  NumMoments = MaxMoments;
  //std::cout << "NumMoments: " << NumMoments << "\n";

  file.close();
	debug_message("Left lDOS::fetch_parameters.\n");
  return result;
}

template <typename U, unsigned DIM>
void ldos<U, DIM>::calculate(){

  Eigen::Matrix<std::complex<U>, -1, -1> LDOS;
  LDOS = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(NumEnergies, NumPositions);

  Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> OrderedMU;
  OrderedMU = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>::Zero(NumEnergies, NumPositions);
  OrderedMU = lMU;
  //std::cout << "lMU: \n" << OrderedMU << "\n";

  omp_set_num_threads(systemInfo->NumThreads);
  //omp_set_num_threads(1);
#pragma omp parallel 
{
#pragma omp critical
{
  int localN = NumMoments/systemInfo->NumThreads;
  //int localN = NumMoments;
  int thread_id = omp_get_thread_num();
  //std::cout << "thread_id: " << thread_id << "\n";
  //std::cout << "NumMoments: " << localN << "\n";
  //thread_id = 0;
  long offset = thread_id*localN*NumPositions;
  Eigen::Map<Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>> locallMU(OrderedMU.data() + offset, localN, NumPositions);

  //std::cout << "locallMU: \n" << locallMU << "\n";

  Eigen::Matrix<std::complex<U>, -1, -1> GammaE;
  GammaE = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(NumEnergies, localN);

  U factor;

  if(kernel == "jackson"){
    std::cout << "Entered jackson kernel\n";
    for(int i = 0; i < NumEnergies; i++){
      for(int m = 0; m < localN; m++){
        factor = 1.0/(1.0 + U((m + thread_id*localN)==0));
        GammaE(i,m) += delta(m + thread_id*localN,energies(i))*kernel_jackson<U>(m + thread_id*localN, NumMoments)*factor;
      }
    }
  }


  if(kernel == "green"){
    std::complex<U> c_energy;
    for(int i = 0; i < NumEnergies; i++){
      c_energy = std::complex<U>(energies(i), kernel_parameter);
      for(int m = 0; m < localN; m++){
        factor = 1.0/(1.0 + U((m + thread_id*localN)==0));
        GammaE(i,m) += -factor*green<std::complex<U>>(m, 1, c_energy).imag();
      }
    }
  }

  //std::cout << "GammaE: \n" << GammaE << "\n";
  //for(int m = 0; m < ; m++){
    //factor = 1.0/(1.0 + U(m==0));
    //for(int i = 0; i < NumEnergies; i++){
      //GammaE(i,m) = delta(m + thread_id*localN, energies(i))*kernel_jackson<U>(m + thread_id*localN, NumMoments)*factor;
    //}
  //}

  Eigen::Matrix<std::complex<U>, -1, -1> localLDOS;
  localLDOS = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(NumEnergies, NumPositions);
  localLDOS = GammaE*locallMU;
//#pragma omp critical
  LDOS += localLDOS;
}
}

  // Save the density of states to a file
  U mult = 1.0/systemInfo->energy_scale;
  std::ofstream myfile;
  double scale = systemInfo->energy_scale;
  double shift = systemInfo->energy_shift;
  for(int i=0; i < NumEnergies; i++){
    myfile.open(filename + std::to_string(energies(i)*scale + shift) + ".dat");
    for(unsigned pos = 0; pos < NumPositions; pos++){
      int x, y, orb;
      x = global_positions(pos,0);
      y = global_positions(pos,1);
      orb = global_positions(pos,2);
      myfile  << x << " " << y << " " << orb << " " << LDOS(i,pos).real()*mult << "\n";
    }
  myfile.close();
  }
}




// Instantiations
template class ldos<float, 1u>;
template class ldos<float, 2u>;
template class ldos<float, 3u>;

template class ldos<double, 1u>;
template class ldos<double, 2u>;
template class ldos<double, 3u>;

template class ldos<long double, 1u>;
template class ldos<long double, 2u>;
template class ldos<long double, 3u>;
