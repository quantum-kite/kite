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
#include "H5Cpp.h"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"

#include "parse_input.hpp"
#include "systemInfo.hpp"
#include "dos.hpp"
#include "functions.hpp"

#include "macros.hpp"

//template <typename T, unsigned DIM>
//dos<T, DIM>::dos(){};

template <typename T, unsigned DIM>
dos<T, DIM>::dos(system_info<T, DIM>& sysinfo, shell_input & vari){
    std::string name = sysinfo.filename;
	file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

    isRequired = false; // was this quantity density of states asked for?
    isPossible = false; // do we have all we need to calculate the density of states?

    NumPoints = -1;
    filename = "dos.dat";

    // retrieve the information about the Hamiltonian
    systemInfo = &sysinfo;

    // retrieve the shell input
    variables = vari;

    // location of the information about the conductivity
    dirName = "/Calculation/dos/";

    // check whether the conductivity_dc was asked for
    // if this quantity exists, so should all the others.
    try{
        H5::Exception::dontPrint();
        get_hdf5(&NumMoments, &file, (char*)(dirName+"NumMoments").c_str());									
        isRequired = true;
    } catch(H5::Exception& e){}


}
	
template <typename T, unsigned DIM>
void dos<T, DIM>::override_parameters(){
    if(variables.DOS_NumEnergies != -1) NEnergies   = variables.DOS_NumEnergies;
    if(variables.DOS_Name != "")        filename    = variables.DOS_Name;
};

template <typename T, unsigned DIM>
void dos<T, DIM>::fetch_parameters(){
	debug_message("Entered conductivit_dc::read.\n");
	//This function reads all the data from the hdf5 file that's needed to 
  //calculate the dc conductivity
	 

  // Check if the data for the DC conductivity exists
  if(!isRequired){
    std::cout << "Data for DC conductivity does not exist. Exiting.\n";
    exit(1);
  }
  
  // Fetch the number of Chebyshev Moments, temperature and number of points
	get_hdf5(&NumMoments, &file, (char*)(dirName+"NumMoments").c_str());	
	get_hdf5(&NumPoints, &file, (char*)(dirName+"NumPoints").c_str());	


  // Check whether the matrices we're going to retrieve are complex or not
  int complex = systemInfo->isComplex;


  NEnergies = NumPoints;
  

  // Retrieve the Gamma Matrix
  std::string MatrixName = dirName + "MU";
  try{
		debug_message("Filling the MU matrix.\n");
		MU = Eigen::Array<std::complex<T>,-1,-1>::Zero(1, NumMoments);
		
		if(complex)
			get_hdf5(MU.data(), &file, (char*)MatrixName.c_str());
		
		if(!complex){
			Eigen::Array<T,-1,-1> MUReal;
			MUReal = Eigen::Array<T,-1,-1>::Zero(1, NumMoments);
			get_hdf5(MUReal.data(), &file, (char*)MatrixName.c_str());
			
			MU = MUReal.template cast<std::complex<T>>();
		}				

    isPossible = true;
  } catch(H5::Exception& e) {debug_message("DOS: There is no MU matrix.\n");}
	



	file.close();
	debug_message("Left DOS::read.\n");
}

template <typename U, unsigned DIM>
void dos<U, DIM>::calculate(){
  if(!isPossible){
    std::cout << "Cannot retrieve DOS since there is not enough information. Exiting.\n";
    exit(0);
  }

  double Emax = 0.99;
  double Emin = -Emax;
  verbose_message("  Number of points: "); verbose_message(NumPoints); verbose_message("\n");
  verbose_message("  Energy range: ["); verbose_message(Emin); verbose_message(",");
      verbose_message(Emax); verbose_message("]\n");
  verbose_message("  Using kernel: Jackson\n");
  verbose_message("  File name: dos.dat\n");


  Eigen::Matrix<U, -1, 1> energies;
  energies = Eigen::Matrix<U, -1, 1>::LinSpaced(NEnergies, Emin, Emax);



  // First perform the part of the product that only depends on the
  // chebyshev polynomial of the first kind
  Eigen::Array<std::complex<U>, -1, -1> GammaE;
  GammaE = Eigen::Array<std::complex<U>, -1, -1>::Zero(NEnergies, 1);

  U mult = 1.0/systemInfo->energy_scale;

  U factor;
  for(int i = 0; i < NEnergies; i++){
    for(int m = 0; m < NumMoments; m++){
      factor = 1.0/(1.0 + U(m==0));
      GammaE(i) += MU(m)*delta(m,energies(i))*kernel_jackson<U>(m, NumMoments)*factor*mult;
    }
  }

  // Save the density of states to a file and find its maximum value
  U max = -1;
  std::ofstream myfile;
  myfile.open(filename);
  for(int i=0; i < NEnergies; i++){
    myfile  << energies(i)*systemInfo->energy_scale << " " << GammaE.real()(i) << " " << GammaE.imag()(i) << "\n";
    if(GammaE.real()(i) > max)
      max = GammaE.real()(i);
  }

  myfile.close();
  
  
  // Check the limits of the density of states
  U threshold = max*0.01;
  U lowest = 2, highest = -2, a, b;
  bool founda = false, foundb = false;
  for(int i = 0; i < NEnergies; i++){
    a = GammaE.real()(i);
    b = GammaE.real()(NEnergies - i - 1);
    
    if(a > threshold && !founda){
      lowest = energies(i);
      if(i > 0)
        lowest = energies(i-1);
      founda = true;
    }
    
    if(b > threshold && !foundb){
      highest = energies(NEnergies - i -1);
      if(i>0)
        highest = energies(NEnergies-i);
      foundb = true;
    }
    

  }

  systemInfo->EnergyLimitsKnown = true;
  systemInfo->minEnergy = lowest;
  systemInfo->maxEnergy = highest;

};

// Instantiations
template class dos<float, 1u>;
template class dos<float, 2u>;
template class dos<float, 3u>;

template class dos<double, 1u>;
template class dos<double, 2u>;
template class dos<double, 3u>;

template class dos<long double, 1u>;
template class dos<long double, 2u>;
template class dos<long double, 3u>;
