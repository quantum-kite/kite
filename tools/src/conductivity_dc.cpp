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
#include "conductivity_dc.hpp"
#include "functions.hpp"

#include "macros.hpp"


template <typename T, unsigned DIM>
conductivity_dc<T, DIM>::conductivity_dc(system_info<T, DIM>& info, shell_input & vari){
  std::string name = info.filename;
	file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

    isRequired = false; // was this quantity (conductivity_dc) asked for?
    isPossible = false; // do we have all we need to calculate the conductivity?

    int NumPoints = -1;
    double temperature = -1;
    double units = unit_scale;
    std::string filename = "condDC.dat";

  // retrieve the information about the Hamiltonian
  systemInfo = info;

  // retrieve the shell input
  variables = vari;

  // location of the information about the conductivity
  dirName = "/Calculation/conductivity_dc/";
  
  // check whether the conductivity_dc was asked for
  try{
    H5::Exception::dontPrint();
    get_hdf5(&direction, &file, (char*)(dirName+"Direction").c_str());									
    isRequired = true;
  } catch(H5::Exception& e){}
  

}
	

template <typename T, unsigned DIM>
void conductivity_dc<T, DIM>::fetch_parameters(){
	debug_message("Entered conductivit_dc::read.\n");
	//This function reads all the data from the hdf5 file that's needed to 
  //calculate the dc conductivity
	 

  // Check if the data for the DC conductivity exists
  if(!isRequired){
    std::cout << "Data for DC conductivity does not exist. Exiting.\n";
    exit(1);
  }
  
  // Fetch the direction of the conductivity and convert it to a string
  get_hdf5(&direction, &file, (char*)(dirName+"Direction").c_str());									
  std::string dirString = num2str2(direction);

  // Fetch the number of Chebyshev Moments, temperature and number of points
	get_hdf5(&NumMoments, &file, (char*)(dirName+"NumMoments").c_str());	
	get_hdf5(&temperature, &file, (char*)(dirName+"Temperature").c_str());	
	get_hdf5(&NumPoints, &file, (char*)(dirName+"NumPoints").c_str());	


  // Check whether the matrices we're going to retrieve are complex or not
  int complex = systemInfo.isComplex;


  
  scat = 0.0032679;
  beta = 1.0/8.6173303*pow(10,5)/temperature;
  NEnergies = 200;
  if(systemInfo.EnergyLimitsKnown){
    minEnergy = systemInfo.minEnergy;
    maxEnergy = systemInfo.maxEnergy;
  }
  else{
    minEnergy = -0.99;
    maxEnergy = 0.99;
    verbose_message("  - Using default limits for the integration in the energies.\n");
    verbose_message("  - For a more precise evaluation, calculate the density of states as well.\n");
  }

  NFermiEnergies = NumPoints;
  //NFermiEnergies = 100;
  minFermiEnergy = -1.0;
  maxFermiEnergy = 1.0;
  //minEnergy = -0.99;
  //maxEnergy = 0.99;


  // Retrieve the Gamma Matrix
  std::string MatrixName = dirName + "Gamma" + dirString;
  try{
		debug_message("Filling the Gamma matrix.\n");
    //Gamma = Eigen::Array<std::complex<T>,-1,-1, Eigen::RowMajor>::Zero(NumMoments, NumMoments);
    Gamma = Eigen::Array<std::complex<T>,-1,-1>::Zero(NumMoments, NumMoments);
		
		if(complex)
			get_hdf5(Gamma.data(), &file, (char*)MatrixName.c_str());
		
		if(!complex){
			Eigen::Array<T,-1,-1> GammaReal;
      //GammaReal = Eigen::Array<T,-1,-1, Eigen::RowMajor>::Zero(NumMoments, NumMoments);
      GammaReal = Eigen::Array<T,-1,-1>::Zero(NumMoments, NumMoments);
      get_hdf5(GammaReal.data(), &file, (char*)MatrixName.c_str());
			
			Gamma = GammaReal.template cast<std::complex<T>>();
		}				

    isPossible = true;
  } catch(H5::Exception& e) {debug_message("Conductivity DC: There is no Gamma matrix.\n");}
	

	file.close();
	debug_message("Left conductivity_dc::read.\n");
}

template <typename U, unsigned DIM>
void conductivity_dc<U, DIM>::override_parameters(){
    if(variables.CondDC_Temp != -1)         temperature     = variables.CondDC_Temp/systemInfo.energy_scale;
    if(variables.CondDC_NumEnergies != -1)  NEnergies       = variables.CondDC_NumEnergies;
    if(variables.CondDC_Scat != -8888)      scat            = variables.CondDC_Scat/systemInfo.energy_scale;
    if(variables.CondDC_FermiMin != -8888)  minFermiEnergy  = variables.CondDC_FermiMin/systemInfo.energy_scale;
    if(variables.CondDC_FermiMax != -8888)  maxFermiEnergy  = variables.CondDC_FermiMax/systemInfo.energy_scale;
    if(variables.CondDC_NumFermi != -1)     NFermiEnergies  = variables.CondDC_NumFermi;
    if(variables.CondDC_Name != "")         filename        = variables.CondDC_Name;
    
    beta = 1.0/8.6173303*pow(10,5)/temperature;
};

template <typename U, unsigned DIM>
void conductivity_dc<U, DIM>::calculate(){
  if(!isPossible){
    std::cout << "Cannot calculate the DC conductivity because there is not enough information. Exiting.\n";
    exit(0);
  }


  verbose_message("  Energy in rescaled units: [-1,1]\n");
  verbose_message("  Beta (1/kT): "); verbose_message(beta); verbose_message("\n");
  //verbose_message("  Fermi energi (in KPM units): "); verbose_message(e_fermi); verbose_message("\n");
  //verbose_message("  Using kernel for delta function: Jackson\n");
  verbose_message("  Broadening parameter for Green's function: ");
    verbose_message(scat); verbose_message("\n");
  verbose_message("  Number of energies: "); verbose_message(NEnergies); verbose_message("\n");
  verbose_message("  Energy range: ["); verbose_message(minEnergy); verbose_message(",");
    verbose_message(maxEnergy); verbose_message("]\n");
  verbose_message("  Number of Fermi energies: "); verbose_message(NFermiEnergies); verbose_message("\n");
  verbose_message("  Fermi energies range: ["); verbose_message(minFermiEnergy); verbose_message(",");
    verbose_message(maxFermiEnergy); verbose_message("]\n");
  verbose_message("  File name: condDC.dat\n");
  energies = Eigen::Matrix<U, -1, 1>::LinSpaced(NEnergies, minEnergy, maxEnergy);

  // First perform the part of the product that only depends on the
  // chebyshev polynomial of the first kind

  Eigen::Matrix<std::complex<U>, -1, -1, Eigen::ColMajor> greenR;
  greenR = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(NumMoments, NEnergies);
  U factor;
  std::complex<U> complexEnergy;
  for(int i = 0; i < NEnergies; i++){
    complexEnergy = std::complex<U>(energies(i), scat);
    for(int m = 0; m < NumMoments; m++){
      factor = -1.0/(1.0 + U(m==0))/M_PI;
      greenR(m, i) = green(m, 1, complexEnergy).imag()*factor;
    }
  }
  int N_threads;

  Eigen::Array<std::complex<U>, -1, -1> GammaE;
  GammaE = Eigen::Array<std::complex<U>, -1, -1>::Zero(NEnergies, 1);
  omp_set_num_threads(systemInfo.NumThreads);
#pragma omp parallel shared(N_threads) firstprivate(greenR)
  {
#pragma omp master
    {
      N_threads = omp_get_num_threads();
    }
#pragma omp barrier

#pragma omp for schedule(static, 1) nowait
    for(int thread = 0; thread < N_threads; thread++){
      int localMoments = NumMoments/N_threads;

      Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> LocalGamma;
      LocalGamma = Gamma.matrix().block(thread*localMoments, 0, localMoments, NumMoments);

      Eigen::Matrix<std::complex<U>, -1, -1, Eigen::ColMajor> GammaEN;
      GammaEN = LocalGamma*greenR;

      std::complex<U> complexEnergyP, complexEnergyN;



      // Matrices of derivatives of Green's functions
      Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> dgreenA, dgreenR;
      dgreenA = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(NEnergies, localMoments);
      dgreenR = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(NEnergies, localMoments);
      for(int i = 0; i < NEnergies; i++){
        complexEnergyP = std::complex<U>(energies(i), scat);
        complexEnergyN = std::complex<U>(energies(i), -scat);
        for(int m = 0; m < localMoments; m++){
          factor = 1.0/(1.0 + U((m + thread*localMoments)==0));
          dgreenR(i, m) = dgreen<U>(m + thread*localMoments,  1, complexEnergyP)*factor;
          dgreenA(i, m) = dgreen<U>(m + thread*localMoments, -1, complexEnergyN)*factor;
        }
      }


      // Now perform the part of the product that depends on both kinds of polynomials
      Eigen::Array<std::complex<U>, -1, -1> LocalGammaE;
      LocalGammaE = Eigen::Array<std::complex<U>, -1, -1>::Zero(NEnergies, 1);

      U den = systemInfo.num_orbitals*systemInfo.spin_degeneracy/systemInfo.unit_cell_area/units; 
      for(int i = 0; i < NEnergies; i++){
        LocalGammaE(i) += std::complex<U>(dgreenR.row(i)*GammaEN.col(i))*den;
        LocalGammaE(i) += -std::complex<U>(dgreenA.row(i)*GammaEN.col(i).conjugate())*den;
      }


#pragma omp critical
      {
        GammaE += LocalGammaE;
      }
    }
#pragma omp barrier
  }

  Eigen::Matrix<std::complex<U>, -1, 1> condDC;
  condDC = Eigen::Matrix<std::complex<U>, -1, 1>::Zero(NFermiEnergies, 1);

  U energy;
  Eigen::Matrix<U, -1, 1> fermiEnergies;
  fermiEnergies = Eigen::Matrix<U, -1, 1>::LinSpaced(NFermiEnergies, minFermiEnergy, maxFermiEnergy);

  Eigen::Matrix<std::complex<U>, -1, 1> integrand;
  integrand = Eigen::Matrix<std::complex<U>, -1, 1>::Zero(NEnergies, 1);

  U fermi;
  for(int i = 0; i < NFermiEnergies; i++){
    fermi = fermiEnergies(i);
    for(int j = 0; j < NEnergies; j++){
      energy = energies(j);
      integrand(j) = GammaE(j)*fermi_function(energy, fermi, beta);
    }
    condDC(i) = integrate(energies, integrand)*std::complex<U>(0.0,1.0);
  }


  std::ofstream myfile;
  myfile.open(filename);
  for(int i=0; i < NFermiEnergies; i++)
    myfile  << fermiEnergies(i)*systemInfo.energy_scale << " " << condDC.real()(i) << " " << condDC.imag()(i) << "\n";
  
  myfile.close();

};

// Instantiations
template class conductivity_dc<float, 1u>;
template class conductivity_dc<float, 2u>;
template class conductivity_dc<float, 3u>;

template class conductivity_dc<double, 1u>;
template class conductivity_dc<double, 2u>;
template class conductivity_dc<double, 3u>;

template class conductivity_dc<long double, 1u>;
template class conductivity_dc<long double, 2u>;
template class conductivity_dc<long double, 3u>;
