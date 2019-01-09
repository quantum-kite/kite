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
    // Class constructor

    std::string name = info.filename;
	file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);
        H5::Exception::dontPrint();

    isRequired = false; // was this quantity (conductivity_dc) asked for?
    isPossible = false; // do we have all we need to calculate the conductivity?

    units = unit_scale;
    std::cout << "unit scale: " << unit_scale << "\n";

    // retrieve the information about the Hamiltonian
    systemInfo = info;

    // retrieve the shell input
    variables = vari;

    // location of the information about the conductivity
    dirName = "/Calculation/conductivity_dc/";
  
    // check whether the conductivity_dc was asked for

    isRequired = is_required();
    if(isRequired){
    set_default_parameters();   // sets a default set of paramters for the calculation
    fetch_parameters();         // finds all the paramters in the .h5 file
    override_parameters();       // overrides paramters with the ones from the shell input
    set_energy_limits();        // sets the energy limits used in the integration
	file.close();

    printDC();                  // Print all the parameters used
    }
    // needs fetch_parameters() to evaluate isPossible
    if(!isPossible and isRequired){
        std::cout << "Couldn't find the Gamma matrix in the .h5 file. Exiting.\n";
        exit(1);
    }
}
template <typename T, unsigned DIM>
bool conductivity_dc<T, DIM>::is_required(){
    bool result;
    try{
        get_hdf5(&direction, &file, (char*)(dirName+"Direction").c_str());
        result = true;
    } catch(H5::Exception& e){}
    return result;
}
	
template <typename T, unsigned DIM>
void conductivity_dc<T, DIM>::set_default_parameters(){
    // Sets default values for the parameters used in the 
    // calculation of the DC conductivity. These are the parameters
    // that will be overwritten by the config file and the
    // shell input parameters

    NFermiEnergies  = 100;
    minFermiEnergy  = -1.0;
    maxFermiEnergy  = 1.0;
    default_NFermi  = true;
    default_mFermi  = true;
    default_MFermi  = true;

    NEnergies           = 512;      // Number of energies used in the energy integration
    default_NEnergies   = true;

    scat            = 0.015;    // scattering parameter of the Green's functions
    default_scat    = true;

    filename            = "condDC.dat"; // Filename to save the final result
    default_filename    = true;

    temperature     = 0.001;
    beta            = 1.0/8.6173303*pow(10,5)/temperature;
    default_temp    = true;
}

template <typename T, unsigned DIM>
void conductivity_dc<T, DIM>::set_energy_limits(){
    // Attempts to find the energy limits from the information
    // about the density of states stored in the systemInfo class.
    // If it cant, uses default limits.


    minEnergy               = -0.99;    // Minimum energy
    maxEnergy               = 0.99;     // Maximum energy
    default_energy_limits   = true;

    if(systemInfo.EnergyLimitsKnown){
        minEnergy = systemInfo.minEnergy;
        maxEnergy = systemInfo.maxEnergy;
        default_energy_limits = false;
    }
}

template <typename T, unsigned DIM>
void conductivity_dc<T, DIM>::fetch_parameters(){
	debug_message("Entered conductivit_dc::read.\n");
	//This function reads all the data from the hdf5 file that's needed to 
    //calculate the dc conductivity
	 
    // Fetch the direction of the conductivity and convert it to a string
    get_hdf5(&direction, &file, (char*)(dirName+"Direction").c_str());
    std::string dirString = num2str2(direction);

    // Fetch the number of Chebyshev Moments
	get_hdf5(&NumMoments, &file, (char*)(dirName+"NumMoments").c_str());	

    // Fetch the temperature from the .h5 file
	get_hdf5(&temperature, &file, (char*)(dirName+"Temperature").c_str());	
    beta = 1.0/8.6173303*pow(10,5)/temperature;
    default_temp = false;

    // Fetch the number of Fermi energies from the .h5 file
	get_hdf5(&NFermiEnergies, &file, (char*)(dirName+"NumPoints").c_str());	
    default_NFermi = false;


    // Check whether the matrices we're going to retrieve are complex or not
    int complex = systemInfo.isComplex;


  

    // Retrieve the Gamma Matrix
    std::string MatrixName = dirName + "Gamma" + dirString;
    try{
        debug_message("Filling the Gamma matrix.\n");
        Gamma = Eigen::Array<std::complex<T>,-1,-1>::Zero(NumMoments, NumMoments);
		
		if(complex)
			get_hdf5(Gamma.data(), &file, (char*)MatrixName.c_str());
		
		if(!complex){
			Eigen::Array<T,-1,-1> GammaReal;
            GammaReal = Eigen::Array<T,-1,-1>::Zero(NumMoments, NumMoments);
            get_hdf5(GammaReal.data(), &file, (char*)MatrixName.c_str());
			
			Gamma = GammaReal.template cast<std::complex<T>>();
		}				

        isPossible = true;
    } catch(H5::Exception& e) {
        debug_message("Conductivity DC: There is no Gamma matrix.\n");
    }
	

	debug_message("Left conductivity_dc::read.\n");
}

template <typename U, unsigned DIM>
void conductivity_dc<U, DIM>::override_parameters(){
    // overrides the current parameters with the ones from the shell input

    if(variables.CondDC_Temp != -1){
        temperature     = variables.CondDC_Temp/systemInfo.energy_scale;
        beta            = 1.0/8.6173303*pow(10,5)/temperature;
        default_temp    = false;
    }

    if(variables.CondDC_NumEnergies != -1){
        NEnergies = variables.CondDC_NumEnergies;
        default_NEnergies = false;
    }

    if(variables.CondDC_Scat != -8888){
        scat            = variables.CondDC_Scat/systemInfo.energy_scale;
        default_scat    = false;
    }

    if(variables.CondDC_FermiMin != -8888){
        minFermiEnergy  = variables.CondDC_FermiMin/systemInfo.energy_scale;
        default_mFermi   = false;
    }

    if(variables.CondDC_FermiMax != -8888){  
        maxFermiEnergy  = variables.CondDC_FermiMax/systemInfo.energy_scale;
        default_MFermi   = false;
    }

    if(variables.CondDC_NumFermi != -1){
        NFermiEnergies  = variables.CondDC_NumFermi;
        default_NFermi   = false;
    }

    if(variables.CondDC_Name != ""){
        filename            = variables.CondDC_Name;
        default_filename    = false;
    }
}


template <typename U, unsigned DIM>
void conductivity_dc<U, DIM>::printDC(){
    double scale = systemInfo.energy_scale;
    std::string energy_range = "[" + std::to_string(minEnergy*scale) + ", " + std::to_string(maxEnergy*scale) + "]";

    // Prints all the information about the parameters
    std::cout << "The DC conductivity will be calculated with these parameters: (eV, Kelvin)\n"
        "   Temperature: "             << temperature*scale        << ((default_temp)?         " (default)":"") << "\n"
        "   Broadening: "              << scat*scale               << ((default_scat)?         " (default)":"") << "\n"
        "   Max Fermi energy: "        << maxFermiEnergy*scale     << ((default_MFermi)?       " (default)":"") << "\n"
        "   Min Fermi energy: "        << minFermiEnergy*scale     << ((default_mFermi)?       " (default)":"") << "\n"
        "   Number Fermi energies: "   << NFermiEnergies           << ((default_NFermi)?       " (default)":"") << "\n"
        "   Filename: "                << filename                 << ((default_filename)?     " (default)":"") << "\n"
        "   Integration range: "       << energy_range             << ((default_energy_limits)?" (default)":" (Estimated from DoS)") << "\n"
        "   Num integration points: "  << NEnergies                << ((default_NEnergies)?    " (default)":"") << "\n"; 
}

template <typename U, unsigned DIM>
void conductivity_dc<U, DIM>::calculate(){


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

  std::cout << "normGreenR: " << greenR.norm() << "\n";

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

      std::cout << "LocalGamma: " << omp_get_thread_num() << " " << LocalGamma.norm() << "\n";
      Eigen::Matrix<std::complex<U>, -1, -1, Eigen::ColMajor> GammaEN;
      GammaEN = LocalGamma*greenR;

      std::cout << "GammaEN: " << GammaEN.norm() << "\n";

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

      std::cout << "R:" << dgreenR.norm() << "\n";
      std::cout << "A:" << dgreenA.norm() << "\n";


      // Now perform the part of the product that depends on both kinds of polynomials
      Eigen::Array<std::complex<U>, -1, -1> LocalGammaE;
      LocalGammaE = Eigen::Array<std::complex<U>, -1, -1>::Zero(NEnergies, 1);

      U den = systemInfo.num_orbitals*systemInfo.spin_degeneracy/systemInfo.unit_cell_area/units; 
      std::cout << "units:" << units << "\n";
      std::cout << "den: " << den<<"\n";
      for(int i = 0; i < NEnergies; i++){
        LocalGammaE(i) += std::complex<U>(dgreenR.row(i)*GammaEN.col(i))*den;
        LocalGammaE(i) += -std::complex<U>(dgreenA.row(i)*GammaEN.col(i).conjugate())*den;
      }

      std::cout << "LocalGammaE norm: " << LocalGammaE.matrix().norm() << "\n";

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
