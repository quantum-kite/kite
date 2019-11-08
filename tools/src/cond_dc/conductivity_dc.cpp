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
#include "../ComplexTraits.hpp"
#include "../myHDF5.hpp"

#include "../parse_input.hpp"
#include "../systemInfo.hpp"
#include "conductivity_dc.hpp"
#include "../functions.hpp"

#include "../macros.hpp"

template <typename T, unsigned DIM>
conductivity_dc<T, DIM>::conductivity_dc(system_info<T, DIM>& info, shell_input & vari){
    /* Class constructor
     * Finds all the information needed to compute the DC conductivity
     * but does not compute it. To calculate the DC conductivity, the method 
     * calculate() needs to be called. If there is not enough data to compute
     * the DC conductivity, the code will exit with an error.
     */

    H5::Exception::dontPrint();
    units = unit_scale;     // units of the DC conductivity
    systemInfo = info;      // retrieve the information about the Hamiltonian
    variables = vari;       // retrieve the shell input

    isPossible = false;         // do we have all we need to calculate the conductivity?
    isRequired = is_required() && variables.CondDC_is_required; // was this quantity (conductivity_dc) asked for?


    // If the DC conductivity was requested, all the necessary parameters for its
    // computation will be set and printed out to std::cout
    if(isRequired){
        set_default_parameters();        // sets a default set of paramters for the calculation
        isPossible = fetch_parameters(); // finds all the paramters in the .h5 file
        override_parameters();           // overrides parameters with the ones from the shell input
        set_energy_limits();             // sets the energy limits used in the integration


        if(isPossible){
          printDC();                  // Print all the parameters used
          calculate2();
        } else {
          std::cout << "ERROR. The DC conductivity was requested but the data "
              "needed for its computation was not found in the input .h5 file. "
              "Make sure KITEx has processed the file first. Exiting.";
          exit(1);
        }
    }
}

template <typename T, unsigned DIM>
bool conductivity_dc<T, DIM>::is_required(){
    // Checks whether the DC conductivity has been requested
    // by analysing the .h5 config file. If it has been requested, 
    // some fields have to exist, such as "Direction"


    // Make sure the config filename has been initialized
    std::string name = systemInfo.filename.c_str();
    if(name == ""){
        std::cout << "ERROR: Filename uninitialized. Exiting.\n";
        exit(1);
    }

    // location of the information about the conductivity
    char dirName[] = "/Calculation/conductivity_dc/Direction";
	H5::H5File file = H5::H5File(name, H5F_ACC_RDONLY);

    // Check if this dataset exists
    bool result = false;
    try{
        get_hdf5(&direction, &file, dirName);
        result = true;
    } catch(H5::Exception& e){}


    file.close();
    return result;
}
	
template <typename T, unsigned DIM>
void conductivity_dc<T, DIM>::set_default_parameters(){
    // Sets default values for the parameters used in the 
    // calculation of the DC conductivity. These are the parameters
    // that will be overwritten by the config file and the
    // shell input parameters. 

    double scale = systemInfo.energy_scale;
    double shift = systemInfo.energy_shift;



    NFermiEnergies      = 100;
    minFermiEnergy      = (-1.0 - shift)/scale;        // -1eV in KPM reduced units [-1,1]
    maxFermiEnergy      = ( 1.0 - shift)/scale;        //  1eV in KPM reduced units [-1,1]
    default_NFermi      = true;
    default_mFermi      = true;
    default_MFermi      = true;

    default_NumMoments  = true;

    // integrate means to use the full energy range in the integration
    full_range          = 1;
    default_full_range  = true;
    NumThreads          = systemInfo.NumThreads;
    default_NumThreads  = true;

    NEnergies           = 512;              // Number of energies used in the energy integration
    default_NEnergies   = true;

    deltascat           = 0.01/scale;       // scattering parameter in the delta function
    scat                = 0.01/scale;       // scattering parameter of 10meV in
    default_scat        = true;             // the Green's functions in KPM reduced units
    default_deltascat   = true;             // the Green's functions in KPM reduced units

    filename            = "condDC.dat";     // Filename to save the final result
    default_filename    = true;

    temperature         = 0.001/scale;      // Temperature of 1mK in KPM reduced units
    beta                = 1.0/8.6173303*pow(10,5)/temperature;
    default_temp        = true;
}

template <typename T, unsigned DIM>
void conductivity_dc<T, DIM>::set_energy_limits(){
    // Attempts to find the energy limits from the information
    // about the density of states stored in the systemInfo class.
    // If it cant, uses default limits.


    minEnergy               = -0.99;    // Minimum energy
    maxEnergy               = 0.99;     // Maximum energy
    default_energy_limits   = true;

    // Choose whether or not to use the limits of integration as
    // computed from the density of states. 
    if(systemInfo.EnergyLimitsKnown and !full_range){
        minEnergy = systemInfo.minEnergy;
        maxEnergy = systemInfo.maxEnergy;
        default_energy_limits = false;
    }
}

template <typename T, unsigned DIM>
bool conductivity_dc<T, DIM>::fetch_parameters(){
	debug_message("Entered conductivit_dc::read.\n");
	//This function reads all the data from the hdf5 file that's needed to 
    //calculate the dc conductivity
	 
	H5::H5File file;
  std::string dirName = "/Calculation/conductivity_dc/";  // location of the information about the conductivity
	file = H5::H5File(systemInfo.filename.c_str(), H5F_ACC_RDONLY);

  // Fetch the direction of the conductivity and convert it to a string
  get_hdf5(&direction, &file, (char*)(dirName+"Direction").c_str());
  std::string dirString = num2str2(direction);

  // Fetch the number of Chebyshev Moments
	get_hdf5(&MaxMoments, &file, (char*)(dirName+"NumMoments").c_str());	

  // Fetch the temperature from the .h5 file
  // The temperature is already in KPM reduced units
	get_hdf5(&temperature, &file, (char*)(dirName+"Temperature").c_str());	
  beta = 1.0/8.6173303*pow(10,5)/temperature;
  default_temp = false;

  // Fetch the number of Fermi energies from the .h5 file
	get_hdf5(&NFermiEnergies, &file, (char*)(dirName+"NumPoints").c_str());	
  default_NFermi = false;

  NumMoments = MaxMoments;
  default_NumMoments = true;


  // Check whether the matrices we're going to retrieve are complex or not
  int complex = systemInfo.isComplex;

  // Retrieve the Gamma Matrix
  std::string MatrixName = dirName + "Gamma" + dirString;
  bool possible = false;
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

    possible = true;
  } catch(H5::Exception& e) {
      debug_message("Conductivity DC: There is no Gamma matrix.\n");
  }
	

  file.close();
	debug_message("Left conductivity_dc::read.\n");
  return possible;
}

template <typename U, unsigned DIM>
void conductivity_dc<U, DIM>::override_parameters(){
    // Overrides the current parameters with the ones from the shell input.
    // These parameters are in eV or Kelvin, so they must scaled down
    // to the KPM units. This includes the temperature

    double scale = systemInfo.energy_scale;
    double shift = systemInfo.energy_shift;

    if(variables.CondDC_Temp != -1){
        temperature     = variables.CondDC_Temp/scale;
        beta            = 1.0/8.6173303*pow(10,5)/temperature;
        default_temp    = false;
    }

    if(variables.CondDC_NumEnergies != -1){
        NEnergies = variables.CondDC_NumEnergies;
        default_NEnergies = false;
    }

    if(variables.CondDC_nthreads != -1){
        NumThreads = variables.CondDC_nthreads;
        default_NumThreads = false;

        if(NumThreads < 1){
          std::cout << "NumThreads cannot be smaller than 1. Aborting.\n";
          assert(NumThreads > 0);
          exit(1);
        }
    }

    if(variables.CondDC_NumMoments != -1){
        NumMoments = variables.CondDC_NumMoments;
        default_NumMoments = false;

        if(NumMoments > MaxMoments){
          std::cout << "NumMoments cannot be larger than the number of Chebyshev ";
          std::cout << "moments computed with KITEx. Aborting.\n";
            
          assert(NumMoments <= MaxMoments);
          exit(1);
        }
    }

    // integrate = true means to use the full energy range in the integration
    if(variables.CondDC_integrate != -1){
        full_range = variables.CondDC_integrate;
        default_full_range = false;
    }


    if(variables.CondDC_Scat != -8888){
        scat            = variables.CondDC_Scat/scale;
        default_scat    = false;
    }

    if(variables.CondDC_deltaScat != -8888){
        deltascat         = variables.CondDC_deltaScat/scale;
        default_deltascat = false;
    } else {
        deltascat = scat;
    }


    if(variables.CondDC_FermiMin != -8888){
        minFermiEnergy  = (variables.CondDC_FermiMin - shift)/scale;
        default_mFermi   = false;
    }

    if(variables.CondDC_FermiMax != -8888){  
        maxFermiEnergy  = (variables.CondDC_FermiMax - shift)/scale;
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
    
  debug_message("Padding Gamma");
  //int N_threads = systemInfo.NumThreads;
  int rest = NumMoments % NumThreads;
  Moments_D = NumMoments - rest + (rest != 0)*NumThreads;
  Moments_G = NumMoments;
  //std::cout << "rest: " << rest << "\n";
  //std::cout << "N_threads: " << N_threads << "\n";
  //std::cout << "D: " << Moments_D << "G: " << Moments_G << "\n";

  Gamma_Padded = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(Moments_D, Moments_G);
  Gamma_Padded.block(0,0,NumMoments, NumMoments) = Gamma.block(0,0,NumMoments, NumMoments);

  //std::cout << "Gamma_Padded: " << Gamma_Padded << "\n";
}


template <typename U, unsigned DIM>
void conductivity_dc<U, DIM>::printDC(){
    double scale = systemInfo.energy_scale;
    double shift = systemInfo.energy_shift;
    std::string energy_range = "[" + std::to_string(minEnergy*scale + shift) + ", " + std::to_string(maxEnergy*scale + shift) + "]";

    // Prints all the information about the parameters
    std::cout << "The DC conductivity will be calculated with these parameters: (eV, Kelvin)\n"
        "   Temperature: "             << temperature*scale             << ((default_temp)?         " (default)":"") << "\n"
        "   Broadening: "              << scat*scale                    << ((default_scat)?         " (default)":"") << "\n"
        "   Delta broadening: "        << deltascat*scale               << ((default_deltascat)?    " (default)":"") << "\n"
        "   Max Fermi energy: "        << maxFermiEnergy*scale + shift  << ((default_MFermi)?       " (default)":"") << "\n"
        "   Min Fermi energy: "        << minFermiEnergy*scale + shift  << ((default_mFermi)?       " (default)":"") << "\n"
        "   Number Fermi energies: "   << NFermiEnergies                << ((default_NFermi)?       " (default)":"") << "\n"
        "   Filename: "                << filename                      << ((default_filename)?     " (default)":"") << "\n"
        "   Integration range: "       << energy_range                  << ((default_energy_limits)?" (default)":" (Estimated from DoS)") << "\n"
        "   Num integration points: "  << NEnergies                     << ((default_NEnergies)?    " (default)":"") << "\n"
        "   Num Chebychev moments: "   << NumMoments                    << ((default_NumMoments)?   " (default)":"") << "\n"
        "   Num threads: "             << NumThreads                    << ((default_NumThreads)?   " (default)":"") << "\n"; 
}



template <typename U, unsigned DIM>
void conductivity_dc<U, DIM>::calculate2(){
  energies = Eigen::Matrix<U, -1, 1>::LinSpaced(NEnergies, minEnergy, maxEnergy);
  fermiEnergies = Eigen::Matrix<U, -1, 1>::LinSpaced(NFermiEnergies, minFermiEnergy, maxFermiEnergy);


  // Imaginary part of the Green's function: Dirac delta (greenR)
  // Derivative of the Green's function: dgreenR
  Eigen::Matrix<std::complex<U>, -1, -1, Eigen::ColMajor> greenR;
  Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> dgreenR;

  // Fill the matrices that are going to be used in the multiplication
  // This is an operation of order ~ (NG + ND) * NE
  // and uses (ND + NG) * NE memory (complex U)
  greenR  = fill_delta();
  dgreenR = fill_dgreenR();


  // Product of all the matrices:
  // dgreenR * Gamma_Padded * greenR
  // This is an operation of order N^2 * NE and is parallelized
  // It uses T * NG * NE   +  NE * ND   +   2 * NE
  Eigen::Matrix<std::complex<U>, -1, -1> GammaE;
  GammaE = triple_product(greenR, dgreenR);


  // integrate over the whole energy range for each Fermi energy
  Eigen::Matrix<std::complex<U>, -1, 1> condDC;
  U den = -systemInfo.num_orbitals*systemInfo.spin_degeneracy/systemInfo.unit_cell_area/units; 
  condDC = calc_cond(GammaE)*den;

  // save to a file
  save_to_file(condDC);
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
