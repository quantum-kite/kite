
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
#include <string>
#include <vector>
#include <omp.h>

#include "../ComplexTraits.hpp"
#include "H5Cpp.h"
#include "../myHDF5.hpp"

#include "../parse_input.hpp"
#include "../systemInfo.hpp"
#include "conductivity_2order.hpp"
#include "../functions.hpp"

#include "../macros.hpp"


template <typename T, unsigned DIM>
int conductivity_nonlinear<T, DIM>::fetch_gamma0(){
  // Retrieve the Gamma0 Matrix. This is the 1D Gamma matrix defined as
  // Gamma0 = Tr[v^abc T_n]. This matrix is not needed for the special case

  int hasGamma0 = 0;
  std::string MatrixName = dirName + "Gamma0" + dirString;
  try{
    debug_message("Filling the Gamma0 matrix.\n");
    Gamma0 = Eigen::Array<std::complex<T>,-1,-1>::Zero(1, NumMoments);
      
    if(complex)
      get_hdf5(Gamma0.data(), &file, (char*)MatrixName.c_str());
      
    if(!complex){
      Eigen::Array<T,-1,-1> Gamma0Real;
      Gamma0Real = Eigen::Array<T,-1,-1>::Zero(1, NumMoments);
      get_hdf5(Gamma0Real.data(), &file, (char*)MatrixName.c_str());
      Gamma0 = Gamma0Real.template cast<std::complex<T>>();
    }				

    hasGamma0 = 1;
  } catch(H5::Exception& e) {
    debug_message("Conductivity nonlinear: There is no Gamma0 matrix.\n");
  }
  return hasGamma0;
}



template <typename T, unsigned DIM>
int conductivity_nonlinear<T, DIM>::fetch_gamma1(){
  // Retrieve the Gamma1 Matrix. This is the 2D Gamma matrix defined as
  // Gamma1 = Tr[v^a Tn v^bc Tm]

  int hasGamma1 = 0;
  std::string MatrixName = dirName + "Gamma1" + dirString;
  try{
    debug_message("Filling the Gamma1 matrix.\n");
    Gamma1 = Eigen::Array<std::complex<T>,-1,-1>::Zero(NumMoments, NumMoments);

    if(complex)
      get_hdf5(Gamma1.data(), &file, (char*)MatrixName.c_str());

    if(!complex){
      Eigen::Array<T,-1,-1> Gamma1Real;
      Gamma1Real = Eigen::Array<T,-1,-1>::Zero(NumMoments, NumMoments);
      get_hdf5(Gamma1Real.data(), &file, (char*)MatrixName.c_str());
      Gamma1 = Gamma1Real.template cast<std::complex<T>>();
    }				

    hasGamma1 = 1;
  } catch(H5::Exception& e) {
    debug_message("Conductivity optical: There is no Gamma1 matrix.\n");
  }

  return hasGamma1;
}


template <typename T, unsigned DIM>
int conductivity_nonlinear<T, DIM>::fetch_gamma2(){
  // Retrieve the Gamma2 Matrix
  int hasGamma2 = 0;
  std::string MatrixName = dirName + "Gamma2" + dirString;
  try{
    debug_message("Filling the Gamma2 matrix.\n");
    Gamma2 = Eigen::Array<std::complex<T>,-1,-1>::Zero(NumMoments, NumMoments);
		
    if(complex)
      get_hdf5(Gamma2.data(), &file, (char*)MatrixName.c_str());
		
    if(!complex){
      Eigen::Array<T,-1,-1> Gamma2Real;
      Gamma2Real = Eigen::Array<T,-1,-1>::Zero(NumMoments, NumMoments);
      get_hdf5(Gamma2Real.data(), &file, (char*)MatrixName.c_str());
      Gamma2 = Gamma2Real.template cast<std::complex<T>>();
    }				

    hasGamma2 = 1;
  } catch(H5::Exception& e) {
    debug_message("Conductivity optical: There is no Gamma2 matrix.\n");
  }

  return hasGamma2;
}

  
template <typename T, unsigned DIM>
int conductivity_nonlinear<T, DIM>::fetch_gamma3(){
  // Retrieve the Gamma3 Matrix. This is the biggest matrix and doesn't need to be
  // calculated when we want hBN because it is identically zero.

  int hasGamma3 = 0;
  std::string MatrixName = dirName + "Gamma3" + dirString;
  try{
    debug_message("Filling the Gamma3 matrix.\n");
    Gamma3 = Eigen::Array<std::complex<T>,-1,-1>::Zero(1, NumMoments*NumMoments*NumMoments);
      
    if(complex)
      get_hdf5(Gamma3.data(), &file, (char*)MatrixName.c_str());
     
    if(!complex){
      Eigen::Array<T,-1,-1> Gamma3Real;
      Gamma3Real = Eigen::Array<T,-1,-1>::Zero(1, NumMoments*NumMoments*NumMoments);
      get_hdf5(Gamma3Real.data(), &file, (char*)MatrixName.c_str());
        
      Gamma3 = Gamma3Real.template cast<std::complex<T>>();
    }				

    hasGamma3 = 1;
  } catch(H5::Exception& e) {
    debug_message("Conductivity DC: There is no Gamma3 matrix.\n");
  }

  return hasGamma3;
}

template <typename T, unsigned DIM>
conductivity_nonlinear<T, DIM>::conductivity_nonlinear(system_info<T, DIM>& info, shell_input & vari){
  // Constructor of the conductivity_nonlinear class. This function simply checks
  // whether the conductivity needs to be calculated or not


    isRequired = false; // was this quantity (conductivity_dc) asked for?

    // KPM parameters
    NumPoints = -1;
    special = 0;

  std::string name = info.filename;                           // name of the hdf5 file
  file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);
  systemInfo = info;                                          // retrieve the information about the Hamiltonian
  variables = vari;                                           // retrieve the shell input
  dirName = "/Calculation/conductivity_optical_nonlinear/";   // location of the information about the conductivity


  // check whether the conductivity_nonlinear was asked for
  try{
    H5::Exception::dontPrint();
    get_hdf5(&direction, &file, (char*)(dirName+"Direction").c_str());									
    isRequired = true;
  } catch(H5::Exception& e){}
  
  if(isRequired and variables.CondOpt2_is_required){
    set_default_parameters();
    isPossible = fetch_parameters(); // do we have all we need to calculate the conductivity?
    override_parameters();
    
    if(isPossible){
      printOpt2();
      calculate();
    } else {
      std::cout << "ERROR. The nonlinear optical conductivity was requested but the data "
        "needed for its computation was not found in the input .h5 file. "
        "Make sure KITEx has processed the file first. Exiting.";
      exit(1);
    }
  }
}
	
template <typename T, unsigned DIM>
void conductivity_nonlinear<T, DIM>::override_parameters(){
  double scale = systemInfo.energy_scale;
  double shift = systemInfo.energy_shift;

    if(variables.CondOpt2_Temp != -8888)     temperature = variables.CondOpt2_Temp/scale;
    if(variables.CondOpt2_NumEnergies != -1) N_energies  = variables.CondOpt2_NumEnergies;
    if(variables.CondOpt2_ratio != 8.8888)   ratio       = variables.CondOpt2_ratio;
    if(variables.CondOpt2_print_all != -2)   print_all   = variables.CondOpt2_print_all;
    if(variables.CondOpt2_FreqMax != -8888)  maxFreq     = variables.CondOpt2_FreqMax/scale;
    if(variables.CondOpt2_FreqMin != -8888)  minFreq     = variables.CondOpt2_FreqMin/scale;
    if(variables.CondOpt2_NumFreq != -1)     N_omegas    = variables.CondOpt2_NumFreq;
    if(variables.CondOpt2_Fermi != -8888)    e_fermi     = (variables.CondOpt2_Fermi - shift)/scale;
    if(variables.CondOpt2_Scat != -8888)     scat        = variables.CondOpt2_Scat/scale;
    if(variables.CondOpt2_Name != "")        filename    = variables.CondOpt2_Name;
    beta = 1.0/8.6173303*pow(10,5)/temperature;
}

template <typename T, unsigned DIM>
void conductivity_nonlinear<T, DIM>::set_default_parameters(){
    
  double scale = systemInfo.energy_scale;
  double shift = systemInfo.energy_shift;
  temperature = 0.01/scale;
  beta        = 1.0/8.6173303*pow(10,5)/temperature;
  default_temperature = true;


  print_all   = 0;
  N_energies  = 512; 
  default_NEnergies = true;

  lim         = 0.995;

  maxFreq     = 7.0/systemInfo.energy_scale; 
  minFreq     = 0.0/systemInfo.energy_scale;
  N_omegas    = 512;
  ratio       = 1.0;

  default_minfreq = true;
  default_maxfreq = true;
  default_Nfreqs = true;
  default_ratio = true;

  e_fermi     = (0.0 - shift)/scale;
  scat        = 0.0166/scale;
  filename    = "nonlinear_cond.dat";  
  
  default_efermi = true;
  default_scat = true;
  default_filename = true;
}

template <typename T, unsigned DIM>
void conductivity_nonlinear<T, DIM>::printOpt2(){
  double scale = systemInfo.energy_scale;
  double shift = systemInfo.energy_shift;
  std::cout << "The second-order optical conductivity will be calculated with these parameters: (eV)\n"
    "   Temperature: "            << temperature*scale      << ((default_temperature)?  " (default)":"") << "\n"
    "   Broadening: "             << scat*scale             << ((default_scat)?         " (default)":"") << "\n"
    "   Fermi energy: "           << e_fermi*scale + shift  << ((default_efermi)?       " (default)":"") << "\n"
    "   Min frequency: "          << minFreq*scale          << ((default_minfreq)?     " (default)":"") << "\n"
    "   Max frequency: "          << maxFreq*scale          << ((default_maxfreq)?     " (default)":"") << "\n"
    "   Number of frequencies: "  << N_omegas               << ((default_Nfreqs)?       " (default)":"") << "\n"
    "   Ratio: "                  << ratio                  << ((default_ratio)?        " (default)":"") << "\n"
    "   Num integration points: " << N_energies             << ((default_NEnergies)?    " (default)":"") << "\n"
    //"   Integration range: "       << energy_range             << ((default_energy_limits)?" (default)":" (Estimated from DoS)") << "\n"
    "   Kernel for Dirac deltas: "<< "jackson"              << " (default)\n"
    "   Filename: "               << filename               << ((default_filename)?     " (default)":"") << "\n";
}

template <typename T, unsigned DIM>
bool conductivity_nonlinear<T, DIM>::fetch_parameters(){
  debug_message("Entered conductivity_nonlinear::read.\n");
  // This function reads all the relevant
  // information from the hdf5 configuration file and uses it to evaluate the parameters
  // needed to calculate the nonlinear conductivity
	 

  // This function should not run if the conductivity is not needed. If, for some reason
  // it is run anyway, the user should be notified that there is not enough data to
  // calculate the conductivity.
  if(!isRequired){
    std::cout << "Data for nonlinear conductivity does not exist. Exiting.\n";
    exit(1);
  }

  imaginary = std::complex<T>(0.0, 1.0);
  
  // Fetch the direction of the conductivity and convert it to a string
  get_hdf5(&direction, &file, (char*)(dirName+"Direction").c_str());
  dirString = num2str3(direction);

  // Fetch the number of Chebyshev Moments, temperature and number of points
  get_hdf5(&NumMoments, &file, (char*)(dirName+"NumMoments").c_str());	
  get_hdf5(&temperature, &file, (char*)(dirName+"Temperature").c_str());	
  get_hdf5(&NumPoints, &file, (char*)(dirName+"NumPoints").c_str());	
  get_hdf5(&special, &file, (char*)(dirName+"Special").c_str());	
  beta = 1.0/8.6173303*pow(10,5)/temperature;   // 1/kT, where k is the Boltzmann constant in eV/K
	
  // Frequency parameters needed to run the simulation
  N_omegas = NumPoints;

  // Check whether the matrices we're going to retrieve are complex or not
  complex = systemInfo.isComplex;

  bool hasGamma0 = false;
  bool hasGamma1 = false;
  bool hasGamma2 = false;
  bool hasGamma3 = false;

  // Fetch the Gamma matrices from the hdf file and return the success value
  if(not special) hasGamma0 = fetch_gamma0(); 
  hasGamma1 = fetch_gamma1();
  hasGamma2 = fetch_gamma2(); 
  if(not special) hasGamma3 = fetch_gamma3();

  // check if we have all the objects that we need
  bool possible = false;
  if(special)     possible = hasGamma1 and hasGamma2;
  if(not special) possible = hasGamma0 and hasGamma1 and hasGamma2 and hasGamma3;

  file.close();
  debug_message("Left conductivity_nonlinear::read.\n");
  return possible;
}


template <typename U, unsigned DIM>
void conductivity_nonlinear<U, DIM>::calculate_photo(){
  Eigen::Matrix<std::complex<U>, -1, -1> omega_energies0, omega_energies1, omega_energies2, omega_energies3, omega_energies4;
  omega_energies0 = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  omega_energies1 = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  omega_energies2 = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  omega_energies3 = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  omega_energies4 = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);

  Eigen::Matrix<std::complex<U>,1,-1> cond;
  Eigen::Matrix<std::complex<U>,1,-1> cond0, cond1, cond2, cond3, cond4;

  cond0 = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);
  cond1 = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);
  cond2 = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);
  cond3 = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);
  cond4 = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);

  cond     = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);

  // Contraction of the Gamma matrices with the delta functions and Green's functions
  omega_energies0 += Gamma0contract();
  omega_energies1 += 0.5*Gamma1contractAandR(); 
  omega_energies2 += Gamma2contractAandR(); 

  if(not special){
    omega_energies4 += Gamma3Contract_RA(); 
    omega_energies3 += Gamma3Contract_RRandAAblocks();
  }

  U freq;
  for(int w = 0; w < N_omegas; w++){
    freq = frequencies(w);  

    // Energy integration
    cond0(w) = integrate(energies, Eigen::Matrix<std::complex<U>,-1,1>(omega_energies0.col(w)));
    cond1(w) = integrate(energies, Eigen::Matrix<std::complex<U>,-1,1>(omega_energies1.col(w)));
    cond2(w) = integrate(energies, Eigen::Matrix<std::complex<U>,-1,1>(omega_energies2.col(w)));
    cond3(w) = integrate(energies, Eigen::Matrix<std::complex<U>,-1,1>(omega_energies3.col(w)));
    cond4(w) = integrate(energies, Eigen::Matrix<std::complex<U>,-1,1>(omega_energies4.col(w)));

    // Divide by the frequencies
    cond0(w) /= -scat*scat - freq*freq; 
    cond1(w) /= -scat*scat - freq*freq; 
    cond2(w) /= -scat*scat - freq*freq; 
    cond3(w) /= -scat*scat - freq*freq; 
    cond4(w) /= -scat*scat - freq*freq; 
  }

  U compat_factor = 2.0;//*systemInfo.energy_scale*systemInfo.energy_scale;
  std::complex<U> factor = imaginary*U(systemInfo.num_orbitals*
      systemInfo.spin_degeneracy/systemInfo.unit_cell_area/systemInfo.energy_scale)*compat_factor;

  cond0 *= factor;
  cond1 *= factor;
  cond2 *= factor;
  cond3 *= factor;
  cond4 *= factor;
  cond = cond1 + cond2 + cond3 + cond4;

  std::ofstream myfile;
  std::ofstream myfile0, myfile1, myfile2, myfile3, myfile4;

  myfile.open(filename);

  if(print_all){
    myfile0.open(filename + "0");
    myfile1.open(filename + "1");
    myfile2.open(filename + "2");
    myfile3.open(filename + "RR_AA");
    myfile4.open(filename + "RA");
  }


  for(int i=0; i < N_omegas; i++){
    freq = std::real(frequencies(i))*systemInfo.energy_scale;
    myfile  << freq << " " << cond.real()(i) << " " << cond.imag()(i) << "\n";
    if(print_all){
      myfile0 << freq << " " << cond0.real()(i) << " " << cond0.imag()(i) << "\n";
      myfile1 << freq << " " << cond1.real()(i) << " " << cond1.imag()(i) << "\n";
      myfile2 << freq << " " << cond2.real()(i) << " " << cond2.imag()(i) << "\n";
      myfile3 << freq << " " << cond3.real()(i) << " " << cond3.imag()(i) << "\n";
      myfile4 << freq << " " << cond4.real()(i) << " " << cond4.imag()(i) << "\n";
    } 
  }
  myfile.close();
  myfile0.close();
  myfile1.close();
  myfile2.close();
  myfile3.close();
  myfile4.close();
}

template <typename U, unsigned DIM>
void conductivity_nonlinear<U, DIM>::calculate_general(){

  Eigen::Matrix<std::complex<U>, -1, -1> omega_energies3shg1, omega_energies3shg2, omega_energies3shg3;
  Eigen::Matrix<std::complex<U>, -1, -1> omega_energies2shg, omega_energies1shg;

  omega_energies1shg  = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  omega_energies2shg  = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  omega_energies3shg1 = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  omega_energies3shg2 = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  omega_energies3shg3 = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);


  Eigen::Matrix<std::complex<U>,1,-1> cond_shg;
  Eigen::Matrix<std::complex<U>,1,-1> cond3shg1, cond3shg2, cond3shg3, cond1shg, cond2shg;

  cond1shg  = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);
  cond2shg  = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);
  cond3shg1 = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);
  cond3shg2 = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);
  cond3shg3 = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);

  cond_shg = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);

  // Contraction of the Gamma matrices with the delta functions and Green's functions

  omega_energies1shg += 0.5*Gamma1shgcontractAandR();
  omega_energies2shg += Gamma2shgcontractAandR();

  if(not special){
    omega_energies3shg1 += Gamma3shgContract_RA();
    omega_energies3shg2 += Gamma3shgContract_RR();
    omega_energies3shg3 += Gamma3shgContract_AA();
  }

  U freq;
  U w1, w2;
  for(int w = 0; w < N_omegas; w++){
    freq = frequencies(w);  
    w1 = frequencies2(w,0);
    w2 = frequencies2(w,1);

    // Energy integration

    cond3shg1(w) = integrate(energies, Eigen::Matrix<std::complex<U>,-1,1>(omega_energies3shg1.col(w)));
    cond3shg2(w) = integrate(energies, Eigen::Matrix<std::complex<U>,-1,1>(omega_energies3shg2.col(w)));
    cond3shg3(w) = integrate(energies, Eigen::Matrix<std::complex<U>,-1,1>(omega_energies3shg3.col(w)));
    cond2shg(w) = integrate(energies, Eigen::Matrix<std::complex<U>,-1,1>(omega_energies2shg.col(w)));
    cond1shg(w) = integrate(energies, Eigen::Matrix<std::complex<U>,-1,1>(omega_energies1shg.col(w)));

    // Divide by the frequencies

    cond3shg1(w) /= (w1 + imaginary*scat)*(w2 + imaginary*scat);
    cond3shg2(w) /= (w1 + imaginary*scat)*(w2 + imaginary*scat);
    cond3shg3(w) /= (w1 + imaginary*scat)*(w2 + imaginary*scat);
    cond2shg(w) /= (w1 + imaginary*scat)*(w2 + imaginary*scat);
    cond1shg(w) /= (w1 + imaginary*scat)*(w2 + imaginary*scat);
  }

  U compat_factor = 2.0;//*systemInfo.energy_scale*systemInfo.energy_scale;
  std::complex<U> factor = imaginary*U(systemInfo.num_orbitals*
      systemInfo.spin_degeneracy/systemInfo.unit_cell_area/systemInfo.energy_scale)*compat_factor;


  cond3shg1 *= factor;
  cond3shg2 *= factor;
  cond3shg3 *= factor;
  cond2shg *= factor;
  cond1shg *= factor;

  cond_shg = cond1shg + cond2shg + cond3shg1 + cond3shg2 + cond3shg3;

  std::ofstream myfile_shg;
  std::ofstream myfile3shg1, myfile3shg2, myfile3shg3, myfile2shg, myfile1shg;

  myfile_shg.open(filename);

  if(print_all){

    myfile3shg1.open(filename + "RA");
    myfile3shg2.open(filename + "RR");
    myfile3shg3.open(filename + "AA");
    myfile2shg.open(filename + "2");
    myfile1shg.open(filename + "1");
  }


  for(int i=0; i < N_omegas; i++){
    freq = std::real(frequencies(i))*systemInfo.energy_scale;
    myfile_shg  << freq << " " << cond_shg.real()(i) << " " << cond_shg.imag()(i) << "\n";
    if(print_all){

      myfile3shg1 << freq << " " << cond3shg1.real()(i) << " " << cond3shg1.imag()(i) << "\n";
      myfile3shg2 << freq << " " << cond3shg2.real()(i) << " " << cond3shg2.imag()(i) << "\n";
      myfile3shg3 << freq << " " << cond3shg3.real()(i) << " " << cond3shg3.imag()(i) << "\n";
      myfile2shg << freq << " " << cond2shg.real()(i) << " " << cond2shg.imag()(i) << "\n";
      myfile1shg << freq << " " << cond1shg.real()(i) << " " << cond1shg.imag()(i) << "\n";
    } 
  }
  myfile_shg.close();
  myfile3shg1.close();
  myfile3shg2.close();
  myfile3shg3.close();
  myfile2shg.close();
  myfile1shg.close();
}

template <typename U, unsigned DIM>
void conductivity_nonlinear<U, DIM>::calculate(){
  debug_message("Entered calc_nonlinear_cond.\n");
  //Calculates the nonlinear conductivity for a set of frequencies in the range [-sigma, sigma].
  //These frequencies are in the KPM scale, that is, the scale where the energy is in the range ]-1,1[.
  //the temperature is already in the KPM scale, but not the broadening or the Fermi Energy

  int photo = 0;
  if(ratio == -1.0) photo = 1;


  energies     = Eigen::Matrix<U, -1, 1>::LinSpaced(N_energies, -lim, lim);
  frequencies  = Eigen::Matrix<U, -1, 1>::LinSpaced(N_omegas, minFreq, maxFreq);
  frequencies2 = Eigen::Matrix<U, -1, 2>::Zero(N_omegas, 2);


  for(int w = 0; w < N_omegas; w++){
    frequencies2(w,0) = frequencies(w);
    frequencies2(w,1) = ratio*frequencies(w);
  }

  if(photo)     calculate_photo();
  if(not photo) calculate_general();
  
  debug_message("Left calc_nonlinear_cond.\n");
}

// Instantiations
template class conductivity_nonlinear<float, 1u>;
template class conductivity_nonlinear<float, 2u>;
template class conductivity_nonlinear<float, 3u>;

template class conductivity_nonlinear<double, 1u>;
template class conductivity_nonlinear<double, 2u>;
template class conductivity_nonlinear<double, 3u>;

template class conductivity_nonlinear<long double, 1u>;
template class conductivity_nonlinear<long double, 2u>;
template class conductivity_nonlinear<long double, 3u>;
