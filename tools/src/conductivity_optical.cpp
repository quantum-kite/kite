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
#include "conductivity_optical.hpp"
#include "functions.hpp"

#include "macros.hpp"

template <typename T, unsigned DIM>
conductivity_optical<T, DIM>::conductivity_optical(system_info<T, DIM>& info, shell_input & vari){
    name = info.filename;

		// Functions to calculate. They will require the objects present in
    // the configuration file
    units = unit_scale;
    systemInfo = info;   // retrieve the information about the Hamiltonian
    variables = vari;    // location of the information about the conductivity
    
    isRequired = is_required() && variables.CondOpt_is_required;  // was this quantity (conductivity_optical) asked for?
    isPossible = false;          // do we have all we need to calculate the conductivity?

    if(isRequired){
        set_default_parameters();   // sets a default set of paramters for the calculation
        isPossible = fetch_parameters();         // finds all the paramters in the .h5 file
        override_parameters();      // overrides parameters with the ones from the shell input

        if(isPossible){
            printOpt();                  // Print all the parameters used

            if(default_Convergence_D and default_Convergence_G){
            calculate();
            } else {
              calculateBlocks();
            }

        } else {
          std::cout << "ERROR. The optical conductivity was requested but the data " "needed for its computation was not found in the input .h5 file. "
              "Make sure KITEx has processed the file first. Exiting.";
          exit(1);

        }
    }
}

template <typename T, unsigned DIM>
void conductivity_optical<T, DIM>::printOpt(){
  double scale = systemInfo.energy_scale;
  std::cout << "The optical conductivity will be calculated with these parameters: (eV, Kelvin)\n"
    "   Temperature: "            << temperature*scale  << ((default_temperature)?  " (default)":"") << "\n"
    "   Broadening: "             << scat*scale         << ((default_scat)?         " (default)":"") << "\n"
    "   Fermi energy: "           << e_fermi*scale      << ((default_efermi)?       " (default)":"") << "\n"
    "   Min frequency: "          << minFreq*scale      << ((default_minfreqs)?     " (default)":"") << "\n"
    "   Max frequency: "          << maxFreq*scale      << ((default_maxfreqs)?     " (default)":"") << "\n"
    "   Number of frequencies: "  << N_omegas           << ((default_Nfreqs)?       " (default)":"") << "\n"
    "   Number of polynomials: "  << NumMoments         << ((default_NumMoments)?   " (from file)":"") << "\n"
    "   Number of delta blocks: " << Convergence_D      << ((default_Convergence_D)?   " (default)":"") << "\n"
    "   Number of green blocks: " << Convergence_G      << ((default_Convergence_G)?   " (default)":"") << "\n"
    "   Num integration points: " << N_energies         << ((default_NEnergies)?    " (default)":"") << "\n"
    "   Kernel for Dirac deltas: "<< "jackson"          << " (default)\n"
    "   Filename: "               << filename           << ((default_filename)?     " (default)":"") << "\n";
  if(!Moments_divisible)
    std::cout << "Warning! The number of moments specified does not divide the number of delta blocks"
      " and the number of green blocks. The Gamma matrix has been padded with zeros so that this requirement"
      " is satisfied. This does not change the output but may cause a slight decrease in efficiency.\n";
}
	
template <typename T, unsigned DIM>
void conductivity_optical<T, DIM>::set_default_parameters(){
    // Sets default values for the parameters used in the 
    // calculation of the density of stats. These are the parameters
    // that will be overwritten by the config file and the
    // shell input parameters. 

    double scale = systemInfo.energy_scale;
    double shift = systemInfo.energy_shift;

    // Temperature
    temperature = 0.001/scale;
    beta = 1.0/8.6173303*pow(10,5)/temperature;
    default_temperature = true;

    Convergence_G = systemInfo.NumThreads;
    Convergence_D = 1;
    default_Convergence_D = true;
    default_Convergence_G = true;


    e_fermi = (0.2 - shift)/scale;
    default_efermi = true;

    scat = 0.0166/scale;
    default_scat = true;

    N_energies = 512;
    default_NEnergies = true;

    N_omegas = 512;
    minFreq = 0.001;
    maxFreq = 1.5;
    default_minfreqs = true;
    default_maxfreqs = true;
    default_Nfreqs = true;

    filename  = "optcond.dat";      // Filename to save final result
    default_filename = true;

    lim = 0.99;
}

template <typename T, unsigned DIM>
bool conductivity_optical<T, DIM>::is_required(){
    // Checks whether the optical conductivity has been requested
    // by analysing the .h5 config file. If it has been requested, 
    // some fields have to exist, such as "Direction"

    // Make sure the config filename has been initialized
    if(name == ""){
        std::cout << "ERROR: Filename uninitialized. Exiting.\n";
        exit(1);
    }

    H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

    std::string dirName = "/Calculation/conductivity_optical/";
    bool result = false;
    try{
        get_hdf5(&NumMoments, &file, (char*)(dirName+"NumMoments").c_str());
        result = true;
    } catch(H5::Exception& e){}


    file.close();
    return result;
}


template <typename T, unsigned DIM>
bool conductivity_optical<T, DIM>::fetch_parameters(){
	debug_message("Entered conductivity_optical::fetch_parameters.\n");
	//This function reads all the data from the hdf5 file that's needed to 
  //calculate the optical conductivity
	 

  H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

  // Check if the data for the optical conductivity exists
  if(!isRequired){
    std::cout << "Data for optical conductivity does not exist. Exiting.\n";
    exit(1);
  }
  
	debug_message("Fetching direction.\n");
  // Fetch the direction of the conductivity and convert it to a string
  std::string dirName = "/Calculation/conductivity_optical/";
  get_hdf5(&direction, &file, (char*)(dirName+"Direction").c_str());									
  std::string dirString = num2str2(direction);

	debug_message("Fetching NumMoments, Temperature and NumPoints.\n");
  // Fetch the number of Chebyshev Moments, temperature and number of points
	get_hdf5(&MaxMoments, &file, (char*)(dirName+"NumMoments").c_str());	
	get_hdf5(&temperature, &file, (char*)(dirName+"Temperature").c_str());	
	get_hdf5(&NumPoints, &file, (char*)(dirName+"NumPoints").c_str());	

  NumMoments = MaxMoments;
  default_NumMoments = true;

  N_omegas = NumPoints;
  default_Nfreqs = false;
  default_temperature = false;
  


  // Check whether the matrices we're going to retrieve are complex or not
  int complex = systemInfo.isComplex;


  

  // Retrieve the Gamma Matrix
  std::string MatrixName = dirName + "Gamma" + dirString;
  bool possibleGamma = false;
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

    possibleGamma = true;
  } catch(H5::Exception& e) {
    debug_message("Conductivity optical: There is no Gamma matrix.\n");
  }




  MatrixName = dirName + "Lambda" + dirString;
  bool possibleLambda = false;
  try{
		debug_message("Filling the Lambda matrix.\n");
		Lambda = Eigen::Array<std::complex<T>,-1,-1>::Zero(NumMoments, 1);
		
		if(complex)
			get_hdf5(Lambda.data(), &file, (char*)MatrixName.c_str());
		
		if(!complex){
			Eigen::Array<T,-1,-1> LambdaReal;
			LambdaReal = Eigen::Array<T,-1,-1>::Zero(NumMoments, 1);
			get_hdf5(LambdaReal.data(), &file, (char*)MatrixName.c_str());
			
			Lambda = LambdaReal.template cast<std::complex<T>>();
		}				

    possibleLambda = true;
  } catch(H5::Exception& e) {
    debug_message("Conductivity optical: There is no Lambda matrix.\n");
  }

	

	file.close();
	debug_message("Left conductivity_optical::fetch_parameters.\n");
  return possibleLambda and possibleGamma;
}


template <typename U, unsigned DIM>
void conductivity_optical<U, DIM>::override_parameters(){
  debug_message("Entered override_parameters");
  double scale = systemInfo.energy_scale;
  double shift = systemInfo.energy_shift;

    if(variables.CondOpt_Temp != -8888){
      debug_message("Overriding temperature");
      temperature = variables.CondOpt_Temp/scale;
      beta = 1.0/8.6173303*pow(10,5)/temperature;
      default_temperature = false;
    }

    if(variables.CondOpt_NumEnergies != -1){
      debug_message("Overriding N_energies");
      N_energies  = variables.CondOpt_NumEnergies;
      default_NEnergies = false;
    }

    if(variables.CondOpt_NumMoments != -1){
      debug_message("Overriding NumMoments");
      if(variables.CondOpt_NumMoments > MaxMoments){
        std::cout << "CondOpt: The numbe of Chebyshev moments specified"
          " cannot be larger than the number of moments calculated by KITEx."
          " Please specify a smaller number. Exiting.\n";
        exit(1);
      }
      NumMoments  = variables.CondOpt_NumMoments;
      default_NumMoments = false;
    }

    if(variables.CondOpt_Convergence_D != -1){
      debug_message("Overriding Convergence_D");
      Convergence_D = variables.CondOpt_Convergence_D; 
      default_Convergence_D = false;
    }

    if(variables.CondOpt_Convergence_G != -1){
      debug_message("Overriding Convergence_G");
      Convergence_G = variables.CondOpt_Convergence_G; 
      default_Convergence_G = false;
    }

    if(variables.CondOpt_FreqMax != -8888){
      debug_message("Overriding FreqMax");
      maxFreq     = variables.CondOpt_FreqMax/scale;
      default_maxfreqs = false;
    }

    if(variables.CondOpt_FreqMin != -8888){
      debug_message("Overriding FreqMin");
      minFreq     = variables.CondOpt_FreqMin/scale;
      default_minfreqs = false;
    }

    if(variables.CondOpt_NumFreq != -1){
      debug_message("Overriding NumFreq");
      N_omegas    = variables.CondOpt_NumFreq;
      default_Nfreqs = false;
    }

    if(variables.CondOpt_Fermi != -8888){
      debug_message("Overriding e_fermi");
      e_fermi     = (variables.CondOpt_Fermi - shift)/scale;
      default_efermi = false;
    }

    if(variables.CondOpt_Scat != -8888){
      debug_message("Overriding scat");
      scat        = variables.CondOpt_Scat/scale;
      default_scat = false;
    }

    if(variables.CondOpt_Name != ""){
      debug_message("Overriding filename");
      filename    = variables.CondOpt_Name;
      default_filename = false;
    }

  // It may happen that the user specified a number of polynomials that is
  // not a multiple of the number of threads. We'll pad Gamma with zeros until it is.
  Moments_D = NumMoments;
  Moments_G = NumMoments;
  Moments_divisible = true;

  if(NumMoments % Convergence_D != 0){
    Moments_divisible = false;
    Moments_D = NumMoments + Convergence_D - NumMoments%Convergence_D;
  }
  if(NumMoments % Convergence_G != 0){
    Moments_divisible = false;
    Moments_G = NumMoments + Convergence_G - NumMoments%Convergence_G;
  }
  //Eigen::Matrix<std::complex<U>, -1, -1> Gamma_Padded;
  debug_message("Padding Gamma");
  Gamma_Padded = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(Moments_D, Moments_G);
  Gamma_Padded.block(0,0,NumMoments, NumMoments) = Gamma.block(0,0,NumMoments, NumMoments);

  debug_message("Padding Lambda");
  Lambda_Padded = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(Moments_D, 1);
  Lambda_Padded.block(0,0,NumMoments,1) = Lambda.block(0,0,NumMoments,1);
  debug_message("Left override_parameters");
}

template <typename U, unsigned DIM>
void conductivity_optical<U, DIM>::calculateBlocks(){
  /* Calculate the optical conductivity by blocks. Instead of using the whole
   * Gamma matrix to perform the calculations, split it into several blocks and 
   * compute the contribution to the conductivity of each individual block. This
   * has the potential to slow down the calculation, but allows for assessment of
   * convergence.*/
	debug_message("Entered conductivity_optical::calculateBlocks.\n");
	std::complex<U> imaginary(0.0, 1.0);
  omp_set_num_threads(systemInfo.NumThreads);

  unsigned Nmax, Mmax;
  Nmax = Convergence_D; // lines, Dirac delta
  Mmax = Convergence_G; // columns, Green's function
  unsigned nmax = Moments_D/Nmax;
  unsigned mmax = Moments_G/Mmax;



  // Array of energies to use in the integration, and frequencies at which to
  // evaluate the optical conductivity
  Eigen::Matrix<U, -1, 1> energies;
  Eigen::Matrix<U, -1, 1> frequencies;
  energies    = Eigen::Matrix<U, -1, 1>::LinSpaced(N_energies, -lim, lim);
  frequencies = Eigen::Matrix<U, -1, 1>::LinSpaced(N_omegas, minFreq, maxFreq);


  
  // Function that is going to be used by the contractor
  int NumMoments1 = Moments_D;
  U beta1 = beta;
  U e_fermi1 = e_fermi;
  std::function<U(int, U)> deltaF = [beta1, e_fermi1, NumMoments1](int n, U energy)->U{
    return delta(n, energy)*U(1.0/(1.0 + U(n==0)))*fermi_function(energy, e_fermi1, beta1)*kernel_jackson<U>(n, NumMoments1);
  };


  // The Dirac delta d(e - H) can be expanded in a series of Chebyshev polynomials. The
  // coefficient of this expansion depends on the Chebyshev index and on the energy, and
  // so it can be cast as a matrix. This is what we call the DeltaMatrix
  Eigen::Matrix<std::complex<U>, -1, -1> DeltaMatrix;
  DeltaMatrix = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, Moments_D);
  for(int n = 0; n < Moments_D; n++)
    for(int e = 0; e < N_energies; e++)
      DeltaMatrix(e,n) = deltaF(n, energies(e)); 



  // declarations and memory allocations
  Eigen::Matrix<std::complex<U>, -1, -1>** matrices;
  Eigen::Matrix<std::complex<U>, -1, -1>** deltas;
  Eigen::Matrix<std::complex<U>, -1, -1>** gammas_EN;
  Eigen::Matrix<std::complex<U>, -1, 1>**  gammas_E;
  Eigen::Matrix<std::complex<U>, -1, -1>** cond;
  matrices  = new Eigen::Matrix<std::complex<U>, -1, -1>*[Nmax*Mmax];
  deltas    = new Eigen::Matrix<std::complex<U>, -1, -1>*[Nmax];
  gammas_EN = new Eigen::Matrix<std::complex<U>, -1, -1>*[Nmax*Mmax];
  gammas_E  = new Eigen::Matrix<std::complex<U>, -1,  1>*[Nmax*Mmax];
  cond      = new Eigen::Matrix<std::complex<U>, -1, -1>*[Nmax*Mmax];

  // Copy the Gamma matrix into several different blocks 
  // that will be used later for parallelization. Each thread will be
  // assigned a block. The Gamma matrix is padded with zeros when the number
  // of polynomials does not divide Mmax and Nmax.
  for(unsigned j = 0; j < Nmax*Mmax; j++){
    matrices[j]  = new Eigen::Matrix<std::complex<U>, -1, -1>;

    // block index j = n + Nmax*m
    int n = j % Nmax;
    int m = j / Nmax;
    //*matrices[j] = Gamma.block(n*nmax, m*mmax, nmax, mmax);
    *matrices[j] = Gamma_Padded.block(n*nmax, m*mmax, nmax, mmax);
  }


  // Copy the DeltaMatrix into blocks
  for(unsigned j = 0; j < Nmax; j++){
    deltas[j]  = new Eigen::Matrix<std::complex<U>, -1, -1>;
    *deltas[j] = DeltaMatrix.block(0, j*nmax, N_energies, nmax);
  }



#pragma omp parallel 
{
#pragma omp for schedule(static, 1) nowait
  for(unsigned j = 0; j < Nmax*Mmax; j++){
    // block index j = n + Nmax*m
    int n = j % Nmax;
    int m = j / Nmax;

    // Matrix product of Gamma and Delta
    gammas_EN[j] = new Eigen::Matrix<std::complex<U>, -1, -1>;
    *gammas_EN[j] = (*deltas[n])*(*matrices[j]);


    gammas_E[j]  = new Eigen::Matrix<std::complex<U>, -1, 1>;
    cond[j]      = new Eigen::Matrix<std::complex<U>, -1, -1>;
    *gammas_E[j] = Eigen::Matrix<std::complex<U>, -1, 1>::Zero(N_energies, 1);
    *cond[j]     = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_omegas, 1);


    U freq_p, freq_n;
    std::complex<U> GammaEp, GammaEn;
    for(unsigned w = 0; w < N_omegas; w++){
      freq_p =  frequencies(w);
      freq_n = -frequencies(w);

      for(int e = 0; e < N_energies; e++){
        GammaEp = 0;
        GammaEn = 0;
        for(unsigned int p = 0; p < mmax; p++){
          GammaEp += (*gammas_EN[j])(e, p)*greenAscat<U>(scat)(m*mmax + p, energies(e) - freq_p);      // contracting with the positive frequencies
          GammaEn += (*gammas_EN[j])(e, p)*greenAscat<U>(scat)(m*mmax + p, energies(e) - freq_n);      // contracting with the negative frequencies
        }
        (*gammas_E[j])(e) = GammaEp + std::conj(GammaEn);
      }
      (*cond[j])(w) = integrate(energies, *gammas_E[j]);
    }


    // free all the unneeded allocated memory. The only
    // object we need to keep now is the conductivity
    delete gammas_E[j];  
    delete gammas_EN[j];  
    delete matrices[j];
  }
  
#pragma omp barrier
}
  delete gammas_E;
  delete matrices;
  delete gammas_EN;




  // divide by the frequency
  std::complex<U> freq;
  for(unsigned j = 0; j < Mmax*Nmax; j++){
    for(unsigned i = 0; i < N_omegas; i++){
      freq = std::complex<U>(frequencies(i), scat);  
      (*cond[j])(i) /= freq;
    }
  }
  
  Eigen::Matrix<std::complex<U>, -1, -1> total_cond;
  total_cond = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_omegas, 1);

  Eigen::Matrix<std::complex<U>, -1, -1>** partial_cond;
  partial_cond = new Eigen::Matrix<std::complex<U>, -1, -1>*[Nmax*Mmax];
  std::complex<U> factor = imaginary*U(systemInfo.num_orbitals*systemInfo.spin_degeneracy/systemInfo.unit_cell_area/units);
  
  for(unsigned n = 0; n < Nmax; n++){
    for(unsigned m = 0; m < Mmax; m++){
      partial_cond[n + m*Nmax] = new Eigen::Matrix<std::complex<U>, -1, -1>;
      (*partial_cond[n + m*Nmax]) = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_omegas, 1);
      for(unsigned n1 = 0; n1 <= n; n1++){
        for(unsigned m1 = 0; m1 <= m; m1++){
          (*partial_cond[n + m*Nmax]) += (*cond[n1 + m1*Nmax])*factor;
        }
      }
    }
  }
  

  // Second part of the conductivity
  std::complex<U>* temp3;
  Eigen::Matrix<std::complex<U>, -1, 1>** Lambda_E;
  temp3    = new std::complex<U>[Nmax];
  Lambda_E = new Eigen::Matrix<std::complex<U>, -1, 1>*[Nmax];

  for(unsigned j = 0; j < Nmax; j++){
    Lambda_E[j] = new Eigen::Matrix<std::complex<U>, -1, 1>;
    (*Lambda_E[j]) = (*deltas[j])*Lambda_Padded.matrix().block(nmax*j,0,nmax,1);
    temp3[j] = integrate(energies, *Lambda_E[j]);
  }  
  
  for(unsigned j = 0; j < Nmax; j++){
    delete deltas[j];  
    delete Lambda_E[j];
  }

  delete deltas;
  delete Lambda_E;

  std::complex<U>* partial_temp3;
  partial_temp3 = new std::complex<U>[Nmax];
  partial_temp3[0] = temp3[0];
  for(unsigned j = 1; j < Nmax; j++){
    partial_temp3[j] = temp3[j] + partial_temp3[j-1];
  }
  delete temp3;



  for(unsigned j = 0; j < Mmax*Nmax; j++){
    for(unsigned i = 0; i < N_omegas; i++){
      freq = std::complex<U>(frequencies(i), scat);  
      (*partial_cond[j])(i) += partial_temp3[j%Nmax]/freq*factor;
    }
  }
  total_cond = *partial_cond[Nmax*Mmax - 1];
  
  delete partial_temp3;

  std::ofstream myfile;
  std::complex<U> cn;
  for(unsigned n = 0; n < Nmax; n++){
    for(unsigned m = 0; m < Mmax; m++){
      myfile.open("g" + std::to_string(mmax*(m+1)) + "_d" + std::to_string(nmax*(n+1)) + "_" + filename);
      for(unsigned int i=0; i < N_omegas; i++){
        freq = std::complex<U>(frequencies(i), scat);  
        cn = (*partial_cond[n + m*Nmax])(i);
        myfile << frequencies.real()(i)*systemInfo.energy_scale << " " << cn.real() << " " << cn.imag() << "\n";
      }
      myfile.close();
    }
  }
	
  
  debug_message("Left calc_optical_cond.\n");
  
  for(unsigned j = 0; j < Mmax*Nmax; j++){
    delete cond[j];  
    delete partial_cond[j];
  }
  delete cond;
  delete partial_cond;


}



template <typename U, unsigned DIM>
void conductivity_optical<U, DIM>::calculate(){
	debug_message("Entered calc_optical_cond.\n");
	//the temperature is already in the KPM scale, but not the broadening or the Fermi Energy


  if(!isPossible){
    std::cout << "Cannot calculate the conductivity because there is no matching Gamma matrix"
      "Exiting\n";
    exit(0);
  }



  Eigen::Matrix<U, -1, 1> energies;
  Eigen::Matrix<U, -1, 1> frequencies;
  energies  = Eigen::Matrix<U, -1, 1>::LinSpaced(N_energies, -lim, lim);
  frequencies = Eigen::Matrix<U, -1, 1>::LinSpaced(N_omegas, minFreq, maxFreq);

	std::complex<U> imaginary(0.0, 1.0);


  
  // Functions that are going to be used by the contractor
  int NumMoments1 = Moments_D;
  U beta1 = beta;
  U e_fermi1 = e_fermi;
  std::function<U(int, U)> deltaF = [beta1, e_fermi1, NumMoments1](int n, U energy)->U{
    return delta(n, energy)*U(1.0/(1.0 + U(n==0)))*fermi_function(energy, e_fermi1, beta1)*kernel_jackson<U>(n, NumMoments1);
  };


  Eigen::Matrix<std::complex<U>, -1, 1> cond; 
  cond = Eigen::Matrix<std::complex<U>, -1, 1>::Zero(N_omegas, 1);
  std::complex<U> temp3 = 0; 


  // Delta matrix of chebyshev moments and energies
  Eigen::Matrix<std::complex<U>,-1, -1> DeltaMatrix;
  DeltaMatrix = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, Moments_D);
  for(int n = 0; n < Moments_D; n++)
    for(int e = 0; e < N_energies; e++)
      DeltaMatrix(e,n) = deltaF(n, energies(e)); 


  int N_threads;
  int thread_num;
  int local_NumMoments;
  omp_set_num_threads(systemInfo.NumThreads);
#pragma omp parallel shared(N_threads) private(thread_num)
{
#pragma omp master
{
  N_threads = omp_get_num_threads();
  // check if each thread will get the same number of moments
  //if(NumMoments%N_threads != 0){
    //std::cout << "The number of Chebyshev moments in the optical conductivity must"
      //"be a multiple of the number of threads\n" << std::flush;
    //exit(1);
  //}
}
#pragma omp barrier


#pragma omp for schedule(static, 1) nowait
  for(int i = 0; i < N_threads; i++){
    local_NumMoments = Moments_G/N_threads;
    thread_num = omp_get_thread_num();

    // The Gamma matrix has been divided among the threads
    // Each thread has one section of that matrix, called local_Gamma
    Eigen::Matrix<std::complex<U>, -1, -1> local_Gamma;
    local_Gamma = Gamma.matrix().block(0, local_NumMoments*thread_num, 
        Moments_D, local_NumMoments);

    // Result of contracting the indices with the delta function
    Eigen::Matrix<std::complex<U>, -1, -1> GammaEM;
    GammaEM = DeltaMatrix*local_Gamma;

    Eigen::Matrix<std::complex<U>,-1, 1> GammaEp;
    Eigen::Matrix<std::complex<U>,-1, 1> GammaEn;
    Eigen::Matrix<std::complex<U>,-1, 1> GammaE;
    Eigen::Matrix<std::complex<U>, -1, 1> local_cond; 
    GammaEp = Eigen::Matrix<std::complex<U>,-1, 1>::Zero(N_energies, 1);
    GammaEn = Eigen::Matrix<std::complex<U>,-1, 1>::Zero(N_energies, 1);
    GammaE  = Eigen::Matrix<std::complex<U>,-1, 1>::Zero(N_energies, 1);
    local_cond = Eigen::Matrix<std::complex<U>, -1, 1>::Zero(N_omegas, 1);

    // Loop over the frequencies
    U freq_p, freq_n;
    for(unsigned int w = 0; w < N_omegas; w++){
      freq_p =  frequencies(w);
      freq_n = -frequencies(w);

      for(int e = 0; e < N_energies; e++){
        GammaEp(e) = 0;
        GammaEn(e) = 0;
        for(int m = 0; m < local_NumMoments; m++){
          GammaEp(e) += GammaEM(e, m)*greenAscat<U>(scat)(local_NumMoments*thread_num + m, energies(e) - freq_p);      // contracting with the positive frequencies
          GammaEn(e) += GammaEM(e, m)*greenAscat<U>(scat)(local_NumMoments*thread_num + m, energies(e) - freq_n);      // contracting with the negative frequencies
        }
      }
    
      GammaE = GammaEp + GammaEn.conjugate();                                            // This is equivalent to the two frequency-dependent terms in the formula
      local_cond(w) = integrate(energies, GammaE);
    }

    
#pragma omp critical
    {
      cond += local_cond;
    }
  }
#pragma omp barrier
}

  
  temp3 = contract1<U>(deltaF, Moments_D, Lambda, energies);

  
  //std::cout << "temp3 regular:" << temp3 << "\n";

  std::complex<U> freq;
  for(unsigned int i = 0; i < N_omegas; i++){
    freq = std::complex<U>(frequencies(i), scat);  
    cond(i) += temp3;
    cond(i) /= freq;
  }
  cond *= imaginary*U(systemInfo.num_orbitals*systemInfo.spin_degeneracy/systemInfo.unit_cell_area/units);

	
  
  //Output to a file
  std::ofstream myfile;
  myfile.open(filename);
  for(unsigned int i=0; i < N_omegas; i++){
    freq = std::complex<U>(frequencies(i), scat);  
    myfile << frequencies.real()(i)*systemInfo.energy_scale << " " << cond.real()(i) << " " << cond.imag()(i) << "\n";
  }
  myfile.close();
  debug_message("Left calc_optical_cond.\n");



}

template class conductivity_optical<float, 1u>;
template class conductivity_optical<float, 2u>;
template class conductivity_optical<float, 3u>;

template class conductivity_optical<double, 1u>;
template class conductivity_optical<double, 2u>;
template class conductivity_optical<double, 3u>;

template class conductivity_optical<long double, 1u>;
template class conductivity_optical<long double, 2u>;
template class conductivity_optical<long double, 3u>;
