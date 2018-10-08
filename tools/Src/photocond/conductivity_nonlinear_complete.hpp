/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

#include "../headers.hpp"

using std::cout;
using namespace H5;
using std::endl;

template <typename T, unsigned DIM>
class conductivity_nonlinear{
	H5::H5File file;
	public:


    bool isRequired = false; // was this quantity (conductivity_dc) asked for?
    bool isPossible = false; // do we have all we need to calculate the conductivity?


    // Parameters
    int direction;
    int NumDisorder;
    int NumMoments;
    int NumPoints = -1;
    int NumRandoms;
    double temperature = -1;
    int special = 0;
    std::string filename = "nonlinear_cond.dat";


    // 1/kT, where k is the Boltzmann constant in eV/K
    T beta;
    T e_fermi;
    T scat;
    
    // Energy parameters needed to run the simulation
    int N_energies;
    double lim;
    Eigen::Matrix<T, -1, 1> energies;

    // Frequency parameters needed to run the simulation
    int N_omegas;
    double minFreq;
    double maxFreq;
    Eigen::Matrix<T, -1, 1> frequencies;
    Eigen::Matrix<T, -1, 2> frequencies2;

    // information about the Hamiltonian
    system_info<T, DIM> systemInfo;

    // Input from the shell to override the configuration file
    shell_input variables;

    // Objects required to successfully calculate the conductivity
		Eigen::Array<std::complex<T>, -1, -1> Gamma0;
		Eigen::Array<std::complex<T>, -1, -1> Gamma1;
		Eigen::Array<std::complex<T>, -1, -1> Gamma2;
		Eigen::Array<std::complex<T>, -1, -1> Gamma3;

	  std::string dirName;


	  std::complex<T> imaginary;
    conductivity_nonlinear(system_info<T, DIM>&, shell_input &);
	void fetch_parameters();
	void override_parameters();
    void calculate();

    Eigen::Matrix<std::complex<T>, -1, -1> Gamma0contract();

    Eigen::Matrix<std::complex<T>, -1, -1> Gamma3Contract_RRandAA();
    Eigen::Matrix<std::complex<T>, -1, -1> Gamma3Contract_RRandAAblocks();
    Eigen::Matrix<std::complex<T>, -1, -1> Gamma3Contract_RA();
    Eigen::Matrix<std::complex<T>, -1, -1> Gamma1contractAandR();
    Eigen::Matrix<std::complex<T>, -1, -1> Gamma2contractAandR();
	
    Eigen::Matrix<std::complex<T>, -1, -1> Gamma1shgcontractAandR();
    Eigen::Matrix<std::complex<T>, -1, -1> Gamma2shgcontractAandR();
    Eigen::Matrix<std::complex<T>, -1, -1> Gamma3shgContract_RA();
    Eigen::Matrix<std::complex<T>, -1, -1> Gamma3shgContract_RR();
    Eigen::Matrix<std::complex<T>, -1, -1> Gamma3shgContract_AA();

};

template <typename T, unsigned DIM>
conductivity_nonlinear<T, DIM>::conductivity_nonlinear(system_info<T, DIM>& info, shell_input & vari){
  // Constructor of the conductivity_nonlinear class. This function simply checks
  // whether the conductivity needs to be calculated or not

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
  

}
	
template <typename T, unsigned DIM>
void conductivity_nonlinear<T, DIM>::override_parameters(){
    if(variables.CondOpt2_Temp != -8888)     temperature = variables.CondOpt2_Temp/systemInfo.energy_scale;
    if(variables.CondOpt2_NumEnergies != -1) N_energies  = variables.CondOpt2_NumEnergies;
    if(variables.CondOpt2_FreqMax != -8888)  maxFreq     = variables.CondOpt2_FreqMax/systemInfo.energy_scale;
    if(variables.CondOpt2_FreqMin != -8888)  minFreq     = variables.CondOpt2_FreqMin/systemInfo.energy_scale;
    if(variables.CondOpt2_NumFreq != -1)     N_omegas    = variables.CondOpt2_NumFreq;
    if(variables.CondOpt2_Fermi != -8888)    e_fermi     = variables.CondOpt2_Fermi/systemInfo.energy_scale;
    if(variables.CondOpt2_Scat != -8888)     scat        = variables.CondOpt2_Scat/systemInfo.energy_scale;
    if(variables.CondOpt2_Name != "")        filename    = variables.CondOpt2_Name;
    beta = 1.0/8.6173303*pow(10,5)/temperature;
};

template <typename T, unsigned DIM>
void conductivity_nonlinear<T, DIM>::fetch_parameters(){
	debug_message("Entered conductivity_nonlinear::read.\n");
  // This function reads all the relevant
  // information from the hdf5 configuration file and uses it to evaluate the parameters
  // needed to calculate the nonlinear conductivity
	 
    variables.printOpt2();

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
  std::string dirString = num2str3(direction);

  // Fetch the number of Chebyshev Moments, temperature and number of points
    get_hdf5(&NumMoments, &file, (char*)(dirName+"NumMoments").c_str());	
	get_hdf5(&temperature, &file, (char*)(dirName+"Temperature").c_str());	
	get_hdf5(&NumPoints, &file, (char*)(dirName+"NumPoints").c_str());	
	get_hdf5(&special, &file, (char*)(dirName+"Special").c_str());	

  // 1/kT, where k is the Boltzmann constant in eV/K
  beta = 1.0/8.6173303*pow(10,5)/temperature;
  e_fermi = 0.0/systemInfo.energy_scale;
  scat = 0.0166/systemInfo.energy_scale;
	
  // Energy parameters needed to run the simulation
  N_energies = 1024;
  lim = 0.995;

  // Frequency parameters needed to run the simulation
	N_omegas = NumPoints;
  minFreq = 0;
  maxFreq = 2.0;

  // Check whether the matrices we're going to retrieve are complex or not
  int complex = systemInfo.isComplex;


  bool hasGamma0 = false; 
  bool hasGamma1 = false;
  bool hasGamma2 = false; 
  bool hasGamma3 = false;


  debug_message("special:" ); debug_message(special); debug_message("\n");
  debug_message("dirstring: "); debug_message(dirString); debug_message("\n");
  
  // Retrieve the Gamma0 Matrix. This is the 1D Gamma matrix defined as
  // Gamma0 = Tr[v^abc T_n]. This matrix is not needed for the special case
  std::string MatrixName = dirName + "Gamma0" + dirString;
  if(special == 0){
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

      hasGamma0 = true;
    } catch(H5::Exception& e) {
      debug_message("Conductivity nonlinear: There is no Gamma0 matrix.\n");
    }
  }





  // Retrieve the Gamma1 Matrix. This is the 2D Gamma matrix defined as
  // Gamma1 = Tr[v^a Tn v^bc Tm]
  MatrixName = dirName + "Gamma1" + dirString;
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

    hasGamma1 = true;
  } catch(H5::Exception& e) {
    debug_message("Conductivity optical: There is no Gamma1 matrix.\n");
  }



  // Retrieve the Gamma2 Matrix
  MatrixName = dirName + "Gamma2" + dirString;
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

    hasGamma2 = true;
  } catch(H5::Exception& e) {
    debug_message("Conductivity optical: There is no Gamma2 matrix.\n");
  }


  

  // Retrieve the Gamma3 Matrix. This is the biggest matrix and doesn't need to be
  // calculated when we want hBN because it is identically zero.
  if(special == 0){
    MatrixName = dirName + "Gamma3" + dirString;
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

      hasGamma3 = true;
    } catch(H5::Exception& e) {
      debug_message("Conductivity DC: There is no Gamma3 matrix.\n");
    }
  }

  // check if we have all the objects that we need
  if(special == 1){
    if(hasGamma1 and hasGamma2){
      isPossible = true;
    }
  } else {
    if(special == 0){
      if(hasGamma0 and hasGamma1 and hasGamma2 and hasGamma3){
        isPossible = true;
      }
    }
  }


	file.close();
	debug_message("Left conductivity_nonlinear::read.\n");
}

template <typename U, unsigned DIM>
void conductivity_nonlinear<U, DIM>::calculate(){
	debug_message("Entered calc_nonlinear_cond.\n");
	//Calculates the nonlinear conductivity for a set of frequencies in the range [-sigma, sigma].
	//These frequencies are in the KPM scale, that is, the scale where the energy is in the range ]-1,1[.
	//the temperature is already in the KPM scale, but not the broadening or the Fermi Energy

  // ########################################################
  // default values for the fermi energy, broadening and frequencies
  // ########################################################


  // Print out some useful information
  verbose_message("  Energy in rescaled units: [-1,1]\n");
  verbose_message("  Beta (1/kT): "); verbose_message(beta); verbose_message("\n");
  verbose_message("  Fermi energy: "); verbose_message(e_fermi); verbose_message("\n");
  verbose_message("  Using kernel for delta function: Jackson\n");
  verbose_message("  Broadening parameter for Green's function: ");
    verbose_message(scat); verbose_message("\n");
  verbose_message("  Number of energies: "); verbose_message(N_energies); verbose_message("\n");
  verbose_message("  Energy range: ["); verbose_message(-lim); verbose_message(",");
    verbose_message(lim); verbose_message("]\n");
  verbose_message("  Number of frequencies: "); verbose_message(N_omegas); verbose_message("\n");
  verbose_message("  Frequency range: ["); verbose_message(minFreq); verbose_message(",");
    verbose_message(maxFreq); verbose_message("]\n");
  verbose_message("  File name: ");verbose_message(filename); verbose_message("\n");

  energies    = Eigen::Matrix<U, -1, 1>::LinSpaced(N_energies, -lim, lim);
  frequencies = Eigen::Matrix<U, -1, 1>::LinSpaced(N_omegas, minFreq, maxFreq);
  frequencies2 = Eigen::Matrix<U, -1, 2>::Zero(N_omegas, 2);

  for(int w = 0; w < N_omegas; w++){
    frequencies2(w,0) = frequencies(w);
    frequencies2(w,1) = frequencies(w);
  }
  
  Eigen::Matrix<std::complex<U>, -1, -1> omega_energies0, omega_energies1, omega_energies2, omega_energies3, omega_energies4;
  Eigen::Matrix<std::complex<U>, -1, -1> omega_energies3shg1, omega_energies3shg2, omega_energies3shg3;
  Eigen::Matrix<std::complex<U>, -1, -1> omega_energies2shg, omega_energies1shg;
  omega_energies0 = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  omega_energies1 = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  omega_energies2 = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  omega_energies3 = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  omega_energies4 = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);

  omega_energies1shg  = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  omega_energies2shg  = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  omega_energies3shg1 = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  omega_energies3shg2 = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);
  omega_energies3shg3 = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, N_omegas);


  Eigen::Matrix<std::complex<U>,1,-1> cond0, cond1, cond2, cond3, cond4, cond;
  Eigen::Matrix<std::complex<U>,1,-1> cond3shg1, cond3shg2, cond3shg3;
  Eigen::Matrix<std::complex<U>,1,-1> cond1shg, cond2shg;
  cond0 = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);
  cond1 = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);
  cond2 = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);
  cond3 = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);
  cond4 = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);


  cond1shg  = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);
  cond2shg  = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);
  cond3shg1 = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);
  cond3shg2 = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);
  cond3shg3 = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);

  cond  = Eigen::Matrix<std::complex<U>, 1, -1>::Zero(1, N_omegas);

  // Contraction of the Gamma matrices with the delta functions and Green's functions
  omega_energies0 += Gamma0contract();
  omega_energies1 += 0.5*Gamma1contractAandR(); 
  omega_energies2 += Gamma2contractAandR(); 

  omega_energies1shg += 0.5*Gamma1shgcontractAandR();
  omega_energies2shg += Gamma2shgcontractAandR();

  //special = 1;
  if(special != 1){
    omega_energies3shg1 += Gamma3shgContract_RA();
    omega_energies3shg2 += Gamma3shgContract_RR();
    omega_energies3shg3 += Gamma3shgContract_AA();
    omega_energies4 += Gamma3Contract_RA(); 
    omega_energies3 += Gamma3Contract_RRandAAblocks();
  }

  U freq;
  U w1, w2;
  for(int w = 0; w < N_omegas; w++){
    freq = frequencies(w);  
    w1 = frequencies2(w,0);
    w2 = frequencies2(w,1);
    cond0(w) = integrate(energies, Eigen::Matrix<std::complex<U>,-1,1>(omega_energies0.col(w)));
    cond1(w) = integrate(energies, Eigen::Matrix<std::complex<U>,-1,1>(omega_energies1.col(w)));
    cond2(w) = integrate(energies, Eigen::Matrix<std::complex<U>,-1,1>(omega_energies2.col(w)));
    cond3(w) = integrate(energies, Eigen::Matrix<std::complex<U>,-1,1>(omega_energies3.col(w)));
    cond4(w) = integrate(energies, Eigen::Matrix<std::complex<U>,-1,1>(omega_energies4.col(w)));

    cond3shg1(w) = integrate(energies, Eigen::Matrix<std::complex<U>,-1,1>(omega_energies3shg1.col(w)));
    cond3shg2(w) = integrate(energies, Eigen::Matrix<std::complex<U>,-1,1>(omega_energies3shg2.col(w)));
    cond3shg3(w) = integrate(energies, Eigen::Matrix<std::complex<U>,-1,1>(omega_energies3shg3.col(w)));
    cond2shg(w) = integrate(energies, Eigen::Matrix<std::complex<U>,-1,1>(omega_energies2shg.col(w)));
    cond1shg(w) = integrate(energies, Eigen::Matrix<std::complex<U>,-1,1>(omega_energies1shg.col(w)));

    cond0(w) /= -scat*scat - freq*freq; 
    cond1(w) /= -scat*scat - freq*freq; 
    cond2(w) /= -scat*scat - freq*freq; 
    cond3(w) /= -scat*scat - freq*freq; 
    cond4(w) /= -scat*scat - freq*freq; 

    cond3shg1(w) /= (w1 + imaginary*scat)*(w2 + imaginary*scat);
    cond3shg2(w) /= (w1 + imaginary*scat)*(w2 + imaginary*scat);
    cond3shg3(w) /= (w1 + imaginary*scat)*(w2 + imaginary*scat);
    cond2shg(w) /= (w1 + imaginary*scat)*(w2 + imaginary*scat);
    cond1shg(w) /= (w1 + imaginary*scat)*(w2 + imaginary*scat);
  }

  std::complex<U> factor = imaginary*U(systemInfo.num_orbitals*
      systemInfo.spin_degeneracy/systemInfo.unit_cell_area/systemInfo.energy_scale);
  cond0 *= factor;
  cond1 *= factor;
  cond2 *= factor;
  cond3 *= factor;
  cond4 *= factor;

  cond3shg1 *= factor;
  cond3shg2 *= factor;
  cond3shg3 *= factor;
  cond2shg *= factor;
  cond1shg *= factor;

  cond =/* cond0 +*/ cond1 + cond2 + cond3 + cond4;

  std::ofstream myfile0, myfile1, myfile2, myfile3, myfile4, myfile;
  std::ofstream myfile3shg1, myfile3shg2, myfile3shg3, myfile2shg, myfile1shg;
  myfile.open(filename);
  myfile0.open ("nonlinear_cond0.dat");
  myfile1.open ("nonlinear_cond1.dat");
  myfile2.open ("nonlinear_cond2.dat");
  myfile3.open ("nonlinear_cond3.dat");
  myfile4.open ("nonlinear_cond4.dat");
  myfile3shg1.open ("nonlinear_cond3shg1.dat");
  myfile3shg2.open ("nonlinear_cond3shg2.dat");
  myfile3shg3.open ("nonlinear_cond3shg3.dat");
  myfile2shg.open ("nonlinear_cond2shg.dat");
  myfile1shg.open ("nonlinear_cond1shg.dat");


  for(int i=0; i < N_omegas; i++){
    freq = std::real(frequencies(i))*systemInfo.energy_scale;
    myfile  << freq << " " << cond.real()(i) << " " << cond.imag()(i) << "\n";
    myfile0 << freq << " " << cond0.real()(i) << " " << cond0.imag()(i) << "\n";
    myfile1 << freq << " " << cond1.real()(i) << " " << cond1.imag()(i) << "\n";
    myfile2 << freq << " " << cond2.real()(i) << " " << cond2.imag()(i) << "\n";
    myfile3 << freq << " " << cond3.real()(i) << " " << cond3.imag()(i) << "\n";
    myfile4 << freq << " " << cond4.real()(i) << " " << cond4.imag()(i) << "\n";

    myfile3shg1 << freq << " " << cond3shg1.real()(i) << " " << cond3shg1.imag()(i) << "\n";
    myfile3shg2 << freq << " " << cond3shg2.real()(i) << " " << cond3shg2.imag()(i) << "\n";
    myfile3shg3 << freq << " " << cond3shg3.real()(i) << " " << cond3shg3.imag()(i) << "\n";
    myfile2shg << freq << " " << cond2shg.real()(i) << " " << cond2shg.imag()(i) << "\n";
    myfile1shg << freq << " " << cond1shg.real()(i) << " " << cond1shg.imag()(i) << "\n";
  }
  myfile.close();
  myfile0.close();
  myfile1.close();
  myfile2.close();
  myfile3.close();
  myfile4.close();

  myfile3shg1.close();
  myfile3shg2.close();
  myfile3shg3.close();
  myfile2shg.close();
  myfile1shg.close();
  debug_message("Left calc_nonlinear_cond.\n");
};

#include "Gamma0.hpp"
#include "Gamma1.hpp"
#include "Gamma2.hpp"
#include "Gamma3.hpp"
#include "Gamma3shg.hpp"
#include "Gamma2shg.hpp"
#include "Gamma1shg.hpp"
