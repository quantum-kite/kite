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

template <typename T, unsigned DIM>
class conductivity_optical{
	H5::H5File file;
	public:


    bool isRequired = false; // was this quantity (conductivity_dc) asked for?
    bool isPossible = false; // do we have all we need to calculate the conductivity?


		// Functions to calculate. They will require the objects present in
    // the configuration file
    int direction;
    int NumDisorder;
    int NumMoments;
    int NumPoints = -1;
    int NumRandoms;
    double temperature = -1;
    double units = unit_scale;
    std::string filename = "optical_cond.dat";

    T beta;
    T e_fermi; 
    T scat; 
	int N_energies; 
	int N_omegas; 
    double minFreq; 
    double maxFreq; 
    T lim;


    // information about the Hamiltonian
    system_info<T, DIM> systemInfo;

    // Objects required to successfully calculate the conductivity
    shell_input variables;

    // Objects required to successfully calculate the conductivity
		Eigen::Array<std::complex<T>, -1, -1> Gamma;
		Eigen::Array<std::complex<T>, -1, -1> Lambda;

	  std::string dirName;


    conductivity_optical(system_info<T, DIM>&, shell_input &);
	void fetch_parameters();
    void override_parameters();
    void calculate();
    void calculate_efficient();
	
};

template <typename T, unsigned DIM>
conductivity_optical<T, DIM>::conductivity_optical(system_info<T, DIM>& info, shell_input & vari){
  std::string name = info.filename;
	file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

  // retrieve the information about the Hamiltonian
  systemInfo = info;

  // location of the information about the conductivity
  variables = vari;

  // location of the information about the conductivity
  dirName = "/Calculation/conductivity_optical/";
  
  // check whether the conductivity_optical was asked for
  try{
    H5::Exception::dontPrint();
    get_hdf5(&direction, &file, (char*)(dirName+"Direction").c_str());									
    isRequired = true;
  } catch(H5::Exception& e){}
  

}
	

template <typename T, unsigned DIM>
void conductivity_optical<T, DIM>::fetch_parameters(){
	debug_message("Entered conductivity_optical::read.\n");
	//This function reads all the data from the hdf5 file that's needed to 
  //calculate the optical conductivity
	 

  // Check if the data for the optical conductivity exists
  // Check if the data for the optical conductivity exists
  if(!isRequired){
    std::cout << "Data for optical conductivity does not exist. Exiting.\n";
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
    debug_message("Conductivity optical: There is no Gamma matrix.\n");
    isPossible = false;
  }




  MatrixName = dirName + "Lambda" + dirString;
  try{
		debug_message("Filling the Lambda matrix.\n");
		Lambda = Eigen::Array<std::complex<T>,-1,-1>::Zero(1, NumMoments);
		
		if(complex)
			get_hdf5(Lambda.data(), &file, (char*)MatrixName.c_str());
		
		if(!complex){
			Eigen::Array<T,-1,-1> LambdaReal;
			LambdaReal = Eigen::Array<T,-1,-1>::Zero(NumMoments, NumMoments);
			get_hdf5(LambdaReal.data(), &file, (char*)MatrixName.c_str());
			
			Lambda = LambdaReal.template cast<std::complex<T>>();
		}				

    isPossible = true;
  } catch(H5::Exception& e) {
    debug_message("Conductivity optical: There is no Lambda matrix.\n");
    isPossible = false;
  }

	
  temperature = 0.001/systemInfo.energy_scale;
  beta = 1.0/8.6173303*pow(10,5)/temperature;
  e_fermi = 0.2/systemInfo.energy_scale;
  scat = 0.0166/systemInfo.energy_scale;


  N_energies = 512;

  N_omegas = NumPoints;
  minFreq = 0.001;
  maxFreq = 1.5;

  lim = 0.99;

	file.close();
	debug_message("Left conductivity_optical::read.\n");
}


template <typename U, unsigned DIM>
void conductivity_optical<U, DIM>::override_parameters(){
    if(variables.CondOpt_Temp != -8888)     temperature = variables.CondOpt_Temp/systemInfo.energy_scale;
    if(variables.CondOpt_NumEnergies != -1) N_energies  = variables.CondOpt_NumEnergies;
    if(variables.CondOpt_FreqMax != -8888)  maxFreq     = variables.CondOpt_FreqMax/systemInfo.energy_scale;
    if(variables.CondOpt_FreqMin != -8888)  minFreq     = variables.CondOpt_FreqMin/systemInfo.energy_scale;
    if(variables.CondOpt_NumFreq != -1)     N_omegas    = variables.CondOpt_NumFreq;
    if(variables.CondOpt_Fermi != -8888)    e_fermi     = variables.CondOpt_Fermi/systemInfo.energy_scale;
    if(variables.CondOpt_Scat != -8888)     scat        = variables.CondOpt_Scat/systemInfo.energy_scale;
    if(variables.CondOpt_Name != "")        filename    = variables.CondOpt_Name;
    beta = 1.0/8.6173303*pow(10,5)/temperature;
	//Calculates the optical conductivity for a set of frequencies in the range [-sigma, sigma].
};
	//These frequencies are in the KPM scale, that is, the scale where the energy is in the range ]-1,1[.
template <typename U, unsigned DIM>
void conductivity_optical<U, DIM>::calculate_efficient(){
	debug_message("Entered calc_optical_cond.\n");
	//the temperature is already in the KPM scale, but not the broadening or the Fermi Energy

  // ########################################################
  if(!isPossible){
    std::cout << "Cannot calculate the conductivity because there is no matching Gamma matrix"
      "Exiting\n";
    exit(0);
  }


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
  verbose_message("  File name: optical_cond.dat\n");


  Eigen::Matrix<U, -1, 1> energies;
  energies  = Eigen::Matrix<U, -1, 1>::LinSpaced(N_energies, -lim, lim);
  Eigen::Matrix<U, -1, 1> frequencies;
  frequencies = Eigen::Matrix<U, -1, 1>::LinSpaced(N_omegas, minFreq, maxFreq);

	std::complex<U> imaginary(0.0, 1.0);


  
  // Functions that are going to be used by the contractor
  int NumMoments1 = NumMoments;
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
  DeltaMatrix = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(N_energies, NumMoments);
  for(int n = 0; n < NumMoments; n++)
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
  if(NumMoments%N_threads != 0){
    std::cout << "The number of Chebyshev moments in the optical conductivity must"
      "be a multiple of the number of threads\n" << std::flush;
    exit(1);
  }
}
#pragma omp barrier


#pragma omp for schedule(static, 1) nowait
  for(int i = 0; i < N_threads; i++){
    local_NumMoments = NumMoments/N_threads;
    thread_num = omp_get_thread_num();

    // The Gamma matrix has been divided among the threads
    // Each thread has one section of that matrix, called local_Gamma
    Eigen::Matrix<std::complex<U>, -1, -1> local_Gamma;
    local_Gamma = Gamma.matrix().block(0, local_NumMoments*thread_num, 
        NumMoments, local_NumMoments);

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
    for(int w = 0; w < N_omegas; w++){
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

  
  temp3 = contract1<U>(deltaF, NumMoments, Lambda, energies);


  std::complex<U> freq;
  for(int i = 0; i < N_omegas; i++){
    freq = std::complex<U>(frequencies(i), scat);  
    cond(i) += temp3;
    cond(i) /= freq;
  }
  cond *= imaginary*U(systemInfo.num_orbitals*systemInfo.spin_degeneracy/systemInfo.unit_cell_area/units);

	
  
  //Output to a file
  std::ofstream myfile;
  myfile.open(filename);
  for(int i=0; i < N_omegas; i++){
    freq = std::complex<U>(frequencies(i), scat);  
    myfile << frequencies.real()(i)*systemInfo.energy_scale << " " << cond.real()(i) << " " << cond.imag()(i) << "\n";
  }
  myfile.close();
  debug_message("Left calc_optical_cond.\n");



};
