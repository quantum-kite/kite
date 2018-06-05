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


    // information about the Hamiltonian
    system_info<T, DIM> systemInfo;

    // Objects required to successfully calculate the conductivity
		Eigen::Array<std::complex<T>, -1, -1> Gamma;
		Eigen::Array<std::complex<T>, -1, -1> Lambda;

	  std::string dirName;


    conductivity_optical(system_info<T, DIM>&);
		void read();
    void calculate();
	
};

template <typename T, unsigned DIM>
conductivity_optical<T, DIM>::conductivity_optical(system_info<T, DIM>& info){
  std::string name = info.filename;
	file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

  // retrieve the information about the Hamiltonian
  systemInfo = info;

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
void conductivity_optical<T, DIM>::read(){
	debug_message("Entered conductivity_optical::read.\n");
	//This function reads all the data from the hdf5 file that's needed to 
  //calculate the optical conductivity
	 

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
    debug_message("Conductivity DC: There is no Gamma matrix.\n");
    isPossible = false;
  }






  // Retrieve the Lambda Matrix
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

	



	file.close();
	debug_message("Left conductivity_optical::read.\n");
}

template <typename U, unsigned DIM>
void conductivity_optical<U, DIM>::calculate(){
	debug_message("Entered calc_optical_cond.\n");
	//Calculates the optical conductivity for a set of frequencies in the range [-sigma, sigma].
	//These frequencies are in the KPM scale, that is, the scale where the energy is in the range ]-1,1[.
	//the temperature is already in the KPM scale, but not the broadening or the Fermi Energy

  // ########################################################
  // default values for the fermi energy, broadening and frequencies
  // ########################################################

  // 1/kT, where k is the Boltzmann constant in eV/K
  U beta = 1.0/8.6173303*pow(10,5)/temperature;
  U e_fermi = 0.1;
  U scat = 0.01;
	
	// Calculate the number of frequencies and energies needed to perform the calculation.
	int N_energies = NumPoints;
  double lim = 0.99;
  Eigen::Matrix<U, -1, 1> energies;
  energies  = Eigen::Matrix<U, -1, 1>::LinSpaced(N_energies, -lim, lim);

	int N_omegas = 200;
  double minFreq = 0.01;
  double maxFreq = 1.5;
  Eigen::Matrix<U, -1, 1> frequencies;
  frequencies = Eigen::Matrix<U, -1, 1>::LinSpaced(N_omegas, minFreq, maxFreq);

  // Print out some useful information
  verbose_message("  Beta (1/kT) (in KPM units): "); verbose_message(beta); verbose_message("\n");
  verbose_message("  Fermi energi (in KPM units): "); verbose_message(e_fermi); verbose_message("\n");
  verbose_message("  Using kernel for delta function: Jackson\n");
  verbose_message("  Using broadening parameter for Green's function (in KPM units): ");
    verbose_message(scat); verbose_message("\n");
  verbose_message("  Number of energies: "); verbose_message(NumPoints); verbose_message("\n");
  verbose_message("  Energy range (in KPM units): ["); verbose_message(-lim); verbose_message(",");
    verbose_message(lim); verbose_message("]\n");
  verbose_message("  Number of frequencies: "); verbose_message(N_omegas); verbose_message("\n");
  verbose_message("  Frequency range (in KPM units): ["); verbose_message(minFreq); verbose_message(",");
    verbose_message(maxFreq); verbose_message("]\n");
  verbose_message("  File name: optical_cond.dat\n");



	std::complex<U> imaginary(0.0, 1.0);

  
  // Functions that are going to be used by the contractor
  int NumMoments1 = NumMoments;
  std::function<U(int, U)> deltaF = [beta, e_fermi, NumMoments1](int n, U energy)->U{
    return delta(n, energy)*U(1.0/(1.0 + U(n==0)))*fermi_function(energy, e_fermi, beta)*kernel_jackson<U>(n, NumMoments1);
  };


  Eigen::Matrix<std::complex<U>, -1, 1> temp1; 
  Eigen::Matrix<std::complex<U>, -1, 1> temp2; 
  Eigen::Matrix<std::complex<U>, -1, 1> cond; 
  std::complex<U> temp3 = 0; 

  temp1 = Eigen::Matrix<std::complex<U>, -1, 1>::Zero(N_omegas, 1);
  temp2 = Eigen::Matrix<std::complex<U>, -1, 1>::Zero(N_omegas, 1);
  cond  = Eigen::Matrix<std::complex<U>, -1, 1>::Zero(N_omegas, 1);


  temp1 = contract2<U>(deltaF, 0, greenAscat<U>(scat), NumMoments, Gamma, energies, -frequencies);
  temp2 = contract2<U>(deltaF, 1, greenRscat<U>(scat), NumMoments, Gamma, energies, frequencies);
  temp3 = contract1<U>(deltaF, NumMoments, Lambda, energies);

  std::complex<U> freq;
  for(int i = 0; i < N_omegas; i++){
    freq = std::complex<U>(frequencies(i), scat);  
    cond(i) += (temp1(i) + temp2(i) + temp3)/freq;
  }
  cond *= imaginary*U(systemInfo.num_orbitals*systemInfo.spin_degeneracy/systemInfo.unit_cell_area/units);

	
  
  //Output to a file
  std::ofstream myfile;
  myfile.open ("optical_cond.dat");
  for(int i=0; i < N_omegas; i++){
    freq = std::complex<U>(frequencies(i), scat);  
    myfile << frequencies.real()(i)*systemInfo.energy_scale << " " << cond.real()(i) << " " << cond.imag()(i) << "\n";
  }
  myfile.close();
  debug_message("Left calc_optical_cond.\n");



};
