#include "headers.hpp"

using std::cout;
using namespace H5;
using std::endl;

template <typename T, unsigned DIM>
class conductivity_nonlinear{
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
    int special = 0;

    // information about the Hamiltonian
    system_info<T, DIM> systemInfo;

    // Objects required to successfully calculate the conductivity
		Eigen::Array<std::complex<T>, -1, -1> Gamma0;
		Eigen::Array<std::complex<T>, -1, -1> Gamma1;
		Eigen::Array<std::complex<T>, -1, -1> Gamma2;
		Eigen::Array<std::complex<T>, -1, -1> Gamma3;

	  std::string dirName;


    conductivity_nonlinear(system_info<T, DIM>&);
		void read();
    void calculate();
	
};

template <typename T, unsigned DIM>
conductivity_nonlinear<T, DIM>::conductivity_nonlinear(system_info<T, DIM>& info){
  std::string name = info.filename;
	file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

  // retrieve the information about the Hamiltonian
  systemInfo = info;

  // location of the information about the conductivity
  dirName = "/Calculation/conductivity_nonlinear/";
  
  // check whether the conductivity_nonlinear was asked for
  try{
    H5::Exception::dontPrint();
    get_hdf5(&direction, &file, (char*)(dirName+"Direction").c_str());									
    isRequired = true;
  } catch(H5::Exception& e){}
  

}
	

template <typename T, unsigned DIM>
void conductivity_nonlinear<T, DIM>::read(){
	debug_message("Entered conductivity_nonlinear::read.\n");
	//This function reads all the data from the hdf5 file that's needed to 
  //calculate the nonlinear conductivity
	 

  // Check if the data for the nonlinear conductivity exists
  if(!isRequired){
    std::cout << "Data for nonlinear conductivity does not exist. Exiting.\n";
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


  bool hasGamma0 = false; 
  bool hasGamma1 = false;
  bool hasGamma2 = false; 
  bool hasGamma3 = false;


  // Retrieve the Gamma0 Matrix. This matrix is not needed for the special case
  std::string MatrixName = dirName + "Gamma0" + dirString;
  if(!special){
    try{
      verbose_message("Filling the Gamma0 matrix.\n");
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





  // Retrieve the Gamma1 Matrix
  MatrixName = dirName + "Gamma1" + dirString;
  try{
		verbose_message("Filling the Gamma1 matrix.\n");
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
  std::string MatrixName = dirName + "Gamma2" + dirString;
  try{
		verbose_message("Filling the Gamma2 matrix.\n");
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
    std::string MatrixName = dirName + "Gamma3" + dirString;
    try{
      verbose_message("Filling the Gamma3 matrix.\n");
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

  // 1/kT, where k is the Boltzmann constant in eV/K
  U beta = 1.0/8.6173303*pow(10,5)/temperature;

  
  U e_fermi = 0.1;
  std::cout << "Using default value for the Fermi energy: " << e_fermi << " in the KPM scale [-1, 1].\n";

  U scat = 0.01;
  //U scat = 0.0032679;
  std::cout << "Using default value for the scattering broadening: " << scat << "in the KPM scale [-1,1].\n"; 

	
	// Calculate the number of frequencies and energies needed to perform the calculation.
	int N_energies = NumPoints;
  double lim = 0.99;
  Eigen::Matrix<U, -1, 1> energies;
  energies  = Eigen::Matrix<U, -1, 1>::LinSpaced(N_energies, -lim, lim);

	int N_omegas = 200;
  double minFreq = 0.01;
  double maxFreq = 1.5;
  Eigen::Matrix<U, -1, 1> frequencies;
  Eigen::Matrix<U, -1, 1> frequencies2;
  frequencies = Eigen::Matrix<U, -1, 1>::LinSpaced(N_omegas, minFreq, maxFreq);
  frequencies2 = Eigen::Matrix<U, -1, 1>::Zero(1,1);
  std::cout << "Using default range of frequencies: " << N_omegas << " points from " << minFreq << " to " << maxFreq;
  std::cout << " in the KPM scale [-1,1]\n";


	std::complex<U> imaginary(0.0, 1.0);

  
  // Functions that are going to be used by the contractor
  int NumMoments1 = NumMoments;
  std::function<U(int, U)> deltaF = [beta, e_fermi, NumMoments1](int n, U energy)->U{
    return delta(n, energy)*U(1.0/(1.0 + U(n==0)))*fermi_function(energy, e_fermi, beta)*kernel_jackson<U>(n, NumMoments1);
  };


  Eigen::Matrix<std::complex<U>, -1, 1> temp1; 
  Eigen::Matrix<std::complex<U>, -1, 1> temp2; 
  Eigen::Matrix<std::complex<U>, -1, 1> temp3; 
  Eigen::Matrix<std::complex<U>, -1, 1> temp4; 
  Eigen::Matrix<std::complex<U>, -1, 1> cond; 

  temp1 = Eigen::Matrix<std::complex<U>, -1, 1>::Zero(N_omegas, 1);
  temp2 = Eigen::Matrix<std::complex<U>, -1, 1>::Zero(N_omegas, 1);
  cond  = Eigen::Matrix<std::complex<U>, -1, 1>::Zero(N_omegas, 1);


  temp1 = contract2<U>(deltaF, 0, greenAscat<U>(scat), NumMoments, Gamma, energies, -frequencies);
  temp2 = contract2<U>(deltaF, 1, greenRscat<U>(scat), NumMoments, Gamma, energies, frequencies);
  temp3 = contract2<U>(deltaF, 0, greenAscat<U>(scat), NumMoments, Gamma, energies, -frequencies2);
  temp4 = contract2<U>(deltaF, 1, greenRscat<U>(scat), NumMoments, Gamma, energies, frequencies2);

  std::complex<U> freq;
  for(int i = 0; i < N_omegas; i++){
    freq = std::complex<U>(frequencies(i), scat);  
    cond(i) += (temp1(i) + temp2(i) + temp3(0) + temp4(0)/freq/freq;
  }
  cond *= -imaginary*U(4.0*systemInfo.num_orbitals*systemInfo.spin_degeneracy/systemInfo.unit_cell_area);

	
  
  //Output to a file
  std::ofstream myfile;
  myfile.open ("nonlinear_cond.dat");
  for(int i=0; i < N_omegas; i++){
    freq = std::complex<U>(frequencies(i), scat);  
    myfile   << frequencies.real()(i)*systemInfo.energy_scale << " " << cond.real()(i) << " " << cond.imag()(i) << "\n";
  }
  myfile.close();
  debug_message("Left calc_nonlinear_cond.\n");



};
