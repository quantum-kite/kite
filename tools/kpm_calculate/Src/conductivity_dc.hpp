#include "headers.hpp"

using std::cout;
using namespace H5;
using std::endl;

template <typename T, unsigned DIM>
class conductivity_dc{
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


    // information about the Hamiltonian
    system_info<T, DIM> systemInfo;

    // Objects required to successfully calculate the conductivity
		Eigen::Array<std::complex<T>, -1, -1> Gamma;

	  std::string dirName;


    conductivity_dc(system_info<T, DIM>&);
		void read();
    void calculate();
	
};

template <typename T, unsigned DIM>
conductivity_dc<T, DIM>::conductivity_dc(system_info<T, DIM>& info){
  std::string name = info.filename;
	file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);

  // retrieve the information about the Hamiltonian
  systemInfo = info;

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
void conductivity_dc<T, DIM>::read(){
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


  

  // Retrieve the Gamma Matrix
  std::string MatrixName = dirName + "Gamma" + dirString;
  try{
		verbose_message("Filling the Gamma matrix.\n");
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
  } catch(H5::Exception& e) {debug_message("Conductivity DC: There is no Gamma matrix.\n");}
	



	file.close();
	debug_message("Left conductivity_dc::read.\n");
}

template <typename U, unsigned DIM>
void conductivity_dc<U, DIM>::calculate(){
  if(!isPossible){
    std::cout << "Cannot calculate the DC conductivity because there is not enough information. Exiting.\n";
    exit(0);
  }

  U scat = 0.0032679;
  U beta = 1.0/8.6173303*pow(10,5)/temperature;
  U minEnergy, maxEnergy;
  //U beta = 20000.0;
  //int NEnergies = NumPoints;
  int NEnergies = 5000;
  Eigen::Matrix<U, -1, 1> energies;
  if(systemInfo.EnergyLimitsKnown){
    minEnergy = systemInfo.minEnergy;
    maxEnergy = systemInfo.maxEnergy;
  }
  else{
    minEnergy = -0.99;
    maxEnergy = 0.99;
    std::cout << "Using default limits for the integration in the energies.\n";
    std::cout << "For a more precise evaluation, calculate the density of states as well.\n";
  }
    minEnergy = -0.99;
    maxEnergy = 0.99;
  energies = Eigen::Matrix<U, -1, 1>::LinSpaced(NEnergies, minEnergy, maxEnergy);

  std::cout << "minEnergy:" << systemInfo.minEnergy << " maxEnergy:" << systemInfo.maxEnergy << "\n";


  // First perform the part of the product that only depends on the
  // chebyshev polynomial of the first kind
  Eigen::Array<std::complex<U>, -1, -1> GammaEN;
  GammaEN = Eigen::Array<std::complex<U>, -1, -1>::Zero(NEnergies, NumMoments);

  U factor;
  std::complex<U> complexEnergyP, complexEnergyN;
  for(int i = 0; i < NEnergies; i++){
    complexEnergyP = std::complex<U>(energies(i), scat);
    for(int n = 0; n < NumMoments; n++){
      for(int m = 0; m < NumMoments; m++){
        //factor = 1.0/(1.0 + U(m==0));
        factor = 1.0/(1.0 + U(m==0));
        //GammaEN(i,n) += Gamma(n,m)*delta(m,energies(i))*kernel_jackson<U>(m, NumMoments)*factor;
        GammaEN(i,n) += Gamma(n,m)*green(m, 1, complexEnergyP).imag()*factor;
      }
    }
  }


  // Now perform the part of the product that depends on both kinds of polynomials
  Eigen::Array<std::complex<U>, -1, -1> GammaE;
  GammaE = Eigen::Array<std::complex<U>, -1, -1>::Zero(NEnergies, 1);

  U energy;
  U den = systemInfo.num_orbitals*systemInfo.spin_degeneracy/systemInfo.unit_cell_area*4.0/2.0; 
  for(int i = 0; i < NEnergies; i++){
    complexEnergyP = std::complex<U>(energies(i), scat);
    complexEnergyN = std::complex<U>(energies(i), -scat);
    for(int n = 0; n < NumMoments; n++){
      factor = 1.0/(1.0 + U(n==0));
      GammaE(i) += (GammaEN(i,n)*dgreen<U>(n, 1, complexEnergyP) - std::conj(GammaEN(i,n))*dgreen<U>(n, -1, complexEnergyN))*factor*den;
      //GammaE(i) += GammaEN(i,n)*green(n, 1, energies(i))*factor*den;
    }
  }



  int NFermiEnergies = 1000;
  double minFermiEnergy = -1.0;
  double maxFermiEnergy = 1.0;
  std::cout << "Using default Fermi energies: " << NFermiEnergies << " points from ";
  std::cout << minFermiEnergy*systemInfo.energy_scale << " to " << maxFermiEnergy*systemInfo.energy_scale << ".\n";

  Eigen::Matrix<std::complex<U>, -1, 1> condDC;
  condDC = Eigen::Matrix<std::complex<U>, -1, 1>::Zero(NFermiEnergies, 1);

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
    condDC(i) = integrate(energies, integrand);
    }
  }


  std::ofstream myfile;
  myfile.open("condDC.dat");
  for(int i=0; i < NFermiEnergies; i++)
    myfile  << fermiEnergies(i)*systemInfo.energy_scale << " " << condDC.real()(i) << " " << condDC.imag()(i) << "\n";
  
  myfile.close();
  std::cout<<"-----Line "<<__LINE__<<"-----\n"<<std::flush;


};
