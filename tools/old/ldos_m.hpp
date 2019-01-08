
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
class ldos{
	H5::H5File file;
	public:


    // Flags that will be accessed from the other classes
    bool isRequired = false;                        // was this quantity local density of states asked for?
    bool isPossible = false;                        // do we have all we need to calculate the density of states?

    // Variables needed to compute the Local density of states
    int NumMoments = -1;                                    // Number of Chebyshev moments
    int NumPositions = -1;                                  // Number of lattice sites in which to compute the LDoS
    int NumEnergies = -1;                                   // Number of energies
    std::string filename = "ldos";                          // Saving results to file with this name
    Eigen::Matrix<unsigned long, -1, -1> ldos_Orbitals;     // Position of the lattice sites
    Eigen::Matrix<float, -1, -1> energies;                  // Energies specified to be calculated
    std::string dirName;                                    // Name of the hdf5 dataset where the ldos moments are saved

    // Aditional variables
    system_info<T, DIM> *systemInfo;            // information about the Hamiltonian
    shell_input variables;                      // Input from the shell to override the configuration file

    // Objects required to successfully calculate the conductivity
    Eigen::Matrix<std::complex<T>, -1, -1> lMU;

    std::string name;

    // Class methods
    ldos(system_info<T, DIM>&, shell_input &);  // Constructor
	void fetch_parameters();                    // Get the parameters from the hdf5 file
	void override_parameters();                 // If shell variables were given, this function overrides the current parameters
    void calculate();                           // Compute the local density of states
	
};

template <typename T, unsigned DIM>
ldos<T, DIM>::ldos(system_info<T, DIM>& sysinfo, shell_input & vari){
    // Class constructor

    systemInfo = &sysinfo;              // retrieve the information about the Hamiltonian
    variables = vari;                   // retrieve the shell input
    dirName = "/Calculation/ldos/";     // location of the information about the conductivity

    // check whether the local density of states was asked for
    // if this quantity exists, so should all the others.
    name = sysinfo.filename;
	file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);
    try{
        H5::Exception::dontPrint();
        get_hdf5(&NumMoments, &file, (char*)(dirName+"NumMoments").c_str());									
        isRequired = true;
    } catch(H5::Exception& e){}
  

    //file.close();
}
	
template <typename T, unsigned DIM>
void ldos<T, DIM>::override_parameters(){
    if(variables.lDOS_Name != "")        filename    = variables.lDOS_Name;
};

template <typename T, unsigned DIM>
void ldos<T, DIM>::fetch_parameters(){
	debug_message("Entered ldos::fetch_parameters.\n");
	//This function reads all the data from the hdf5 file that's needed to 
    //calculate the LDoS
	 
    // Check if the data for the ldos exists
    if(!isRequired){
        std::cout << "Data for LDoS does not exist. Exiting.\n";
        exit(1);
    }
  
    H5::DataSet * dataset;
    H5::DataSpace * dataspace;
    hsize_t dim[1];

	//file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);
    dataset            = new H5::DataSet(file.openDataSet("/Calculation/ldos/Orbitals")  );
    dataspace          = new H5::DataSpace(dataset->getSpace());
    dataspace -> getSimpleExtentDims(dim, NULL);
    dataspace->close(); delete dataspace;
    dataset->close();   delete dataset;
    NumPositions = dim[0];

    dataset            = new H5::DataSet(file.openDataSet("/Calculation/ldos/Energy")  );
    dataspace          = new H5::DataSpace(dataset->getSpace());
    dataspace -> getSimpleExtentDims(dim, NULL);
    dataspace->close(); delete dataspace;
    dataset->close();   delete dataset;
    NumEnergies = dim[1];

    ldos_Orbitals = Eigen::Matrix<unsigned long, -1, -1>::Zero(NumPositions,1);
    energies = Eigen::Matrix<float, -1, -1>::Zero(NumEnergies,1);

     //Fetch the relevant parameters from the hdf file
    get_hdf5(&NumMoments, &file, (char*)(dirName+"NumMoments").c_str());	
    get_hdf5(ldos_Orbitals.data(), &file, (char*)"/Calculation/ldos/Orbitals");
    get_hdf5(energies.data(), &file, (char*)"/Calculation/ldos/Energy");


  // Check whether the matrices we're going to retrieve are complex or not
  int complex = systemInfo->isComplex;

  // Retrieve the lmu Matrix
  std::string MatrixName = dirName + "lMU";
  try{
		debug_message("Filling the lMU matrix.\n");
		lMU = Eigen::Matrix<std::complex<T>,-1,-1>::Zero(NumMoments, NumPositions);
		
		if(complex)
			get_hdf5(lMU.data(), &file, (char*)MatrixName.c_str());
		if(!complex){
			Eigen::Matrix<T,-1,-1> lMUReal; 
			lMUReal = Eigen::Matrix<T,-1,-1>::Zero(NumMoments, NumPositions); 
			get_hdf5(lMUReal.data(), &file, (char*)MatrixName.c_str()); 

			
            lMU = lMUReal.template cast<std::complex<T>>();
		}				

    isPossible = true;
  } catch(H5::Exception& e) {debug_message("DOS: There is no lMU matrix.\n");}
	
    file.close();
	debug_message("Left lDOS::fetch_parameters.\n");
}

template <typename U, unsigned DIM>
void ldos<U, DIM>::calculate(){
  if(!isPossible){
    std::cout << "Cannot retrieve lDOS since there is not enough information. Exiting.\n";
    exit(0);
  }

  verbose_message("  Number of energies: "); verbose_message(NumEnergies); verbose_message("\n");
  verbose_message("  Using kernel: Jackson\n");
  verbose_message("  File name: dos.dat\n");




  Eigen::Matrix<std::complex<U>, -1, -1> LDOS;
  LDOS = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(NumEnergies, NumPositions);
  omp_set_num_threads(systemInfo->NumThreads);

  Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> OrderedMU;
  OrderedMU = Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>::Zero(NumEnergies, NumPositions);
  OrderedMU = lMU;

#pragma omp parallel 
{
  int localN = NumMoments/systemInfo->NumThreads;
  int thread_id = omp_get_thread_num();
  long offset = thread_id*localN*NumPositions;
  Eigen::Map<Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>> locallMU(OrderedMU.data() + offset, localN, NumPositions);
  //std::cout << "thread_id: " << thread_id << "\n";
  //std::cout << "local LMU\n" << locallMU << "\n";
  // First perform the part of the product that only depends on the
  // chebyshev polynomial of the first kind
  Eigen::Matrix<std::complex<U>, -1, -1> GammaE;
  GammaE = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(NumEnergies, localN);


  U factor;
  for(int m = 0; m < localN; m++){
    factor = 1.0/(1.0 + U(m==0));
    for(int i = 0; i < NumEnergies; i++){
      GammaE(i,m) = delta(m + thread_id*localN, energies(i))*kernel_jackson<U>(m + thread_id*localN, NumMoments)*factor;
    }
  }

  Eigen::Matrix<std::complex<U>, -1, -1> localLDOS;
  localLDOS = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(NumEnergies, NumPositions);
  localLDOS = GammaE*locallMU;
#pragma omp critical
    LDOS += localLDOS;
}

  // Save the density of states to a file
  U mult = 1.0/systemInfo->energy_scale;
  std::ofstream myfile;
  for(unsigned pos = 0; pos < NumPositions; pos++){
    myfile.open(filename + std::to_string(pos) + ".dat");
    for(int i=0; i < NumEnergies; i++){
    myfile  << energies(i)*systemInfo->energy_scale << " " << LDOS(i).real()*mult << "\n";
    }
  myfile.close();
  }

  
};
