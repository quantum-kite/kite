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
#include "arpes.hpp"

#include "functions.hpp"
#include "macros.hpp"

template <typename T, unsigned DIM>
arpes<T, DIM>::arpes(system_info<T, DIM>& sysinfo, shell_input & vari){
    // Class constructor
    
    systemInfo  = &sysinfo;              // retrieve the information about the Hamiltonian
    variables   = vari;                   // retrieve the shell input
    dirName     = "/Calculation/arpes/";     // location of the information about ARPES
    
    isRequired = is_required() && variables.ARPES_is_required;         // was ARPES requested?
    isPossible = false;                 // do we have all we need to calculate ARPES?
    if(isRequired){
        set_default_parameters();
        isPossible = fetch_parameters();
        override_parameters();

      if(isPossible){
          printARPES();                  // Print all the parameters used
          calculate();
      } else {
        std::cout << "ERROR. ARPES was requested but the data "
            "needed for its computation was not found in the input .h5 file. "
            "Make sure KITEx has processed the file first. Exiting.";
        exit(1);
      }
    }
}

template <typename T, unsigned DIM>
void arpes<T, DIM>::printARPES(){
  double scale = systemInfo->energy_scale;
  double shift = systemInfo->energy_shift;
    std::cout << "ARPES will be calculated with the following parameters:\n"
        "   Min energy: "           << Emin*scale + shift << ((default_energies)?" (default)":"") << "\n"
        "   Max energy: "           << Emax*scale + shift << ((default_energies)?" (default)":"") << "\n"
        "   Number of energies: "   << NumEnergies        << ((default_energies)?" (default)":"") << "\n"
        "   Number of moments: "    << NumMoments         << ((default_NumMoments)?       " (default)":"") << "\n";

        if(calculate_full_arpes){
          std::cout << 
        "   Incident vector: "      << incident_vector.transpose()    << ((default_incident)?" (default)":"") << "\n"
        "   Incident frequency: "   << freq               << ((default_freq)?" (default)":"") << "\n"
        "   Fermi energy: "         << fermi               << ((default_fermi)?" (default)":"") << "\n";
        } else {
          std::cout << "Calculating spectral function instead of ARPES\n";
        }

          std::cout << 
        "   Number of k-vectors: "  << NumVectors         << "\n"
        "   Filename: "             << filename           << ".dat" << ((default_filename)?" (default)":"") << "\n"
        "   Kernel: "               << kernel           << ((default_kernel)?           " (default)":"") << "\n";
    if(kernel == "green"){
        std::cout << "   Kernel parameter: "     << kernel_parameter*scale << ((default_kernel_parameter)? " (default)":"") << "\n";
    }
}

template <typename T, unsigned DIM>
bool arpes<T, DIM>::is_required(){
    // check whether the local density of states was asked for
    // if this quantity exists, so should all the others.

    name = systemInfo->filename;
	H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);
    bool result = false;
    try{
        H5::Exception::dontPrint();
        get_hdf5(&NumMoments, &file, (char*)(dirName+"NumMoments").c_str());									
        result = true;
    } catch(H5::Exception& e){}
  

    file.close();

    return result;
}
	
template <typename T, unsigned DIM>
void arpes<T, DIM>::override_parameters(){
  debug_message("Entered override_parameters.\n");
  double scale = systemInfo->energy_scale;
  double shift = systemInfo->energy_shift;

    if(variables.ARPES_Name != ""){
        filename          = variables.ARPES_Name;
        default_filename  = false;
    }

    if(variables.ARPES_Temp != -8888){
        temperature       = variables.ARPES_Temp;
        beta              = 1.0/8.6173303*pow(10,5)/temperature;
        default_temp      = false;
    }

    if(variables.ARPES_freq != -8888){
        freq              = (variables.ARPES_freq - shift)/scale;
        default_freq      = false;
    }

    if(variables.ARPES_Fermi != -8888){
        fermi          = (variables.ARPES_Fermi - shift)/scale;
        default_fermi  = false;
    }

    if(variables.ARPES_vec.matrix().norm() != 0.0){
        incident_vector  = variables.ARPES_vec;
        default_incident = false;
    }

    if(variables.ARPES_Emin != 8888){
        Emin              = (variables.ARPES_Emin - shift)/scale;
        default_energies  = false;
    }
    if(variables.ARPES_Emax != -8888){
        Emax              = (variables.ARPES_Emax - shift)/scale;
        default_energies  = false;
    }
    if(variables.ARPES_NumEnergies != -1){
        NumEnergies      = variables.ARPES_NumEnergies;
        default_energies = false;
    }
    if(variables.ARPES_NumMoments != -1){
        NumMoments      = variables.ARPES_NumMoments;
        default_NumMoments = false;
    }

    if(variables.ARPES_kernel != ""){
        kernel         = variables.ARPES_kernel;
        default_kernel = false;
    }
    if(kernel == "green"){
      if(variables.ARPES_kernel_parameter != -8888.8){
        kernel_parameter = variables.ARPES_kernel_parameter/systemInfo->energy_scale;
        default_kernel_parameter = false;
      }
    }

    if(default_energies == false)
      energies = Eigen::Matrix<float, -1, 1>::LinSpaced(NumEnergies, Emin, Emax);

    if(variables.ARPES_calculate_full_arpes == false){
      calculate_full_arpes = false;
    }


  debug_message("Left override_parameters.\n");
}

template <typename T, unsigned DIM>
void arpes<T, DIM>::set_default_parameters(){
  /* Sets default values for all the parameters needed for the computation of ARPES */
  debug_message("Entered set_default_parameters.\n");

  // Name of the output file created by KITE-tools
  filename = "arpes";
  name = systemInfo->filename;
  default_filename = true;

  // KPM scaling parameters
  double scale = systemInfo->energy_scale;
  double shift = systemInfo->energy_shift;

  // default values for the range of energies in which to evaluate the spectral function
  Emin              = -1.0;
  Emax              = 1.0;
  NumEnergies       = 1000;
  default_energies  = true;
  energies = Eigen::Matrix<float, -1, 1>::LinSpaced(NumEnergies, Emin, Emax);

  // default temperature
  temperature       = 0.1/scale;
  beta              = 1.0/8.6173303*pow(10,5)/temperature;
  default_temp      = true;

  // default Fermi energy of the system
  fermi             = (0.0 - shift)/scale;
  default_fermi     = true;

  // default momentum of the incident wave
  incident_vector = Eigen::Array<double, -1, -1>::Zero(DIM, 1);
  if(DIM == 2) incident_vector << 1.0, 0.0; 
  if(DIM == 3) incident_vector << 1.0, 0.0, 0.0; 
  default_incident  = true;

  // default frequency of the incident wave
  freq          = 0.0;
  default_freq  = true;

  // The user may want to only calculate the spectral function
  calculate_full_arpes = true;

  kernel = "jackson";
  default_kernel = true;
  default_kernel_parameter = true;
}


template <typename T, unsigned DIM>
bool arpes<T, DIM>::fetch_parameters(){
	/* Reads all the data from the hdf5 file that's needed to calculate ARPES */
	debug_message("Entered arpes::fetch_parameters.\n");
	 
  // Check if the data for the arpes exists
  if(!isRequired){
    std::cout << "Data for ARPES does not exist. Exiting.\n";
    exit(1);
  }

  H5::DataSet * dataset;
  H5::DataSpace * dataspace;
  
  hsize_t dim[2];
  
  H5::H5File file = H5::H5File(name.c_str(), H5F_ACC_RDONLY);
  

  dataset            = new H5::DataSet(file.openDataSet("/Calculation/arpes/k_vector")  );
  dataspace          = new H5::DataSpace(dataset->getSpace());
  dataspace -> getSimpleExtentDims(dim, NULL);
  dataspace->close(); delete dataspace;
  dataset->close();   delete dataset;
  NumVectors = dim[0];

  
  arpes_k_vectors = Eigen::Matrix<double, -1, -1>::Zero(DIM,NumVectors);
  //std::cout << "Energies: " << Emin << " " << Emax << "\n";
  energies = Eigen::Matrix<float, -1, 1>::LinSpaced(NumEnergies, Emin, Emax);

  
   //Fetch the momentum vectors from the hdf file
  get_hdf5(&NumMoments, &file, (char*)(dirName+"NumMoments").c_str());	
  get_hdf5(arpes_k_vectors.data(), &file, (char*)"/Calculation/arpes/k_vector");
  
  // Check whether the matrices we're going to retrieve are complex or not
  int complex = systemInfo->isComplex;

  bool result = false;
  // Retrieve the kmu Matrix
  std::string MatrixName = dirName + "kMU";
  try{
    debug_message("Filling the kMU matrix.\n");
    kMU = Eigen::Matrix<std::complex<T>,-1,-1>::Zero(NumMoments, NumVectors);
    
    if(complex)
      get_hdf5(kMU.data(), &file, (char*)MatrixName.c_str());
    if(!complex){
      Eigen::Matrix<T,-1,-1> kMUReal; 
      kMUReal = Eigen::Matrix<T,-1,-1>::Zero(NumMoments, NumVectors); 
      get_hdf5(kMUReal.data(), &file, (char*)MatrixName.c_str()); 

      
      kMU = kMUReal.template cast<std::complex<T>>();
    }				

    result = true;
  } catch(H5::Exception& e) {debug_message("ARPES: There is no kMU matrix.\n");}

  file.close();
  debug_message("Left ARPES::fetch_parameters.\n");
  
  return result;
}

template <typename U, unsigned DIM>
void arpes<U, DIM>::calculate(){
  debug_message("Entered arpes::calculate\n");
    

  Eigen::Matrix<std::complex<U>, -1, -1> ARPES = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(NumEnergies, NumVectors);
  Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor> OrderedMU;
  OrderedMU = kMU;

  
  debug_message("starting parallelization\n");
  omp_set_num_threads(systemInfo->NumThreads);
#pragma omp parallel 
{
  int localN = NumMoments/systemInfo->NumThreads;
  int thread_id = omp_get_thread_num();
  long offset = thread_id*localN*NumVectors;
  Eigen::Map<Eigen::Matrix<std::complex<U>, -1, -1, Eigen::RowMajor>> localkMU(OrderedMU.data() + offset, localN, NumVectors);

    


    

  // Calculate the part that depends on the energies
  Eigen::Matrix<std::complex<U>, -1, -1> GammaE;
  GammaE = Eigen::Matrix<std::complex<U>, -1, -1>::Zero(NumEnergies, localN);
  U factor, kern, ferm, del;
  if(kernel =="jackson"){
      for(int m = 0; m < localN; m++){
        factor = 1.0/(1.0 + U((m + thread_id*localN)==0));
        kern   = kernel_jackson<U>(m + thread_id*localN, NumMoments)*factor;
        for(int i = 0; i < NumEnergies; i++){
          
          ferm = 1.0;
          if(calculate_full_arpes){
            ferm      = fermi_function(energies(i) - freq, fermi, beta);
          }

          del         = delta(m + thread_id*localN, energies(i) - freq);
          GammaE(i,m) = del*ferm*kern; 
        }
      }
  }

  if(kernel == "green"){
      std::complex<U> c_energy;
      for(int m = 0; m < localN; m++){
        factor = 1.0/(1.0 + U((m + thread_id*localN)==0));
        for(int i = 0; i < NumEnergies; i++){

          ferm = 1.0;
          if(calculate_full_arpes){
            ferm      = fermi_function(energies(i) - freq, fermi, beta);
          }

          c_energy    = std::complex<U>(energies(i) - freq, kernel_parameter);
          del         = -green<std::complex<U>>(m + thread_id*localN, 1, c_energy).imag()/M_PI;
          GammaE(i,m) = del*factor;
        }
      }
  }
    
#pragma omp critical
  ARPES += GammaE*localkMU;
}

  // Save the density of states to a file
  double scale = systemInfo->energy_scale;
  float shift = systemInfo->energy_shift;

  // Generate a vector which which to shift the values of the energies
  Eigen::Matrix<float, -1, -1> shifts = Eigen::Matrix<float, -1, -1>::Zero(NumEnergies, 1);
  for(int ii = 0; ii < NumEnergies; ii++){
    shifts(ii) = shift;
  }


  if(calculate_full_arpes){

    Eigen::Matrix<double, -1, -1> incidents = Eigen::Matrix<double, -1, -1>::Zero(DIM, NumVectors);
    for(int ii = 0; ii < NumVectors; ii++){
      incidents.col(ii) = incident_vector;
    }

    Eigen::Matrix<double, -1, -1> modulation = Eigen::Matrix<double, -1, -1>::Zero(1, NumVectors);
    modulation = incident_vector.matrix().transpose()*(incidents + arpes_k_vectors.matrix())/incident_vector.matrix().norm();


    for(int ii = 0; ii < NumVectors; ii++){
      ARPES.col(ii) *= modulation(ii)*modulation(ii);
    }
  }



  std::ofstream myfile;
  myfile.open(filename + ".dat");
  myfile << "k-vectors:\n";
  for(unsigned int i = 0; i < NumVectors; i++)  myfile << arpes_k_vectors(0,i) << " " << arpes_k_vectors(1,i) << "\n";
  myfile << "Energies:\n";
  myfile << energies*scale + shifts << "\n";
  myfile << "ARPES:\n";
  myfile << ARPES.real() << "\n";
  myfile.close();
  debug_message("Left arpes::calculate\n");
}




// Instantiations
template class arpes<float, 1u>;
template class arpes<float, 2u>;
template class arpes<float, 3u>;

template class arpes<double, 1u>;
template class arpes<double, 2u>;
template class arpes<double, 3u>;

template class arpes<long double, 1u>;
template class arpes<long double, 2u>;
template class arpes<long double, 3u>;
