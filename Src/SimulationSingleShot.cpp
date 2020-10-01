/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2020, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/


#include "Generic.hpp"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"
#include "Global.hpp"
#include "Random.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
template <typename T, unsigned D>
class Hamiltonian;
template <typename T, unsigned D>
class KPM_Vector;
#include "queue.hpp"
#include "Simulation.hpp"
#include "Hamiltonian.hpp"
#include "KPM_VectorBasis.hpp"
#include "KPM_Vector.hpp"

std::complex<double> green(int n, int sigma, std::complex<double> energy){
  const std::complex<double> i(0.0,1.0); 
  std::complex<double> sq = sqrt(1.0 - energy*energy);
  return 2.0*sigma/sq*i*exp(-sigma*n*1.0*acos(energy)*i);
}

template <typename T, unsigned D>
void Simulation<T,D>::calc_singleshot() {
  Eigen::Array<double, -1, 1> energies;
  Eigen::Array<double, -1, 1> gammas;
  Eigen::Array<int, -1, 1> preserve_disorders;
  Eigen::Array<int, -1, 1> moments;
  int NDisorder, NRandom, direction;
  std::string direction_string;

    // Make sure that all the threads are ready before opening any files
    // Some threads could still be inside the Simulation constructor
    // This barrier is essential
#pragma omp barrier
  int calculate_singleshot_local = false;
#pragma omp master
{
  H5::H5File * file1 = new H5::H5File(name, H5F_ACC_RDONLY);
  Global.calculate_singleshot = false;
  try{
    debug_message("single_shot dc checking if we need to calculate it.\n");
    get_hdf5<int>(&direction, file1, (char *)   "/Calculation/singleshot_conductivity_dc/Direction");
    Global.calculate_singleshot = true;
    
  } catch(H5::Exception& e) {debug_message("singleshot dc: no need to calculate it.\n");}
  file1->close();
  delete file1;
  
}
#pragma omp barrier
#pragma omp critical
  calculate_singleshot_local = Global.calculate_singleshot;
#pragma omp barrier


  if(calculate_singleshot_local){
#pragma omp master
      {
        std::cout << "Calculating SingleShot.\n";
      }
#pragma omp barrier
#pragma omp critical
{
      H5::H5File * file = new H5::H5File(name, H5F_ACC_RDONLY);
    get_hdf5<int>(&direction, file, (char *)   "/Calculation/singleshot_conductivity_dc/Direction");
    get_hdf5<int>(&NRandom, file, (char *)   "/Calculation/singleshot_conductivity_dc/NumRandoms");
    get_hdf5<int>(&NDisorder, file, (char *)   "/Calculation/singleshot_conductivity_dc/NumDisorder");
       
    if(direction == 0)
      direction_string = "x,x";
    else if(direction == 1)
      direction_string = "y,y";
    else if(direction == 2)
      direction_string = "z,z";
    else{
      std::cout << "Invalid singleshot direction. Has to be xx, yy or zz. Exiting.\n";
      exit(1);
    }
       
    // We also need to determine the number of energies that we need to calculate
    H5::DataSet * dataset_energy     	= new H5::DataSet(file->openDataSet("/Calculation/singleshot_conductivity_dc/Energy"));
    H5::DataSpace * dataspace_energy 	= new H5::DataSpace(dataset_energy->getSpace());
    hsize_t dims_out[2];		
    dataspace_energy->getSimpleExtentDims(dims_out, NULL);	
    energies = Eigen::Array<double, -1, 1>::Zero(dims_out[0]*dims_out[1], 1);	
    delete dataspace_energy;
    delete dataset_energy;
    get_hdf5<double>(energies.data(),  	file, (char *)   "/Calculation/singleshot_conductivity_dc/Energy");
      
    // We also need to determine the number of gammas that we need to calculate
    H5::DataSet * dataset_gamma     	= new H5::DataSet(file->openDataSet("/Calculation/singleshot_conductivity_dc/Gamma"));
    H5::DataSpace * dataspace_gamma 	= new H5::DataSpace(dataset_gamma->getSpace());
    dataspace_gamma->getSimpleExtentDims(dims_out, NULL);	
    gammas = Eigen::Array<double, -1, 1>::Zero(dims_out[0]*dims_out[1], 1);	
    delete dataspace_gamma;
    delete dataset_gamma;
    get_hdf5<double>(gammas.data(),  	file, (char *)   "/Calculation/singleshot_conductivity_dc/Gamma");

    // We also need to determine the number of preserve disorders that we need to calculate
    H5::DataSet * dataset_preserve_disorder     	= new H5::DataSet(file->openDataSet("/Calculation/singleshot_conductivity_dc/PreserveDisorder"));
    H5::DataSpace * dataspace_preserve_disorder 	= new H5::DataSpace(dataset_preserve_disorder->getSpace());
    dataspace_preserve_disorder->getSimpleExtentDims(dims_out, NULL);	
    preserve_disorders = Eigen::Array<int, -1, 1>::Zero(dims_out[0]*dims_out[1], 1);	
    delete dataspace_preserve_disorder;
    delete dataset_preserve_disorder;
    get_hdf5<int>(preserve_disorders.data(),  	file, (char *)   "/Calculation/singleshot_conductivity_dc/PreserveDisorder");

    // We also need to determine the number of moments that we need to calculate
    H5::DataSet * dataset_moments     	= new H5::DataSet(file->openDataSet("/Calculation/singleshot_conductivity_dc/NumMoments"));
    H5::DataSpace * dataspace_moments 	= new H5::DataSpace(dataset_moments->getSpace());
    dataspace_moments->getSimpleExtentDims(dims_out, NULL);	
    moments = Eigen::Array<int, -1, 1>::Zero(dims_out[0]*dims_out[1], 1);	
    delete dataspace_moments;
    delete dataset_moments;
    get_hdf5<int>(moments.data(),  	file, (char *)   "/Calculation/singleshot_conductivity_dc/NumMoments");

    file->close();
  delete file;
}

  singleshot(energies, gammas, preserve_disorders, moments, NDisorder, NRandom, direction_string);
  }
}

template <typename T, unsigned D>
void Simulation<T,D>::singleshot(Eigen::Array<double, -1, 1> energies,
  Eigen::Array<double, -1, 1> gammas,
  Eigen::Array<int, -1, 1> preserve_disorders,
  Eigen::Array<int, -1, 1> moments,
  int NDisorder, int NRandomV, const std::string & direction_string){
  // Calculate the longitudinal dc conductivity for a single value of the energy
    
  T tmp = assign_value(0.,0.);
  debug_message("Entered Single_Shot\n");
    
  // Obtain the relevant quantities from the queue
  std::string indices_string = direction_string;
  std::string name_dataset = "/Calculation/singleshot_conductivity_dc/SingleShot";
  int N_energies = energies.rows();
  double EnergyScale;

#pragma omp critical
{
  H5::H5File * fetchfile         = new H5::H5File(name, H5F_ACC_RDONLY);
  get_hdf5<double>(&EnergyScale,  fetchfile, (char *)   "/EnergyScale");
  fetchfile->close();
  delete fetchfile;
}
#pragma omp barrier
  double EScale = EnergyScale;


    
  // Fixing the factor
  double unit_cell_area = fabs(r.rLat.determinant());
  unsigned int number_of_orbitals = r.Orb; 	// This is necessary because the normalization factor inside the random 
  // vectors is not the number of lattice sites N but the number of states, 
  // which is N*number_of_orbitals
  unsigned int spin_degeneracy = 1;
    
  double factor = -2.0*spin_degeneracy*number_of_orbitals/unit_cell_area;	// This is in units of sigma_0, hence the 4
    
  // process the string with indices and verify if the demanded
  // calculation is meaningful. For that to be true, this has to be a 
  // longitudinal conductivity
  std::vector<std::vector<unsigned>> indices = process_string(indices_string);
  if(indices.at(0).at(0) != indices.at(1).at(0)){
    std::cout << "SingleShot is only meaningful for the longitudinal conductivity.";
    std::cout << "Please use directions 'x,x', 'y,y' or 'z,z'. Exiting.\n";
    exit(1);
  }


  // initialize the kpm vectors necessary for this calculation
  typedef typename extract_value_type<T>::value_type value_type;
	
  // initialize the conductivity array

#if (SSPRINT == 0)
  KPM_Vector<T,D> phi (2, *this);
#elif (SSPRINT != 0)
  KPM_Vector<T,D> phir1 (2, *this);
  KPM_Vector<T,D> phir2 (2, *this);
#endif


  KPM_Vector<T,D> phi0(1, *this);
  KPM_Vector<T,D> phi1(1, *this);
  // iteration over each energy

#if (SSPRINT != 0)
  KPM_Vector<T,D> phi2(1, *this);
#endif
  // iteration over disorder and the number of random vectors
  // initialize the conductivity array
  Eigen::Array<T, -1, -1> cond_array;
  cond_array = Eigen::Array<T, -1, -1>::Zero(1, N_energies);
#if (SSPRINT != 0)
#pragma omp master
  {
    std::cout << "Study of the convergence. To disable these messages"
      " set the flag SSPRINT to 0 in the main.cpp file\n";
  }
#pragma omp barrier 
#endif
  long average = 0;
  double job_energy, job_gamma, job_preserve_disorder;
  int job_NMoments;
  for(int disorder = 0; disorder < NDisorder; disorder++){
    h.generate_disorder();
    h.build_velocity(indices.at(0),0u);
    h.build_velocity(indices.at(1),1u);

    // iteration over each energy and gammma
    for(int job_index = 0; job_index < N_energies; job_index++){
      job_energy = energies(job_index);
      job_gamma = gammas(job_index);
      job_preserve_disorder = preserve_disorders(job_index);
      job_NMoments = moments(job_index);
      std::complex<double> energy(job_energy, job_gamma);
        
      if(job_preserve_disorder == 0){
        h.generate_disorder();
        h.build_velocity(indices.at(0),0u);
        h.build_velocity(indices.at(1),1u);
      }


      long average_R = average;
      // iteration over disorder and the number of random vectors
      for(int randV = 0; randV < NRandomV; randV++){

#if (SSPRINT == 0)
        debug_message("Started SingleShot calculation for SSPRINT=0\n");
        // initialize the random vector
        phi0.initiate_vector();					
        phi0.Exchange_Boundaries(); 	
        phi1.v.col(0).setZero();

            
        // calculate the left KPM vector
        phi.set_index(0);				
        generalized_velocity(&phi, &phi0, indices, 0);      // |phi> = v |phi_0>


        for(int n = 0; n < job_NMoments; n++){		
          if(n!=0) cheb_iteration(&phi, n-1);

          phi1.v.col(0) += phi.v.col(phi.get_index())
            *green(n, 1, energy).imag()/(1.0 + int(n==0));
        }
          
        // multiply phi1 by the velocity operator again. 
        // We need a temporary vector to mediate the operation, which will be |phi>
        phi.v.col(0) = phi1.v.col(0);
        phi.set_index(0);
        phi.Exchange_Boundaries();
        generalized_velocity(&phi1, &phi, indices, 1);
        phi1.empty_ghosts(0);
        
          
        phi.set_index(0);			
        phi.v.col(0) = phi0.v.col(0);
        phi0.v.col(0).setZero(); 

        for(int n = 0; n < job_NMoments; n++){		
          if(n!=0) cheb_iteration(&phi, n-1);

          phi0.v.col(0) += phi.v.col(phi.get_index())
            *green(n, 1, energy).imag()/(1.0 + int(n==0));
        }
          
          
        // finally, the dot product of phi1 and phi0 yields the conductivity
        tmp *= 0.;
        for(std::size_t ii = 0; ii < r.Sized ; ii += r.Ld[0])
          tmp += T(phi1.v.col(0).segment(ii,r.Ld[0]).adjoint() * phi0.v.col(0).segment(ii,r.Ld[0]));
        cond_array(job_index) += (tmp - cond_array(job_index))/value_type(average_R+1);						
        debug_message("Concluded SingleShot calculation for SSPRINT=0\n");
#elif (SSPRINT != 0)
#pragma omp master
        {
          std::cout << "   Random vector " << randV << "\n";
        }
#pragma omp barrier
        debug_message("Started SingleShot calculation for SSPRINT!=0\n");
        // initialize the random vector
        phi0.initiate_vector();					
        phi0.Exchange_Boundaries(); 	
        phi1.v.col(0).setZero();

        // chebyshev recursion vectors
        phir1.set_index(0);				
        phir2.set_index(0);				

        phir2.v.col(0) = phi0.v.col(0);
        generalized_velocity(&phir1, &phi0, indices, 0);      // |phi> = v |phi_0>
        // from here on, phi0 is free to be used elsewhere, it is no longer needed
        phi0.v.col(0).setZero();

        for(int nn = 0; nn < SSPRINT; nn++){

          for(int n = nn*job_NMoments/SSPRINT; n < job_NMoments/SSPRINT*(nn+1); n++){	
            if(n!=0) cheb_iteration(&phir1, n-1);

            phi1.v.col(0) += phir1.v.col(phir1.get_index())
              *green(n, 1, energy).imag()/(1.0 + int(n==0));
          }
            
          // multiply phi1 by the velocity operator again. 
          // We need a temporary vector to mediate the operation, which will be |phi>
          // If the SSPRINT flag is set to true, we are going to need the phi1 vector again
          // so the product of phi1 by the velocity is stored in phi2 instead
          phi2.set_index(0);
          generalized_velocity(&phi2, &phi1, indices, 1);
          phi2.empty_ghosts(0);
          
            
          for(int n = nn*job_NMoments/SSPRINT; n < job_NMoments/SSPRINT*(nn+1); n++){		
            if(n!=0) cheb_iteration(&phir2, n-1);

            phi0.v.col(0) += phir2.v.col(phir2.get_index())
              *green(n, 1, energy).imag()/(1.0 + int(n==0));
          }
           
          // This is the conductivity for a smaller number of chebyshev moments
          // if you want to add it to the average conductivity, you have yo wait
          // until all the moments have been summed. otherwise the result would be wrong
          T temp;
          temp *= 0.;
          for(std::size_t ii = 0; ii < r.Sized ; ii += r.Ld[0])
            temp += phi2.v.col(0).segment(ii,r.Ld[0]).adjoint() * phi0.v.col(0).segment(ii,r.Ld[0]);


          if(nn == SSPRINT-1){
            cond_array(job_index) += (temp - cond_array(job_index))/value_type(average_R+1);
          }
            
#pragma omp master
          {
            std::cout << "   energy: " << (energy*EScale).real() << " broadening: "
                      << (energy*EScale).imag() << " moments: "; 
            std::cout << job_NMoments/SSPRINT*(nn+1) << " SS_Cond: " << temp*factor*(1.0*omp_get_num_threads()) << "\n" << std::flush;
            if(nn == SSPRINT-1)
              std::cout << "\n";
          }
#pragma omp barrier
        }
#endif
        average_R++;
        debug_message("Concluded SingleShot calculation for SSPRINT!=0\n");
      }
#if (SSPRINT!=0)
#pragma omp master
      {
        std::cout << "Average over " << NRandomV << " random vectors: " << cond_array(job_index)*factor*(1.0*omp_get_num_threads()) << "\n\n";
      }
#pragma omp barrier
#endif
    }
    average += NRandomV;
  }
  // finished calculating the longitudinal DC conductivity for all the energies
  // Now let's store the gamma matrix. Now we're going to use the 
  // property that gamma is hermitian: gamma_nm=gamma_mn*
    
#pragma omp master
  { 
    Global.singleshot_cond = Eigen::Matrix<T, -1, -1> :: Zero(1, N_energies);
  }
#pragma omp barrier
    
    
  //std::cout << "IMPORTANT ! ! !:\n V is not hermitian. Make sure you take this into account\n";
  // in this case there's no problem. both V are anti-hermitic, so the minus signs cancel
#pragma omp critical
  {
    Global.singleshot_cond += cond_array;			
  }
#pragma omp barrier
    
    
#pragma omp master
  {
			
    Global.singleshot_cond *= factor;
      
    // Create array to store the data
    Eigen::Array<double, -1, -1> store_data;
    store_data = Eigen::Array<double, -1, -1>::Zero(4, energies.rows());
      
    for(int ener = 0; ener < N_energies; ener++)
      {
        store_data(0, ener) = energies(ener)*EScale;
        store_data(1, ener) = gammas(ener)*EScale;
        store_data(2, ener) = preserve_disorders(ener);
        store_data(3, ener) = Global.singleshot_cond.real()(ener);
      }
      
      
			
    H5::H5File * file = new H5::H5File(name, H5F_ACC_RDWR);
    write_hdf5(store_data, file, name_dataset);
    delete file;
      
    // make sure the global matrix is zeroed
    Global.singleshot_cond.setZero();
    debug_message("Left single_shot");
  }
#pragma omp barrier
}

template class Simulation<float ,1u>;
template class Simulation<double ,1u>;
template class Simulation<long double ,1u>;
template class Simulation<std::complex<float> ,1u>;
template class Simulation<std::complex<double> ,1u>;
template class Simulation<std::complex<long double> ,1u>;

template class Simulation<float ,3u>;
template class Simulation<double ,3u>;
template class Simulation<long double ,3u>;
template class Simulation<std::complex<float> ,3u>;
template class Simulation<std::complex<double> ,3u>;
template class Simulation<std::complex<long double> ,3u>;

template class Simulation<float ,2u>;
template class Simulation<double ,2u>;
template class Simulation<long double ,2u>;
template class Simulation<std::complex<float> ,2u>;
template class Simulation<std::complex<double> ,2u>;
template class Simulation<std::complex<long double> ,2u>;
