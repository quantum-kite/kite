/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/
#include "Generic.hpp"
#include "ComplexTraits.hpp"
#include "Global.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
#include "myHDF5.hpp"
#include "Random.hpp"
template <typename T, unsigned D>
class Hamiltonian;
template <typename T, unsigned D>
class KPM_Vector;
#include "queue.hpp"
#include "Simulation.hpp"
#include "SimulationGlobal.hpp"
#include "Hamiltonian.hpp"
template <typename T,unsigned D>
GlobalSimulation<T,D>::GlobalSimulation( char *name ) : rglobal(name)
{
  debug_message("Entered global_simulation\n");
  Global.ghosts.resize( rglobal.get_BorderSize() );
  std::fill(Global.ghosts.begin(), Global.ghosts.end(), 0);
    
  // Regular quantities to calculate, such as DOS and CondXX
  H5::H5File * file         = new H5::H5File(name, H5F_ACC_RDONLY);

  // Fetch the energy scale and the magnetic field, if it exists
  get_hdf5<double>(&EnergyScale,  file, (char *)   "/EnergyScale");
  delete file;

    
    
    
    
  // This function reads the h5 configuration file and checks which regular
  // functions need to be calculated. Then, it places the requests in a queue
  std::vector<measurement_queue> queue = fill_queue(name); 
  std::vector<singleshot_measurement_queue> ss_queue = fill_singleshot_queue(name);

	
    
  omp_set_num_threads(rglobal.n_threads);
  debug_message("Starting parallelization\n");
  int calculate_wavepacket = 0;
#pragma omp parallel default(shared)
  {
    Simulation<T,D> simul(name, Global);
      
    // Measure the average time it takes to run a multiplication
    // This will allow us to obtain an estimate for the time it'll take
    // for the program to run

    // Maybe put this inside an #if statement??
#pragma omp master
    {
      if(ESTIMATE_TIME == 1){
        std::cout << "------------------- TIME ESTIMATE --------------------\n";
        std::cout << "Estimate Chebyshev recursion time.\n";
        std::cout << "To disable this feature set the flag ESTIMATE_TIME=0.\n";
      }
    }
#pragma omp barrier
      
    if(ESTIMATE_TIME == 1){
      Global.kpm_iteration_time = simul.time_kpm(100);
    }
      
#pragma omp master 
    {
      if(ESTIMATE_TIME == 1){

        double queue_time = 0;
        double ss_queue_time = 0;
        // obtain the times for the singlehsot queue
        for(unsigned int i = 0; i < ss_queue.size(); i++){
          ss_queue.at(i).embed_time(Global.kpm_iteration_time);
          ss_queue_time += ss_queue.at(i).time_length;
        }
          
        // obtain the times for the normal queue
        for(unsigned int i = 0; i < queue.size(); i++){
          queue.at(i).embed_time(Global.kpm_iteration_time);
          queue_time += queue.at(i).time_length;
        }

        std::cout << "Estimated run time: ";
        std::cout << print_time(queue_time + ss_queue_time);
        std::cout << "\n------------------------------------------------------\n\n";
      }
    }
#pragma omp barrier

    verbose_message("-------------------------- CALCULATIONS --------------------------\n");
    // execute the singleshot queue
    for(unsigned int i = 0; i < ss_queue.size(); i++){
      verbose_message("Calculating SingleShot. This will take around ");
      verbose_message(print_time(ss_queue.at(i).time_length));
      verbose_message("\n");
      simul.Single_Shot(EnergyScale, ss_queue.at(i)); 
    }
      
    // execute the regular queue
    for(unsigned int i = 0; i < queue.size(); i++){
      verbose_message("Calculating ");
      verbose_message(queue.at(i).label);
      verbose_message(". This will take around ");
      verbose_message(print_time(queue.at(i).time_length));
      verbose_message("\n");
      simul.Measure_Gamma(queue.at(i));			
    }
    verbose_message("------------------------------------------------------------------\n\n");
#pragma omp barrier



    // Check if the Gaussian_Wave_Packet needs to be calculated
    //H5::Exception::dontPrint();
#pragma omp master
    {
      try{
        H5::H5File * file = new H5::H5File(name, H5F_ACC_RDONLY);
        int dummy_var;
        get_hdf5<int>(&dummy_var, file, (char *) "/Calculation/gaussian_wave_packet/NumDisorder");
        file->close();  
        delete file;
        calculate_wavepacket = 1;
      } catch(H5::Exception& e) {debug_message("Wavepacket: no need to calculate.\n");}
    }
      
    // Now calculate it
#pragma omp barrier
    if(calculate_wavepacket)
      simul.Gaussian_Wave_Packet();
  }
  debug_message("Left global_simulation\n");
};

template class GlobalSimulation<float ,1u>;
template class GlobalSimulation<double ,1u>;
template class GlobalSimulation<long double ,1u>;
template class GlobalSimulation<std::complex<float> ,1u>;
template class GlobalSimulation<std::complex<double> ,1u>;
template class GlobalSimulation<std::complex<long double> ,1u>;

template class GlobalSimulation<float ,3u>;
template class GlobalSimulation<double ,3u>;
template class GlobalSimulation<long double ,3u>;
template class GlobalSimulation<std::complex<float> ,3u>;
template class GlobalSimulation<std::complex<double> ,3u>;
template class GlobalSimulation<std::complex<long double> ,3u>;

template class GlobalSimulation<float ,2u>;
template class GlobalSimulation<double ,2u>;
template class GlobalSimulation<long double ,2u>;
template class GlobalSimulation<std::complex<float> ,2u>;
template class GlobalSimulation<std::complex<double> ,2u>;
template class GlobalSimulation<std::complex<long double> ,2u>;
