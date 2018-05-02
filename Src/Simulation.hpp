#ifndef _SIMULATION_HPP
#define _SIMULATION_HPP
#include "queue.hpp"


std::complex<double> green(int n, int sigma, std::complex<double> energy){
  const std::complex<double> i(0.0,1.0); 
  std::complex<double> sq = sqrt(1.0 - energy*energy);
  return pow(-1,n)*2.0*sigma/sq*i*exp(-sigma*n*1.0*acos(energy)*i);
}



template <typename T,unsigned D>
class GlobalSimulation {
private:
  GLOBAL_VARIABLES <T> Global;
  LatticeStructure <D> rglobal;
  
  
  // Regular quantities to calculate, such as DOS and CondXX
  Eigen::Array<double, -1, 1> singleshot_energies;
  double EnergyScale;

public:
  GlobalSimulation( char *name ) : rglobal(name)
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
#pragma omp parallel default(shared)
    {
      Simulation<T,D> simul(name, Global);
      
      // Measure the average time it takes to run a multiplication
      // This will allow us to obtain an estimate for the time it'll take
      // for the program to run
      Global.kpm_iteration_time = simul.time_kpm(100);
#pragma omp master 
      {
      verbose_message("On average, a KPM iteration takes ");
      verbose_message(Global.kpm_iteration_time);
      verbose_message(" seconds.\n");
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

      verbose_message("The entire calculation will take around ");
      verbose_message(print_time(queue_time + ss_queue_time));
      verbose_message("\n");
      }
#pragma omp barrier

      // execute the singleshot queue
      for(unsigned int i = 0; i < ss_queue.size(); i++){
        simul.Single_Shot(EnergyScale, ss_queue.at(i)); 
      }
      
      // execute the regular queue
			for(unsigned int i = 0; i < queue.size(); i++){
        simul.Measure_Gamma(queue.at(i));			
      }



    }
    
    debug_message("Left global_simulation\n");
  };
  
  
};
#endif


template <typename T,unsigned D>
class Simulation  {
public:
  KPMRandom <T>          rnd;
  std::vector<T>         ghosts;
  LatticeStructure <D>   r;      
  GLOBAL_VARIABLES <T> & Global;
  char                 * name;
  Hamiltonian<T,D>       h;
  Simulation(char *filename, GLOBAL_VARIABLES <T> & Global1): r(filename),  Global(Global1), name(filename), h(*this)  {
    rnd.init_random();
    ghosts.resize(Global.ghosts.size()/r.n_threads);
  };
  
  

  void Measure_Gamma(measurement_queue queue) {
    debug_message("Entered Measure_Gamma.\n");
    /* Calculates the Gamma matrix of parameters given by 'indices_string'. These matrices
       are used to calculate various properties of the quantum system, such as the density
       of states and the optical conductivity. The Gamma matrix is a multi-dimensional matrix 
       defined by:
		
       Gamma{i1, i2, ... , iN}(n1, n2, ..., nN) = < v^i1  T_n1(H)  v^i2  T_n2(H) ... v^iN  T_nN(H) >
			
       The first argument is the parameters and the second argument is the matrix indices. Here's 
       a quick rundown of what each of these objects means:
		 
       T_n(H) is the n-th order Chebyshev polynomial of the Hamiltonian matrix H
			
       Note that each i in the parameters actually stands for a set of letters, such as x, or xx, or xyyxxy.
       Unwinding those, v^{ij...mn} represents the  nested commutator of the position operator:
		
       v^{x} = [x,H],  v^{xy} = [x,[y,H]],  v^{yyy} = [y,[y,[y,H]]],  etc...
			
       Note that these velocity operators do not have the imaginary number in their definition, so they are not
       hermitian. When the parameter of v is empty, it is treated as the identity operator: v^{} = 1.		
       Some valid parameters for the Gamma matrix are: 
		
       ""   		  ->   G^{}(n) = < Tn(H) > 
       "x" 		  	->   G^{x}(n) = < v^x  Tn(H) >
       "xyyyxyx"	->   G^{xyyyxyx}(n) = < v^xyyyxyx  Tn(H) >
       "xy,y"		  ->	 G^{xy,y}(n,m) = < v^x  Tn(H)  v^y  Tm(H) >
       ","		  	->	 G^{,}(n,m) = < Tn(H)  Tm(H) >
       "y,,y,"		->	 G^{y,,y,}(n,m,p,q) = < v^y  Tn(H)  Tm(H)  v^y  Tp(H)  Tq(H) >
       etc...	  	->	 etc...
		
       These parameters should NOT contain anything else other than the characters ',', 'x' and 'y',
       not even whitespace characters. The dimension of the Gamma matrix is the number of commas + 1 */
		
		


    // Obtain the quantities needed for the simulation from the queue
	  int NRandomV = queue.NRandom;
    int NDisorder = queue.NDisorder;
    std::vector<int> N_moments = queue.NMoments;
    std::string indices_string = queue.direction_string;
    std::string name_dataset = queue.label;
 
    debug_message("Indices: "); debug_message(indices_string); debug_message(".\n");
		


    // First of all, we need to process the indices_string into something the program can use
    // Each element of this vector is a list of the indices of a generalized velocity operator
#pragma omp barrier
    std::vector<std::vector<int>> indices = process_string(indices_string);
    int dim = indices.size();
		
		
		
    // Check if the dimensions match
    if(dim != int(N_moments.size())){
      std::cout << "Dimension of the Gamma matrix does not match the number of chebyshev moments. Aborting.\n";
      exit(0);
    }
			
    // Determine the size of the gamma matrix we want to calculate
    int size_gamma = 1;
    for(int i = 0; i < dim; i++){
      if(N_moments.at(i) % 2 != 0){
	std::cout << "The number of moments must be an even number, due to limitations of the program. Aborting\n";
	exit(0);
      }
      size_gamma *= N_moments.at(i);
    }
		
    // Estimate of the time it'll take to run this function. 
    // It doesn't take into account parallelization or velocity matrix products
    /*debug_message("This will take around ");
      debug_message(size_gamma*Global.kpm_iteration_time);
      debug_message(" seconds\n");*/
		
		
		
    // Initialize the KPM vectors that will be needed to run the program 
    std::vector<KPM_Vector<T,D>*> kpm_vector(dim+1);
    kpm_vector.at(0) = new KPM_Vector<T,D> (1, *this);
    for(int i = 0; i < dim; i++)
		kpm_vector.at(i+1) = new KPM_Vector<T,D> (2, *this);
		
		//kpm_vector.push_back(KPM_Vector<T,D> (2, *this));
		
		
    // Define some pointers to make the code easier to read
    KPM_Vector<T,D> *kpm0 = kpm_vector.at(0);
    KPM_Vector<T,D> *kpm1 = kpm_vector.at(1);
    T * kpm0data = kpm0->v.col(0).data();
    T * kpm1data = kpm1->v.col(kpm1->get_index()).data();
    int axis0, axis1;
			
			
    // Make sure the local gamma matrix is zeroed
    Eigen::Array<T, -1, -1> gamma = Eigen::Array<T, -1, -1 >::Zero(1, size_gamma);

    long average = 0;
    for(int disorder = 0; disorder < NDisorder; disorder++){
      h.generate_disorder();
      for(int randV = 0; randV < NRandomV; randV++){
	
	
	kpm0->initiate_vector();			// original random vector
	kpm1->set_index(0);
	kpm1->v.col(0) = kpm0->v.col(0);
	kpm1->Exchange_Boundaries();
	// Check which generalized velocity operator needs to be calculated. 
	// This replaces the original random vector |0> by v|0> 
				
	switch(indices.at(0).size()){
	case 0:
	  verbose_message("case 0");
	  break;
	case 1:
						
	  verbose_message("case 1");
	  axis0 = indices.at(0).at(0);
	  kpm0->Velocity(kpm0data, kpm1data, axis0);  
	  kpm0->empty_ghosts(0);
	  kpm0->v.col(0) = -kpm0->v.col(0); // This minus sign is due to the fact that this Velocity operator is not self-adjoint
						
						
	  break;
					
	case 2:
	  verbose_message("case 2");
	  axis0 = indices.at(0).at(0);
	  axis1 = indices.at(0).at(1);
						
	  kpm0->Velocity2(kpm0data, kpm1data, axis0, axis1);  
	  kpm0->empty_ghosts(0);
	  break;
						
	case 3:
	  std::cout << "The matrix you're trying to calculate requires an operator that is not yet implemented: Velocity3.\n";
	  exit(0);
						
	default:
	  std::cout << "The matrix you're trying to calculate requires an operator that is not yet implemented.\n";
	  exit(0);
		
	}
					
	long index_gamma = 0;
	recursive_KPM(1, dim, N_moments, &average, &index_gamma, indices, &kpm_vector, &gamma);
				
				
	average++;
      }
    } 
		
		
    store_gamma(&gamma, N_moments, indices, name_dataset);
		
	// delete the kpm_vector
	delete kpm_vector.at(0);
    for(int i = 0; i < dim; i++)
		delete kpm_vector.at(i+1);
		
		
    debug_message("Left Measure_Gamma\n");
  }
	
  void recursive_KPM(int depth, int max_depth, std::vector<int> N_moments, long *average, long *index_gamma, std::vector<std::vector<int>> indices, std::vector<KPM_Vector<T,D>*> *kpm_vector, Eigen::Array<T, -1, -1> *gamma){
    verbose_message("Entered recursive_KPM\n");
    typedef typename extract_value_type<T>::value_type value_type;
		
		
    if(depth != max_depth){
      KPM_Vector<T,D> *kpm1 = kpm_vector->at(depth);
      KPM_Vector<T,D> *kpm2 = kpm_vector->at(depth + 1);
			
      T * kpm1data;
      T * kpm2data;
      int axis1, axis2;
			
      //std::cout << "first branch. Depth: " << depth << "\n" << std::flush;
			
			
      for(int p = 0; p < N_moments.at(depth - 1); p++){
	kpm2->set_index(0);
	switch(indices.at(depth).size()){
	case 0:
	  //std::cout << "case0\n"<< std::flush;
	  break;
	case 1:
	  //std::cout << "case1\n"<< std::flush;
	  axis1 = indices.at(depth).at(0);
	  kpm1data = kpm1->v.col(kpm1->get_index()).data();
	  kpm2data = kpm2->v.col(kpm2->get_index()).data();
						
	  //std::cout << "indices.at(" << depth << "): " << axis1 << "\n"<< std::flush;
	  kpm2->Velocity(kpm2data, kpm1data, axis1); 
						
	  break;
						
	case 2:
	  //std::cout << "case2\n"<< std::flush;
	  axis1 = indices.at(depth).at(0);
	  axis2 = indices.at(depth).at(1);
	  kpm1data = kpm1->v.col(kpm1->get_index()).data();
	  kpm2data = kpm2->v.col(kpm2->get_index()).data();
						
	  kpm2->Velocity2(kpm2data, kpm1data, axis1, axis2); 
						
	  break;
						
	default:
	  //std::cout << "The matrix you're trying to calculate requires an operator that is not yet implemented.\n";
	  exit(0);
						
	}
				
				
	recursive_KPM(depth + 1, max_depth, N_moments, average, index_gamma, indices, kpm_vector, gamma);
	//std::cout << "left second branch\n" << std::flush;
	if(p == 0){
	  //std::cout << "p=0\n" << std::flush;
	  kpm1->template Multiply<0>(); 
	}
	else if(p < N_moments.at(depth-1) - 1){
	  //std::cout << "p!=0\n" << std::flush;
	  kpm1->template Multiply<1>(); 
	}
			
      }
			
    } else {
      KPM_Vector<T,D> *kpm0 = kpm_vector->at(0);
      KPM_Vector<T,D> *kpm1 = kpm_vector->at(depth);
			
      //std::cout << "second branch. Depth: " << depth << "\n" << std::flush<< std::flush;
      kpm1->template Multiply<0>();		
      //std::cout << "You failed at multiplication\n" << std::flush<< std::flush;	
      gamma->matrix().block(0,*index_gamma,1,2) += (kpm0->v.adjoint() * kpm1->v - gamma->matrix().block(0,*index_gamma,1,2))/value_type(*average + 1);			
      *index_gamma += 2;
	
      for(int m = 2; m < N_moments.at(depth - 1); m += 2)
	{
	  kpm1->template Multiply<1>();
	  kpm1->template Multiply<1>();
	  gamma->matrix().block(0, *index_gamma,1,2) += (kpm0->v.adjoint() * kpm1->v - gamma->matrix().block(0,*index_gamma,1,2))/value_type(*average + 1);
	  //std::cout << "product: " << kpm0->v.adjoint() * kpm1->v << "\n";
				
	  *index_gamma += 2;
	}
			
    }
		
    verbose_message("Left recursive_KPM\n");
  }
	
	
  std::vector<std::vector<int>> process_string(std::string indices_string){
    // First of all, split the indices string by commas ','
    std::vector<std::string> strings;
    int end_pos = 0;
    while(end_pos != -1){
      end_pos = indices_string.find(',');
      strings.push_back(indices_string.substr(0, end_pos));
      indices_string = indices_string.substr(end_pos + 1, indices_string.size() - end_pos - 1);
    }
		
    int dim = strings.size();
    /*
      for(int i = 0; i < dim; i++)
      std::cout << strings.at(i) << "\n";*/
		
		
		
		
		
    std::vector<std::vector<int>> indices;
		
    for(int i = 0; i < dim; i++){
			
      int len_str = strings.at(i).size();
      int single_digit;
      std::vector<int> temp;
			
      for(int j = 0; j < len_str; j++){
	char single_char = strings.at(i)[j];
				
	if(single_char == 'x'){
	  single_digit = 0;
	} else {
	  if(single_char == 'y'){
	    single_digit = 1;
	  } else {
	    // This block should never run
	    std::cout << "Please enter a valid expression.\n";
	    exit(0);
	  }
	} 
	temp.push_back(single_digit);
      }
				
      indices.push_back(temp);
			   
    }
			
			
    return indices;
  }

  void store_gamma(Eigen::Array<T, -1, -1> *gamma, std::vector<int> N_moments, std::vector<std::vector<int>> indices, std::string name_dataset){
    debug_message("Entered store_gamma\n");
    /* Depending on the type of Gamma matrix we're calculating, there may be some symmetries
     * among the matrix entries that could be taken into account.
     * 
     * */



    long int size_gamma = gamma->cols();
    int dim = indices.size();

		
    // Number of commutators inside the Gamma matrix. 
    // V^{x}  = [x,H]		-> one commutator
    // V^{xy} = [x,[y,H]]	-> two commutators
    // This is important because the commutator is anti-hermitian. So, an odd number of commutators
    // means that the conjugate of the Gamma matrix has an overall minus sign
    int num_velocities = 0;
    for(int i = 0; i < int(indices.size()); i++)
      num_velocities += indices.at(i).size();
		
    int factor = 1 - (num_velocities % 2)*2;
		
    //std::cout << "gamma: " << *gamma << "\n";
    
    switch(dim)
      {
      case 2: {
	//std::cout << "gamma_matrix dimension: " << dim << "\n" << std::flush;
	Eigen::Array<T,-1,-1> general_gamma = Eigen::Map<Eigen::Array<T,-1,-1>>(gamma->data(), N_moments.at(0), N_moments.at(1));
	//std::cout << "after map?\n" << std::flush;
#pragma omp master
	Global.general_gamma = Eigen::Array<T, -1, -1 > :: Zero(N_moments.at(0), N_moments.at(1));
#pragma omp barrier
	// Gather the data from all the threads, one by one.
	// There are some additional symmetry operations that we can take advantage of.
	// In the case of two indices: Gamma{x,y}(n,m) = Gamma{x,y}(n,m)*, that is, it's self-adjoint.
	// The factor is -1 when the matrix is anti-hermitian and 1 when hermitian
#pragma omp critical
	Global.general_gamma.matrix() += (general_gamma.matrix() + factor*general_gamma.matrix().adjoint())/2.0;
#pragma omp barrier
	break;
      }
      case 1: {
	Eigen::Array<T,-1,-1> general_gamma = Eigen::Map<Eigen::Array<T,-1,-1>>(gamma->data(), 1, size_gamma);
	
#pragma omp master
	Global.general_gamma = Eigen::Array<T, -1, -1 > :: Zero(1, size_gamma);
#pragma omp barrier
	// Gather the data from all the threads, one by one.
	
#pragma omp critical
	Global.general_gamma += general_gamma;
#pragma omp barrier
	break;
      }
	/*
	  case 3: 
	  break;
	*/
      default:
	{
	  std::cout << "You're trying to store a matrix that is not expected by the program. Exiting.\n";
	  exit(0);
	}
      }
    

    
#pragma omp master
    {
      H5::H5File * file = new H5::H5File(name, H5F_ACC_RDWR);
      write_hdf5(Global.general_gamma, file, name_dataset);
      delete file;
    }
#pragma omp barrier    

    
    debug_message("Left store_gamma\n");
  }
  
  double time_kpm(int N_average){
     //This function serves to provide an estimate of the time it takes for each kpm iteration
      
     
#pragma omp barrier
    KPM_Vector<T,D> kpm0(1, *this);
    KPM_Vector<T,D> kpm1(2, *this);
		
		
    kpm0.initiate_vector();
    kpm1.set_index(0);
    kpm1.v.col(0) = kpm0.v.col(0);
    kpm1.template Multiply<0>(); 
		
    auto t0 =  std::chrono::high_resolution_clock::now();
    for(int i = 0; i < N_average; i++)
      kpm1.template Multiply<1>(); 
    auto t1 =  std::chrono::high_resolution_clock::now();
		
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    return time_span.count()/N_average;
		
#pragma omp barrier
  }

  void Single_Shot(double EScale, singleshot_measurement_queue queue) {
    // Calculate the longitudinal dc conductivity for a single value of the energy

    
    // Obtain the relevant quantities from the queue
    int NRandomV = queue.NRandom;
    int NDisorder = queue.NDisorder;
    int N_cheb_moments = queue.NMoments;
    Eigen::Array<double, -1, 1> energy_array = queue.singleshot_energies;
    double finite_gamma = queue.single_gamma;
    std::string indices = queue.direction_string;
    std::string name_dataset = queue.label;
    
		 
		 
    // Process the indices
    debug_message("Entered Single_Shot");
    std::vector<int> first_indices, second_indices;
		
    int comma_location = indices.find_first_of(',');
	
    std::string first_string, second_string;
    first_string = indices.substr(0, comma_location);
    second_string = indices.substr(comma_location+1);		
    first_indices.resize(first_string.size()); 
    second_indices.resize(second_string.size());
	  
    debug_message("Strings: "); debug_message(first_string); 
    debug_message(" ");         debug_message(second_string);debug_message(".\n");

    for(unsigned int i = 0; i < first_string.size(); i++){
      if(first_string[i]=='y')
	      first_indices.at(i) = 1;
      else
	      first_indices.at(i) = 0;
    }
		
    for(unsigned  int i = 0; i < second_string.size(); i++){
      if(second_string[i]=='y')
	      second_indices.at(i) = 1;
      else
	      second_indices.at(i) = 0;
    }
   
    
    //Done processing the indices
		
		
    typedef typename extract_value_type<T>::value_type value_type;
    KPM_Vector<T,D> phi0(1, *this);
    KPM_Vector<T,D> phi (2, *this);
    
    KPM_Vector<T,D> phi1(2, *this);
    KPM_Vector<T,D> phi2(2, *this);
		
	
    //int a = 4;
    //int b = 4;
		
    Eigen::Array<T, -1, -1> cond_array;
    int N_energies = energy_array.cols()*energy_array.rows(); //one of them is one. 
    cond_array = Eigen::Array<T, -1, -1>::Zero(1, N_energies);
    debug_message("Starting calculating.\n");

    for(int ener = 0; ener < N_energies; ener++){
      std::complex<double> energy(energy_array(ener), finite_gamma);
      debug_message("Finished setting complex energy.\n");
      long average = 0;
      for(int disorder = 0; disorder < NDisorder; disorder++){
        debug_message("Before disorder.\n");
	      h.generate_disorder();
        debug_message("After disorder.\n");
    	  for(int randV = 0; randV < NRandomV; randV++){
          debug_message("Started calculating the first vector.\n");
	      phi0.initiate_vector();					
	      phi0.Exchange_Boundaries(); 	
		      
	      // calculate the left vector
	      phi.set_index(0);				
		      
        
	      // |phi> = v |phi_0>
	      phi.Velocity(phi.v.col(0).data(), phi0.v.col(0).data(), first_indices.at(0));
        debug_message("Multiplied by velocity.\n");
	      phi.Exchange_Boundaries();	
	      phi1.v.col(0) = phi.v.col(phi.get_index())*green(0, 1, energy).imag()/2.0;
		      
	      phi.template Multiply<0>();		
	      phi1.v.col(0) += phi.v.col(1)*green(1, 1, energy).imag();
		     
        debug_message("Finished first cheb.\n");
	      for(int n = 2; n < N_cheb_moments; n++){		
      		phi.template Multiply<1>();
      		phi1.v.col(0) += phi.v.col(phi.get_index())*green(n, 1, energy).imag();
	      }
		      
	      // multiply phi1 by the velocity operator again. 
	      // We need a temporary vector to mediate the operation, which will be |phi>
	      phi.v.col(0) = phi1.v.col(0);
	      phi.Velocity(phi1.v.col(0).data(), phi.v.col(0).data(), second_indices.at(0));
	      phi1.empty_ghosts(0);
		  
        debug_message("Finished calculating the first vector.\n");
	      
        // calculate the right vector
	      phi.set_index(0);			
	      phi.v.col(0) = phi0.v.col(0);
		      
	      phi2.v.col(0) = phi.v.col(phi.get_index())*green(0, 1, energy).imag()/2.0;
        
	      phi.template Multiply<0>();		
	      phi2.v.col(0) += phi.v.col(1)*green(1, 1, energy).imag();
		      
	      for(int n = 2; n < N_cheb_moments; n++){		
		      phi.template Multiply<1>();
      		phi2.v.col(0) += phi.v.col(phi.get_index())*green(n, 1, energy).imag();
	      }
        
	      cond_array(ener) += (T(phi1.v.col(0).adjoint()*phi2.v.col(0)) - cond_array(ener))/value_type(average+1);						
	      average++;
        debug_message("Finished calculating the second vector.\n");
	    }
	}
    }
    
    // the gamma matrix has been calculated. Now we're going to use the 
    // property that gamma is hermitian: gamma_nm=gamma_mn*
    
#pragma omp master
    { 
      Global.singleshot_cond = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > :: Zero(1, energy_array.rows()*energy_array.cols());
    }
#pragma omp barrier
    
	  
    //std::cout << "IMPORTANT ! ! !:\n V is not hermitian. Make sure you take this into account\n";
    // in this case there's no problem. both V are anti-hermitic, so the minus signs cancel
#pragma omp critical
    Global.singleshot_cond += cond_array;			
#pragma omp barrier

    
#pragma omp master
    {
			
      // Fixing the factor
      double unit_cell_area = fabs(r.rLat.determinant());
      unsigned int number_of_orbitals = r.Orb; 	// This is necessary because the normalization factor inside the random 
      // vectors is not the number of lattice sites N but the number of states, 
      // which is N*number_of_orbitals
      unsigned int spin_degeneracy = 1;
      
      double factor = -4.0*spin_degeneracy*number_of_orbitals/M_PI/unit_cell_area;	// This is in units of sigma_0, hence the 4
      Global.singleshot_cond *= factor;
      
      // Create array to store the data
      Eigen::Array<double, -1, -1> store_data;
      store_data = Eigen::Array<double, -1, -1>::Zero(2, Global.singleshot_cond.cols());
      
      for(int ener = 0; ener < N_energies; ener++){
        store_data(0, ener) = energy_array.real()(ener)*EScale;
        store_data(1, ener) = Global.singleshot_cond.real()(ener);
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
	
  
};


  
