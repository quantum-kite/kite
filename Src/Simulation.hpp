/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

#ifndef _SIMULATION_HPP
#define _SIMULATION_HPP
#include "queue.hpp"


std::complex<double> green(int n, int sigma, std::complex<double> energy){
  const std::complex<double> i(0.0,1.0); 
  std::complex<double> sq = sqrt(1.0 - energy*energy);
  return 2.0*sigma/sq*i*exp(-sigma*n*1.0*acos(energy)*i);
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

    // Fetch the energy scale
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
      
      if(ESTIMATE_TIME == 1)
        Global.kpm_iteration_time = simul.time_kpm(100);
      
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
  
  

  void cheb_iteration(KPM_Vector<T,D>* kpm, long int current_iteration){
    // Performs a chebyshev iteration
    if(current_iteration == 0){
      kpm->template Multiply<0>(); 
    } else {
      kpm->template Multiply<1>(); 
    }
  };
  

  void generalized_velocity(KPM_Vector<T,D>* kpm0, KPM_Vector<T,D>* kpm1, 
			    std::vector<std::vector<unsigned>> indices, int pos){
    // Check which generalized velocity operator needs to be calculated. 
    // reads from kpm1 and writes on kpm0

    T * kpm1data = kpm1->v.col(kpm1->get_index()).data();
    T * kpm0data = kpm0->v.col(kpm0->get_index()).data();

    switch(indices.at(pos).size()){
    case 0:
      break;
    default:
      kpm0->Velocity(kpm0data, kpm1data, pos); 											
      break;
    }
  }


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
    std::vector<std::vector<unsigned>> indices = process_string(indices_string);
    int dim = indices.size();
		
		
		
    // Check if the dimensions match
    if(dim != int(N_moments.size())){
      std::cout << "Dimension of the Gamma matrix does not match the number of chebyshev moments. Aborting.\n";
      exit(1);
    }
    
    if(dim == 2 and MEMORY > 2){
      if(N_moments.at(0)%MEMORY!=0 or N_moments.at(1)%MEMORY!=0){
        std::cout << "The number of Chebyshev moments ("<< N_moments.at(0)<<","<< N_moments.at(1)<<")"; 
        std::cout << "has to be a multiple of MEMORY ("<< MEMORY <<"). Exiting.\n";
        exit(1);
      }
      Gamma2D(NRandomV, NDisorder, N_moments, indices, name_dataset);
    } else {
      if(dim == 3){
        Gamma3D(NRandomV, NDisorder, N_moments, indices, name_dataset);
      } else {
        GammaGeneral(NRandomV, NDisorder, N_moments, indices, name_dataset);
      }
    }


    debug_message("Left Measure_Gamma\n");
  }
	
  void Gamma2D(int NRandomV, int NDisorder, std::vector<int> N_moments, 
	       std::vector<std::vector<unsigned>> indices, std::string name_dataset){
    Eigen::Matrix<T, MEMORY, MEMORY> tmp;
    // This function calculates all kinds of two-dimensional gamma matrices such
    // as Tr[V^a Tn v^b Tm] = G_nm
    //
    // The matrices are stored as
    //
    // | G_00   G_01   G_02   ...   G_0M | 
    // | G_10   G_11   G_12   ...   G_1M | 
    // | G_20   G_21   G_22   ...   G_2M | 
    // | ...    ...    ...    ...   ...  |
    // | G_N0   G_N1   G_N2   ...   G_NM | 
    //
    // For example, a 3x3 matrix would be represented as
    //
    // | G_00   G_01   G_02 | 
    // | G_10   G_11   G_12 | 
    // | G_20   G_21   G_22 | 
    //
    //
    //

    typedef typename extract_value_type<T>::value_type value_type;

    int num_velocities = 0;
    for(int i = 0; i < int(indices.size()); i++)
      num_velocities += indices.at(i).size();
    int factor = 1 - (num_velocities % 2)*2;

    //  --------- INITIALIZATIONS --------------
    
    KPM_Vector<T,D> kpm0(1, *this);      // initial random vector
    KPM_Vector<T,D> kpm1(2, *this); // left vector that will be Chebyshev-iterated on
    KPM_Vector<T,D> kpm2(MEMORY, *this); // right vector that will be Chebyshev-iterated on
    KPM_Vector<T,D> kpm3(MEMORY, *this); // kpm1 multiplied by the velocity

    // initialize the local gamma matrix and set it to 0
    int size_gamma = 1;
    for(int i = 0; i < 2; i++){
      if(N_moments.at(i) % 2 != 0){
        std::cout << "The number of moments must be an even number, due to limitations of the program. Aborting\n";
        exit(1);
      }
      size_gamma *= N_moments.at(i);
    }

    Eigen::Array<T, -1, -1> gamma = Eigen::Array<T, -1, -1 >::Zero(1, size_gamma);
 
    // finished initializations


    
    
    // start the kpm iteration
    long average = 0;
    for(int disorder = 0; disorder < NDisorder; disorder++){
      h.generate_disorder();
      for(unsigned it = 0; it < indices.size(); it++)
        h.build_velocity(indices.at(it), it);
      for(int randV = 0; randV < NRandomV; randV++){
        

        kpm0.initiate_vector();			// original random vector. This sets the index to zero
        kpm0.Exchange_Boundaries();
	kpm1.set_index(0);

        generalized_velocity(&kpm1, &kpm0, indices, 0);
        
        // run through the left loop MEMORY iterations at a time
        for(int n = 0; n < N_moments.at(0); n+=MEMORY){
          
          // Iterate MEMORY times. The first time this occurs, we must exclude the zeroth
          // case, because it is already calculated, it's the identity
          for(int i = n; i < n + MEMORY; i++){
            if(i!=0){
              cheb_iteration(&kpm1, i-1);
            }

            kpm3.set_index(i%MEMORY);
            generalized_velocity(&kpm3, &kpm1, indices, 1);
            kpm3.empty_ghosts(i%MEMORY);
          }
          
          // copy the |0> vector to |kpm2>
          kpm2.set_index(0);
          kpm2.v.col(0) = kpm0.v.col(0);
          for(int m = 0; m < N_moments.at(1); m+=MEMORY){

            // iterate MEMORY times, just like before. No need to multiply by v here
            for(int i = m; i < m + MEMORY; i++){
              if(i!=0){
                cheb_iteration(&kpm2, i-1);
              }
            }
            
            // Finally, do the matrix product and store the result in the Gamma matrix
            // Eigen::Matrix<T, -1, -1> kpm_product;
            // kpm_product = kpm3.v.adjoint() * kpm2.v; 
	    
	    tmp.setZero();
	    for(std::size_t ii = 0; ii < r.Sized ; ii += r.Ld[0])
	      tmp += kpm3.v.block(ii,0, r.Ld[0], MEMORY).adjoint() * kpm2.v.block(ii, 0, r.Ld[0], MEMORY);
	    
	    
            T flatten;
            long int ind;
            for(int j = 0; j < MEMORY; j++)
              for(int i = 0; i < MEMORY; i++){
                flatten = tmp(i,j);
                ind = (m+j)*N_moments.at(0) + n+i;
                gamma(ind) += (flatten - gamma(ind))/value_type(average + 1);			
              }
          }
        }
        average++;
      }
    } 
    gamma = gamma*factor;

    store_gamma(&gamma, N_moments, indices, name_dataset);
  };


  void Gamma3D(int NRandomV, int NDisorder, std::vector<int> N_moments, 
	       std::vector<std::vector<unsigned>> indices, std::string name_dataset){
    Eigen::Matrix<T, MEMORY, MEMORY> tmp;
    // This calculates all the kinds of three-dimensional gamma matrices
    // such as Tr[v^a Tn v^b Tm v^c Tp] = G_nmp. The output is a 2D matrix 
    // organized as follows:
    // 
    //    p=0     p=1     p=2    ...    p=P
    // 
    // | G_000   G_001   G_002   ...   G_00P | 
    // | G_100   G_101   G_102   ...   G_10P |  
    // | G_200   G_201   G_202   ...   G_20P |   m=0
    // |  ...     ...     ...    ...    ...  |  
    // | G_N00   G_N01   G_N02   ...   G_N1P |_____
    // | G_010   G_011   G_012   ...   G_01P | 
    // | G_110   G_111   G_112   ...   G_11P |  
    // | G_210   G_211   G_212   ...   G_21P |   m=1
    // |  ...     ...     ...    ...    ...  |  
    // | G_N10   G_N11   G_N12   ...   G_N0P |_____ 
    // |  ...     ...     ...    ...    ...  |
    // |  ...     ...     ...    ...    ...  |
    // |  ...     ...     ...    ...    ...  |_____
    // | G_0M0   G_0M1   G_0M2   ...   G_0MP | 
    // | G_1M0   G_1M1   G_1M2   ...   G_1MP |  
    // | G_2M0   G_2M1   G_2M2   ...   G_2MP |   m=M
    // |  ...     ...     ...    ...    ...  |  
    // | G_NM0   G_NM1   G_NM2   ...   G_NMP |_____ 
    //
    // To exemplify, in the case of a 3x3x3 matrix, 
    //
    // | G_000   G_001   G_002 |
    // | G_100   G_101   G_102 |
    // | G_200   G_201   G_202 |
    // | G_010   G_011   G_012 |
    // | G_110   G_111   G_112 |
    // | G_210   G_211   G_212 |
    // | G_020   G_021   G_022 |
    // | G_120   G_121   G_122 |
    // | G_220   G_221   G_222 |
    //
    //

    typedef typename extract_value_type<T>::value_type value_type;

    //  --------- INITIALIZATIONS --------------
    
    KPM_Vector<T,D> kpm0(1, *this);           // initial random vector
    KPM_Vector<T,D> kpm_Vn(2, *this);          // left vector that will be Chebyshev-iterated on
    KPM_Vector<T,D> kpm_VnV(MEMORY, *this);    // kpmL multiplied by the velocity
    KPM_Vector<T,D> kpm_p(2, *this);          // right-most vector that will be Chebyshev-iterated on
    KPM_Vector<T,D> kpm_pVm(MEMORY, *this);         // middle vector that will be Chebyshev-iterated on

    // initialize the local gamma matrix and set it to 0
    int size_gamma = 1;
    for(int i = 0; i < 3; i++){
      if(N_moments.at(i) % 2 != 0){
        std::cout << "The number of moments must be an even number, due to limitations of the program. Aborting\n";
        exit(1);
      }
      size_gamma *= N_moments.at(i);
    }
    Eigen::Array<T, -1, -1> gamma = Eigen::Array<T, -1, -1 >::Zero(1, size_gamma);
 
    // finished initializations
    
    // start the kpm iteration
    long average = 0;
    for(int disorder = 0; disorder < NDisorder; disorder++){

      // Distribute the disorder and update the velocity matrices
      h.generate_disorder();
      for(unsigned it = 0; it < indices.size(); it++)
        h.build_velocity(indices.at(it), it);

      for(int randV = 0; randV < NRandomV; randV++){
        

        kpm0.initiate_vector();			// original random vector. This sets the index to zero
        kpm0.Exchange_Boundaries();
	kpm_Vn.set_index(0);

        generalized_velocity(&kpm_Vn, &kpm0, indices, 0);
        
        for(int n = 0; n < N_moments.at(0); n+=MEMORY){

          // Calculation of the left kpm vector
          for(int ni = n; ni < n + MEMORY; ni++){
            if(ni!=0) cheb_iteration(&kpm_Vn, ni-1);
           
            kpm_VnV.set_index(ni%MEMORY);
            generalized_velocity(&kpm_VnV, &kpm_Vn, indices, 1);
            kpm_VnV.empty_ghosts(ni%MEMORY);
          }
          
          // Calculation of the right kpm vector
          kpm_p.set_index(0);
          kpm_p.v.col(0) = kpm0.v.col(0);
          for(int p = 0; p < N_moments.at(2); p++){
            if(p!=0) cheb_iteration(&kpm_p, p-1);
            
            kpm_pVm.set_index(0);
            generalized_velocity(&kpm_pVm, &kpm_p, indices, 2);
            for(int m = 0; m < N_moments.at(1); m += MEMORY){
              for(int mi = m; mi < m + MEMORY; mi++)
                if(mi != 0) cheb_iteration(&kpm_pVm, mi-1);

              
	    
	      tmp.setZero();
	      for(std::size_t ii = 0; ii < r.Sized ; ii += r.Ld[0])
		tmp += kpm_VnV.v.block(ii,0, r.Ld[0], MEMORY).adjoint() * kpm_pVm.v.block(ii, 0, r.Ld[0], MEMORY);
              //Eigen::Matrix<T, -1, -1> kpm_product;
	      
	      // kpm_product = Eigen::Matrix<T, -1, -1>::Zero(MEMORY, MEMORY); // this line is not necessary
              //kpm_product = kpm_VnV.v.adjoint() * kpm_pVm.v; 
	      
              long int index;
              for(int i = 0; i < MEMORY; i++)
                for(int j = 0; j < MEMORY; j++){
                  index = p*N_moments.at(1)*N_moments.at(0) + (m+j)*N_moments.at(0) + n+i;
                  gamma(index) += (tmp(i, j) - gamma(index))/value_type(average + 1);
                }
            }
          }
        }
        average++;
      }
    } 
    store_gamma(&gamma, N_moments, indices, name_dataset);
  };

  void GammaGeneral(int NRandomV, int NDisorder, std::vector<int> N_moments,
		    std::vector<std::vector<unsigned>> indices, std::string name_dataset){
    
    int dim = indices.size();
		
    // Check if the dimensions match
    if(dim != int(N_moments.size())){
      std::cout << "Dimension of the Gamma matrix does not match the number of chebyshev moments. Aborting.\n";
      exit(1);
    }
			
    // Determine the size of the gamma matrix we want to calculate
    int size_gamma = 1;
    for(int i = 0; i < dim; i++){
      if(N_moments.at(i) % 2 != 0){
        std::cout << "The number of moments must be an even number, due to limitations of the program. Aborting\n";
        exit(1);
      }
      size_gamma *= N_moments.at(i);
    }
		
		
    // Initialize the KPM vectors that will be needed to run the program 
    std::vector<KPM_Vector<T,D>*> kpm_vector(dim+1);
    kpm_vector.at(0) = new KPM_Vector<T,D> (1, *this);
    for(int i = 0; i < dim; i++)
      kpm_vector.at(i+1) = new KPM_Vector<T,D> (2, *this);
		
		
    // Define some pointers to make the code easier to read
    KPM_Vector<T,D> *kpm0 = kpm_vector.at(0);
    KPM_Vector<T,D> *kpm1 = kpm_vector.at(1);
			
    // Make sure the local gamma matrix is zeroed
    Eigen::Array<T, -1, -1> gamma = Eigen::Array<T, -1, -1 >::Zero(1, size_gamma);

    long average = 0;
    for(int disorder = 0; disorder < NDisorder; disorder++){
      h.generate_disorder();
      for(unsigned it = 0; it < indices.size(); it++)
        h.build_velocity(indices.at(it), it);

      for(int randV = 0; randV < NRandomV; randV++){
        
        kpm0->initiate_vector();			// original random vector
        kpm1->set_index(0);
        kpm1->v.col(0) = kpm0->v.col(0);
        kpm1->Exchange_Boundaries();

        // replace <0| by  <0|v. Note that v is not self-adjoint in this formulation
        generalized_velocity(kpm0, kpm1, indices, 0);
        int factor = 1 - (indices.at(0).size() % 2)*2;		
        kpm0->v.col(0) = factor*kpm0->v.col(0); // This factor is due to the fact that this Velocity operator is not self-adjoint

        kpm0->empty_ghosts(0);
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
	
  }

  void recursive_KPM(int depth, int max_depth, std::vector<int> N_moments, long *average, long *index_gamma, 
		     std::vector<std::vector<unsigned>> indices, std::vector<KPM_Vector<T,D>*> *kpm_vector, Eigen::Array<T, -1, -1> *gamma){
    debug_message("Entered recursive_KPM\n");
    typedef typename extract_value_type<T>::value_type value_type;
    Eigen::Matrix < T, 1, 2> tmp =  Eigen::Matrix < T, 1, 2> ::Zero();		
		
    if(depth != max_depth){
      KPM_Vector<T,D> *kpm1 = kpm_vector->at(depth);
      KPM_Vector<T,D> *kpm2 = kpm_vector->at(depth + 1);
			
      T * kpm1data;
      T * kpm2data;
			
      //std::cout << "first branch. Depth: " << depth << "\n" << std::flush;
			
			
      for(int p = 0; p < N_moments.at(depth - 1); p++){
        kpm2->set_index(0);
        switch(indices.at(depth).size()){
        case 0:
          break;
        default:
          kpm1data = kpm1->v.col(kpm1->get_index()).data();
          kpm2data = kpm2->v.col(kpm2->get_index()).data();
                  
          kpm2->Velocity(kpm2data, kpm1data, max_depth - depth); 											
        }
				
        recursive_KPM(depth + 1, max_depth, N_moments, average, index_gamma, indices, kpm_vector, gamma);
        if(p == 0){
          kpm1->template Multiply<0>(); 
        }
        else if(p < N_moments.at(depth-1) - 1){
          kpm1->template Multiply<1>(); 
        }
      }
			
    } else {
      KPM_Vector<T,D> *kpm0 = kpm_vector->at(0);
      KPM_Vector<T,D> *kpm1 = kpm_vector->at(depth);
			
      kpm1->template Multiply<0>();

      tmp.setZero();
      for(std::size_t ii = 0; ii < r.Sized ; ii += r.Ld[0])
	tmp += kpm0->v.block(ii,0, r.Ld[0], 1).adjoint() * kpm1->v.block(ii, 0, r.Ld[0], 2);
      
      gamma->matrix().block(0,*index_gamma,1,2) += (tmp - gamma->matrix().block(0,*index_gamma,1,2))/value_type(*average + 1);			
      *index_gamma += 2;
	
      for(int m = 2; m < N_moments.at(depth - 1); m += 2){
        kpm1->template Multiply<1>();
        kpm1->template Multiply<1>();
	
	tmp.setZero();
	for(std::size_t ii = 0; ii < r.Sized ; ii += r.Ld[0])
	  tmp += kpm0->v.block(ii,0, r.Ld[0], 1).adjoint() * kpm1->v.block(ii, 0, r.Ld[0], 2);
	
        gamma->matrix().block(0, *index_gamma,1,2) += (tmp - gamma->matrix().block(0,*index_gamma,1,2))/value_type(*average + 1);
            
        *index_gamma += 2;
      }
    }
		
    debug_message("Left recursive_KPM\n");
  }
	
	
  std::vector<std::vector<unsigned>> process_string(std::string indices_string){
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
		
		
		
		
		
    std::vector<std::vector<unsigned>> indices;
		
    for(int i = 0; i < dim; i++){
			
      int len_str = strings.at(i).size();
      int single_digit;
      std::vector<unsigned> temp;
			
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
	    exit(1);
	  }
	} 
	temp.push_back(single_digit);
      }
				
      indices.push_back(temp);
			   
    }
			
			
    return indices;
  }

  void store_gamma(Eigen::Array<T, -1, -1> *gamma, std::vector<int> N_moments, 
		   std::vector<std::vector<unsigned>> indices, std::string name_dataset){
    debug_message("Entered store_gamma\n");
    // The whole purpose of this function is to take the Gamma matrix calculated by
    // Gamma2D, Gamma3D or Gamma_general, check if there are any symmetries among the 
    // matrix elements and then store the matrix in an HDF file.

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

		
    switch(dim){
    case 2:{
      Eigen::Array<T,-1,-1> general_gamma = Eigen::Map<Eigen::Array<T,-1,-1>>(gamma->data(), N_moments.at(0), N_moments.at(1));
#pragma omp master
      Global.general_gamma = Eigen::Array<T, -1, -1 > :: Zero(N_moments.at(0), N_moments.at(1));
#pragma omp barrier
#pragma omp critical
      Global.general_gamma.matrix() += (general_gamma.matrix() + factor*general_gamma.matrix().adjoint())/2.0;
#pragma omp barrier
      break;
    }
    case 1:{
      Eigen::Array<T,-1,-1> general_gamma = Eigen::Map<Eigen::Array<T,-1,-1>>(gamma->data(), 1, size_gamma);
#pragma omp master
      Global.general_gamma = Eigen::Array<T, -1, -1 > :: Zero(1, size_gamma);
#pragma omp barrier
#pragma omp critical
      Global.general_gamma += general_gamma;
#pragma omp barrier
      break;
    }
    case 3:{
      int N0 = N_moments.at(0);
      int N1 = N_moments.at(1);
      int N2 = N_moments.at(2);

      Eigen::Array<T,-1,-1> general_gamma = Eigen::Map<Eigen::Array<T,-1,-1>>(gamma->data(), N0*N1, N2);
#pragma omp master
      Global.general_gamma = Eigen::Array<T, -1, -1 > :: Zero(N0*N1, N2);
#pragma omp barrier

#pragma omp critical
      {
        // Check if all the directions are the same. In this case, there are 
        // six symmetries that may be taken advantage of ('*' denotes complex conjugation)
        // G_nmp  = G_mpn = G_pnm 
        // G_nmp* = G_pmn = G_mnp = G_npm 
        if(indices.at(0) == indices.at(1) and indices.at(0) == indices.at(2)){
          for(int n = 0; n < N0; n++){
            for(int m = 0; m < N1; m++){
              for(int p = 0; p < N2; p++){
                Global.general_gamma(n + N0*m,p) += general_gamma(n + N0*m,p)/T(6.0);
                Global.general_gamma(n + N0*m,p) += general_gamma(m + N0*p,n)/T(6.0);
                Global.general_gamma(n + N0*m,p) += general_gamma(p + N0*n,m)/T(6.0);

		Global.general_gamma(n + N0*m,p) += T(factor/6.0)*myconj(general_gamma(p + N0*m,n));
		Global.general_gamma(n + N0*m,p) += T(factor/6.0)*myconj(general_gamma(n + N0*p,m));
		Global.general_gamma(n + N0*m,p) += T(factor/6.0)*myconj(general_gamma(m + N0*n,p));
		
              }
            }
          }
          
        } 
        
        // Now check if any two directions are the same
        // Check if the two first directions are the same but different from the third
        if(indices.at(0) == indices.at(1) and indices.at(0) != indices.at(2) and N1 == N2){
          for(int n = 0; n < N0; n++){
            for(int m = 0; m < N1; m++){
              for(int p = 0; p < N2; p++){
                Global.general_gamma(n + N0*m,p) += general_gamma(n + N0*m,p)/T(2.0);
		Global.general_gamma(n + N0*m,p) += T(factor/2.0)*myconj(general_gamma(n + N0*p,m));
	      }	
	    }
          }
        }

        // Check if the first and last directions are the same but different from the second
        if(indices.at(0) == indices.at(2) and indices.at(0) != indices.at(1) and N0 == N2){
          for(int n = 0; n < N0; n++){
            for(int m = 0; m < N1; m++){
              for(int p = 0; p < N2; p++){
                Global.general_gamma(n + N0*m,p) += general_gamma(n + N0*m,p)/T(2.0);
		Global.general_gamma(n + N0*m,p) += T(factor/2.0)*myconj(general_gamma(m + N0*n,p));
              }
            }
          }
        }

        // Check if the last two directions are the same but different from the first
        if(indices.at(2) == indices.at(1) and indices.at(0) != indices.at(2) and N0 == N1){
          for(int n = 0; n < N0; n++){
            for(int m = 0; m < N1; m++){
              for(int p = 0; p < N2; p++){
                Global.general_gamma(n + N0*m,p) += general_gamma(n + N0*m,p)/T(2.0);
		Global.general_gamma(n + N0*m,p) += T(factor/2.0)*myconj(general_gamma(p + N0*m, n));

              }
            }
          }
        }

        // Check if all the directions are different
        if(indices.at(0) != indices.at(1) and indices.at(0) != indices.at(2) and indices.at(1) != indices.at(2)){
          Global.general_gamma += general_gamma;
        }


      }
#pragma omp barrier

      break;
    }
    default:
      std::cout << "You're trying to store a matrix that is not expected by the program. Exiting.\n";
      exit(1);
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

  void Single_Shot(double EScale, singleshot_measurement_queue queue)
  {
    // Calculate the longitudinal dc conductivity for a single value of the energy
    T tmp;
    debug_message("Entered Single_Shot\n");
    
    // Obtain the relevant quantities from the queue
    int NRandomV = queue.NRandom;
    int NDisorder = queue.NDisorder;
    Eigen::Array<double, -1, -1> jobs = queue.singleshot_energiesgammas;
    std::string indices_string = queue.direction_string;
    std::string name_dataset = queue.label;
    int N_energies = jobs.rows();
    
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
    if(indices.at(0).at(0) != indices.at(1).at(0))
      {
	std::cout << "SingleShot is only meaningful for the longitudinal conductivity.";
	std::cout << "Please use directions 'x,x' or 'y,y'. Exiting.\n";
	exit(1);
      }
    
    
    // initialize the kpm vectors necessary for this calculation
    typedef typename extract_value_type<T>::value_type value_type;
    
    // if the SSPRINT flag is true, we need one kpm vector for the right vector
    // and one for the left vector. Otherwise, we can just recycle it
#if (SSPRINT == 0)
    KPM_Vector<T,D> phi (2, *this);
#elif (SSPRINT != 0)
    KPM_Vector<T,D> phir1 (2, *this);
    KPM_Vector<T,D> phir2 (2, *this);
#endif
    
    // right and left vectors
    KPM_Vector<T,D> phi0(1, *this);
    KPM_Vector<T,D> phi1(1, *this);
    
    // if SSPRINT is true, we need a temporary vector to store v|phi>
#if (SSPRINT != 0)
    KPM_Vector<T,D> phi2(1, *this);
#endif
    
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
    for(int disorder = 0; disorder < NDisorder; disorder++)
      {
	h.generate_disorder();
	h.build_velocity(indices.at(0),0u);
	h.build_velocity(indices.at(1),1u);
	
	// iteration over each energy and gammma
	for(int job_index = 0; job_index < N_energies; job_index++)
	  {
	    job_energy = jobs(job_index, 0);
	    job_gamma = jobs(job_index, 1);
	    job_preserve_disorder = jobs(job_index, 2);
	    job_NMoments = int(jobs(job_index,3));
	    std::complex<double> energy(job_energy, job_gamma);
	    
	    if(job_preserve_disorder == 0.0)
	      {
		h.generate_disorder();
		h.build_velocity(indices.at(0), 0u);
		h.build_velocity(indices.at(1), 1u);
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
	      
	      // initialize the kpm vectors that will be used in the kpm recursion
	      // the left vector is multiplied by the velocity
	      phi.set_index(0);				
	      generalized_velocity(&phi, &phi0, indices, 0);      // |phi> = v |phi_0>
	      
	      
	      for(int n = 0; n < job_NMoments; n++)
		{		
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
	      
	      // calculate the right KPM vector. Now we may reuse phi0
	      phi.set_index(0);			
	      phi.v.col(0) = phi0.v.col(0);
	      phi0.v.col(0).setZero(); 
	      
	      for(int n = 0; n < job_NMoments; n++)
		{		
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
	      
	      for(int nn = 0; nn < SSPRINT; nn++)
		{
		  for(int n = nn*job_NMoments/SSPRINT; n < job_NMoments/SSPRINT*(nn+1); n++)
		    {	
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
		  
		  
		  for(int n = nn*job_NMoments/SSPRINT; n < job_NMoments/SSPRINT*(nn+1); n++)
		    {		
		      if(n!=0) cheb_iteration(&phir2, n-1);
		      
		      phi0.v.col(0) += phir2.v.col(phir2.get_index())
			*green(n, 1, energy).imag()/(1.0 + int(n==0));
		    }
		  
		  // This is the conductivity for a smaller number of chebyshev moments
		  // if you want to add it to the average conductivity, you have yo wait
		  // until all the moments have been summed. otherwise the result would be wrong
		  // Alterar

		  
		  if(nn == SSPRINT-1)
		    {
		      T tmp ;
		      tmp *= 0.;
		      for(std::size_t ii = 0; ii < r.Sized ; ii += r.Ld[0])
			tmp += phi2.v.col(0).segment(ii,r.Ld[0]).adjoint() * phi0.v.col(0).segment(ii,r.Ld[0]);
		      cond_array(job_index) += (tmp - cond_array(job_index))/value_type(average_R + 1);
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
      store_data = Eigen::Array<double, -1, -1>::Zero(4, jobs.rows());
      
      for(int ener = 0; ener < N_energies; ener++)
	{
	  store_data(0, ener) = jobs.real()(ener, 0)*EScale;
	  store_data(1, ener) = jobs.real()(ener, 1)*EScale;
	  store_data(2, ener) = jobs(ener, 3);
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
};



