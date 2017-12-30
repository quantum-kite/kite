#ifndef _SIMULATION_HPP
#define _SIMULATION_HPP

#define DEBUG 0

int custom_find(int *arr, int arr_size, int value_to_search){
	/* This function searches array 'arr' for any occurence of number 'value to search'.
	 * If it exists, it exits and returns the index where it occurs.
	 * Otherwise, returns -1
	 */
	for(int i = 0; i < arr_size; i++){
		if(*(arr+i) == value_to_search)
			return i;
	}
	
	return -1;	
}


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
  std::vector <int> Quantities, NMoments, NRandomV, NDisorder; 
  
  // Other quantities that require special care, such as SingleShotXX
  std::vector <int> Quantities_special, NMoments_special, NRandomV_special, NDisorder_special, EnergyScale, MagneticField; 
  std::vector <double> gamma_special;
  Eigen::Array<double, -1, -1> singleshot_energies;
  
public:
  GlobalSimulation( char *name ) : rglobal(name) {
    Global.ghosts.resize( rglobal.get_BorderSize() );
    std::fill(Global.ghosts.begin(), Global.ghosts.end(), 0);
    
    /*
     *   /Calculation/FunctionNum:
     *
     *   DOS - denstity of states == 1,
     *   CondXX - conductivity in xx direction == 2,
     *   CondXY - conductivity in xy direction == 3,
     *   OptCond - optical conductivity == 4
     *   SpinCond - spin conductivity == 5

     * /Calculation/NumRandoms   : number of different random vector realisations,
     * /Calculation/NumMoments   : number of moments for the calculation,
     * /Calculation/NumDisorder  : number of disorder realisations.
     */
    
    // Regular quantities to calculate, such as DOS and CondXX
    H5::H5File * file         = new H5::H5File(name, H5F_ACC_RDONLY);
    H5::DataSet * dataset     = new H5::DataSet(file->openDataSet("/Calculation/FunctionNum"));
    H5::DataSpace * dataspace = new H5::DataSpace(dataset->getSpace());
    size_t NQuantities        = dataspace->getSimpleExtentNpoints();
    
    Quantities.resize (NQuantities);
    NMoments.resize   (NQuantities);
    NRandomV.resize   (NQuantities);
    NDisorder.resize  (NQuantities);
    EnergyScale.resize  (1);
    MagneticField.resize  (1);

    get_hdf5<int>(Quantities.data(), file, (char *)   "/Calculation/FunctionNum");
    get_hdf5<int>(NRandomV.data(),   file, (char *)   "/Calculation/NumRandoms");
    get_hdf5<int>(NMoments.data(),   file, (char *)   "/Calculation/NumMoments");
    get_hdf5<int>(NDisorder.data(),  file, (char *)   "/Calculation/NumDisorder");
    get_hdf5<int>(EnergyScale.data(),  file, (char *)   "/EnergyScale");
    MagneticField.at(0) = 0;
		try {
			H5::Exception::dontPrint();
			get_hdf5<int>(MagneticField.data(),  file, (char *)   "/Hamiltonian/MagneticField");
		}
		catch (H5::Exception& e){}
    
    delete dataspace;
    delete dataset;
    delete file;
    
    // What quantities from this list do we need to calculate? Let's find out
    int DOS 	 = custom_find(Quantities.data(), NQuantities, 1);
    int CondXX = custom_find(Quantities.data(), NQuantities, 2);
    int CondXY = custom_find(Quantities.data(), NQuantities, 3);

   
      
    
	// Other quantities that require special care, such as SingleShotXX
	int SingleShotXX = -1;
	int SingleShotXY = -1;
	H5::H5File * file_special	= new H5::H5File(name, H5F_ACC_RDONLY);
	
	try{
		if(DEBUG)std::cout << "start try\n";fflush(stdout);
		H5::Exception::dontPrint();
		
		// This is here just to determine the number of quantities we need to calculate
		H5::DataSet * dataset_special     	= new H5::DataSet(file_special->openDataSet("/Calculation/Calculation_spec/FunctionNum"));
		H5::DataSpace * dataspace_special 	= new H5::DataSpace(dataset_special->getSpace());
		size_t NQuantities_special  				= dataspace_special->getSimpleExtentNpoints();
		delete dataspace_special;
		delete dataset_special;
		
		Quantities_special.resize (NQuantities_special);
		NMoments_special.resize   (NQuantities_special);
		NRandomV_special.resize   (NQuantities_special);
		NDisorder_special.resize  (NQuantities_special);
		gamma_special.resize 	  (NQuantities_special);
		
		
		
		
		// We also need to determine the number of energies that we need to calculate
		H5::DataSet * dataset_energy     	= new H5::DataSet(file_special->openDataSet("/Calculation/Calculation_spec/Energy"));
		H5::DataSpace * dataspace_energy 	= new H5::DataSpace(dataset_energy->getSpace());
		hsize_t dims_out[2];		
		dataspace_energy->getSimpleExtentDims(dims_out, NULL);	
		singleshot_energies = Eigen::Array<double, -1, -1>::Zero(dims_out[0], dims_out[1]);	
		delete dataspace_energy;
		delete dataset_energy;
			
			
			
			// Now let's fetch those quantities
			get_hdf5<int>(Quantities_special.data(),	file, (char *)   "/Calculation/Calculation_spec/FunctionNum");
			get_hdf5<int>(NRandomV_special.data(),  	file, (char *)   "/Calculation/Calculation_spec/NumRandoms");
			get_hdf5<int>(NMoments_special.data(),  	file, (char *)   "/Calculation/Calculation_spec/NumMoments");
			get_hdf5<int>(NDisorder_special.data(), 	file, (char *)   "/Calculation/Calculation_spec/NumDisorder");
			get_hdf5<double>(gamma_special.data(),  	file, (char *)   "/Calculation/Calculation_spec/Gamma");
			get_hdf5<double>(singleshot_energies.data(),  	file, (char *)   "/Calculation/Calculation_spec/Energy");
			
			// Make sure the moments are a power of 2
			for(unsigned i = 0; i < NMoments_special.size(); i++)
				NMoments_special.at(i) = 2*(NMoments_special.at(i)/2);

			
			
			SingleShotXX = custom_find(Quantities_special.data(), NQuantities_special, 6);
			SingleShotXY = custom_find(Quantities_special.data(), NQuantities_special, 7);
			if(DEBUG)std::cout << "ended try\n";fflush(stdout);
    }
    catch(H5::Exception& e) {
			if(DEBUG)std::cout << "exception \n";fflush(stdout);
    }
    
		delete file_special;
			
    omp_set_num_threads(rglobal.n_threads);
#pragma omp parallel default(shared)
    {
      
      Simulation<T,D> simul(name, Global);	  
            
      if(SingleShotXX != -1){
				if(DEBUG)std::cout << "calculating of singleshotxx\n";fflush(stdout);
				Global.singleshot_cond = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > :: Zero(1, singleshot_energies.cols());
				simul.Single_Shot(EnergyScale.at(0), NRandomV_special.at(SingleShotXX), NDisorder_special.at(SingleShotXX), NMoments_special.at(SingleShotXX), singleshot_energies.row(SingleShotXX), gamma_special.at(SingleShotXX), "x,x", "SingleShotXX");
				
				if(DEBUG)std::cout << "ended singleshotxx\n";fflush(stdout);
				
			}
       
      if(SingleShotXY != -1){
				if(DEBUG)std::cout << "calculating of singleshotxy\n";fflush(stdout);
				Global.singleshot_cond = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > :: Zero(1, singleshot_energies.cols());
				simul.Single_Shot(EnergyScale.at(0), NRandomV_special.at(SingleShotXY), NDisorder_special.at(SingleShotXY), NMoments_special.at(SingleShotXY), singleshot_energies.row(SingleShotXY), gamma_special.at(SingleShotXY), "x,y", "SingleShotXY");
				if(DEBUG)std::cout << "ended singleshotxx\n";fflush(stdout);
				
			}
      if(DOS != -1){
				if(DEBUG)std::cout << "calculating of DOS\n";fflush(stdout);
				if(DEBUG)std::cout << "NRandomV: " << NRandomV.at(DOS) << "\n";
				if(DEBUG)std::cout << "NDisorder: " << NDisorder.at(DOS) << "\n";
				if(DEBUG)std::cout << "NMoments: " << NMoments.at(DOS) << "\n";
				
				Global.mu = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > :: Zero(1, NMoments.at(DOS));
				simul.Measure_Dos(NRandomV.at(DOS), NDisorder.at(DOS) );
				if(DEBUG)std::cout << "ended DOS\n";fflush(stdout);
			}
			
			if(CondXX != -1){
				if(DEBUG)std::cout << "calculating of condxx\n";fflush(stdout);
				Global.gamma = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > :: Zero(NMoments.at(CondXX), NMoments.at(CondXX));
				Global.lambda = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > :: Zero(1, NMoments.at(CondXX));
				simul.Measure_Cond(NRandomV.at(CondXX), NDisorder.at(CondXX), "x,x", "GammaXX");
				simul.Measure_Lambda(NRandomV.at(CondXX), NDisorder.at(CondXX), "xx", "LambdaXX");
				if(DEBUG)std::cout << "ended condxx\n";fflush(stdout);
			}
			
			if(CondXY != -1){
				if(DEBUG)std::cout << "calculating of condxy\n";fflush(stdout);
				Global.gamma = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > :: Zero(NMoments.at(CondXY), NMoments.at(CondXY));
				Global.lambda = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > :: Zero(1, NMoments.at(CondXY));
				simul.Measure_Cond(NRandomV.at(CondXY), NDisorder.at(CondXY), "x,y", "GammaXY");
				simul.Measure_Lambda(NRandomV.at(CondXY), NDisorder.at(CondXY), "xy", "LambdaXY");
				if(DEBUG)std::cout << "ended condxy\n";fflush(stdout);
			}
			
			
			
			
    }
    
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
  
  
  void Measure_Dos(int & NRandomV, int & NDisorder) {
    typedef typename extract_value_type<T>::value_type value_type;
    Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> mu = Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic > :: Zero(Global.mu.rows(), Global.mu.cols());
    KPM_Vector<T,D>  phi0(1, *this);
    KPM_Vector<T,D>  phi (2, *this);

    //    phi0.test_boundaries_system();

#pragma omp barrier
    long average = 0;

    for(int disorder = 0; disorder <  NDisorder; disorder++)
      {
	h.generate_disorder();
	for(int randV = 0; randV < NRandomV; randV++)
	  {
	    phi0.initiate_vector();
	    phi.set_index(0);
	    phi.v.col(0) = phi0.v.col(0);
	    phi.Exchange_Boundaries();
	    
	    phi.template Multiply<0>();
	    mu.matrix().block(0,0,1,2) +=  (phi0.v.adjoint() * phi.v - mu.matrix().block(0,0,1,2))/value_type(average + 1);
	    
	    for(int m = 2; m < mu.cols(); m += 2)
	      {	    
		phi.template Multiply<1>();
		phi.template Multiply<1>();
		mu.matrix().block(0,m,1,2) +=  (phi0.v.adjoint() * phi.v - mu.matrix().block(0,m,1,2))/value_type(average + 1);
	      }
	    average++;
	  }
	  
	  if(DEBUG)std::cout << "Finished chb iteration in DOS.\n";fflush(stdout);

#pragma omp critical
	Global.mu += mu;
#pragma omp barrier
	
	
#pragma omp master
	{
		if(DEBUG)std::cout << "before creating file with name: " << name << "\n";fflush(stdout);
	  H5::H5File * file1 = new H5::H5File(name, H5F_ACC_RDWR);
	  if(DEBUG)std::cout << "before write to file\n";fflush(stdout);
	  write_hdf5(Global.mu, file1, "MU");
	  if(DEBUG)std::cout << "before delete file\n";fflush(stdout);
	  delete file1;
	  Global.mu.setZero();
	}
#pragma omp barrier	
      }
#pragma omp barrier	

	if(DEBUG)std::cout << "Left calculation of DOS\n";fflush(stdout);
  }
  
  
  
  void Measure_Cond(int & NRandomV, int & NDisorder, std::string indices, std::string name_dataset) {
		/*
		 * Calculate the gamma matrix gammaxx = Tr[v^x T_n v^x T_m] 
		 * used to calculate the longitudinal conductivity
		 */
		 
		// Process the indices
		if(DEBUG)std::cout << "Calculating cond\n";
		std::vector<int> first_indices, second_indices;
		
		int comma_location = indices.find_first_of(',');
	
		std::string first_string, second_string;
		first_string = indices.substr(0, comma_location);
		second_string = indices.substr(comma_location+1);		
		first_indices.resize(first_string.size()); 
		second_indices.resize(second_string.size());
		
		if(DEBUG)std::cout << "strings: " << first_string << " "<< second_string << "\n";
		
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
			/*
		for(unsigned  int i = 0; i < first_string.size(); i++)
			std::cout << "first_indices: " << first_indices.at(i) << "\n";
			
		for(unsigned  int i = 0; i < second_string.size(); i++)
			std::cout << "second_indices: " << second_indices.at(i) << "\n";
		*/
		
		
		//Done processing the indices
		
		
		
		
		 
		
    typedef typename extract_value_type<T>::value_type value_type;
    KPM_Vector<T,D>  phi0(1, *this);
    KPM_Vector<T,D> phi_n(2, *this);
    KPM_Vector<T,D> phi_m(2, *this);
		
		//make sure the local gamma matrix is zeroed
		Eigen::Array<T, -1, -1> gamma = Eigen::Array<T, -1, -1 > :: Zero(Global.gamma.rows(), Global.gamma.cols());
		//gamma.setZero();
		
		
		/* We want to calculate <0|vmvn|0>
		 * Define |n> = n|0> and store it. 
		 * After having |n> stored, define |n'> = v |n>, which will act as the new |0> for m
		 * Now, having |n'>, define |m> = m |n'>. These may all be calculated using 
		 * Chebyshev's recursion relations. So, for a fixed |n>, we have found all the |m>
		 * Since |m> = m v n |0>, what remains is to do the dot product <0|v|m>, which
		 * gives us the Gamma-matrix-element we sought. After this, we go back to |n>. Using 
		 * the recursion relations, we may find the next |n+1>, and repeat the process all over again
		 * */
		
		int N_cheb_moments = gamma.cols(); //assuming cols=rows
		//std::cout << "calculating vmvn\n";fflush(stdout);
	
    long average = 0;
    for(int disorder = 0; disorder < NDisorder; disorder++){
			h.generate_disorder();
			for(int randV = 0; randV < NRandomV; randV++){
				
				phi0.initiate_vector();					
				phi_n.set_index(0);
				phi_n.v.col(0) = phi0.v.col(0); 
				phi_n.Exchange_Boundaries(); 	
				//phi_n.Velocity(phi0.v.col(0).data(), phi_n.v.col(phi_n.get_index()).data(), 0);
				//replace |phi0> by v |phi0>
				
				
				// DONT FORGET V IS NOT HERMITIAN, so multiplication by <n| gives a minus sign
				if(first_indices.size() == 1){
					//std::cout << indices << "first index has 1\n";
					phi_n.Velocity(phi0.v.col(0).data(), phi_n.v.col(phi_n.get_index()).data(), first_indices.at(0));
					phi0.v.col(0) = -phi0.v.col(0); // compensating the fact that V is anti-hermitian
				}
				
				if(first_indices.size() == 2){
					//std::cout << indices << "first index has 2\n";
					phi_n.Velocity2(phi0.v.col(0).data(), phi_n.v.col(phi_n.get_index()).data(), first_indices.at(0), first_indices.at(1));
					//this one is hermitian, no problem here
				}
				
				phi0.empty_ghosts(0);
				
				
				// actually we should have multiplied by the adjoint of the velocity. 
				// Since we forgot the "i" in V, adjoint(Velocity) = -Velocity
				// The only difference is an overall minus sign
				
				
				for(int n = 0; n < N_cheb_moments; n += 1){
					phi_m.set_index(0);
					
					if(second_indices.size() == 1){
						//std::cout << indices <<"second index has 1\n";
						phi_n.Velocity(phi_m.v.col(phi_m.get_index()).data(), phi_n.v.col(phi_n.get_index()).data(), second_indices.at(0)); 
					}
					
					if(second_indices.size() == 2){
						//std::cout << indices << "second index has 2\n";
						phi_n.Velocity2(phi_m.v.col(phi_m.get_index()).data(), phi_n.v.col(phi_n.get_index()).data(), second_indices.at(0), second_indices.at(1)); 
					}
					
					phi_m.Exchange_Boundaries(); 
					phi_m.template Multiply<0>();					
					gamma.matrix().block(n,0,1,2) += (phi0.v.adjoint() * phi_m.v - gamma.matrix().block(n,0,1,2))/value_type(average + 1);
					//std::cout << "(" << n << " " << 0 << ")\t" << gamma.matrix().block(n,0,1,2) << "\n";
					
					
					for(int m = 2; m < N_cheb_moments; m += 2)
					{
						phi_m.template Multiply<1>();
						phi_m.template Multiply<1>();
						
						gamma.matrix().block(n,m,1,2) += (phi0.v.adjoint() * phi_m.v - gamma.matrix().block(n,m,1,2))/value_type(average + 1);
						//std::cout << "(" << n << " " << m << ")\t" << gamma.matrix().block(n,m,1,2) << "\n";
					}
					
					// all done with |m>, now we return to |n>
					
					if(n == 0)
						phi_n.template Multiply<0>(); //first cheb iteration, multiply by the hamiltonian
					else if(n < N_cheb_moments - 1)
						phi_n.template Multiply<1>(); 
				}
				average++;
			}
		}
		
		//std::cout << "gamma_00\n" << gamma.matrix().block(0,0,3,1) << "\n";fflush(stdout);
		
		// the gamma matrix has been calculated. Now we're going to use the 
		// property that gamma is hermitian: gamma_nm=gamma_mn*
	
		
#pragma omp critical

		//std::cout << "IMPORTANT ! ! !:\n V is not hermitian. Make sure you take this into account\n";
		// in this case there's no problem. both V are anti-hermitic, so the minus signs cancel
		if(first_indices.size()==1 && second_indices.size()==1)
			Global.gamma.matrix() += (gamma.matrix() + gamma.matrix().adjoint())/2.0;
			
		// not here, though
		else
			Global.gamma.matrix() += (gamma.matrix() - gamma.matrix().adjoint())/2.0;
		
#pragma omp barrier

    
#pragma omp master
    {
      H5::H5File * file = new H5::H5File(name, H5F_ACC_RDWR);
			write_hdf5(Global.gamma, file, name_dataset);
      delete file;
      
      // make sure the global matrix is zeroed
      Global.gamma.setZero();
    }
#pragma omp barrier
  }



  void Measure_Lambda(int & NRandomV, int & NDisorder, std::string indices, std::string filename_dataset) {
		
		/* This function calculates Tr[ V^ab T_n ], the analogue in tight binding
		 * of the diamagnetic term. The calculation is almost identical with the
		 * calculation of the density of states (go figure). The only difference resides
		 * in multiplying the left-most phi0 by Velocity2 before doing the dot product
		 * with phi. Velocity2 = v^ab
		 * */
		 
		if(DEBUG)std::cout << "entered LAMBDA\n";fflush(stdout);
		
		// Process the indices
		std::vector<int> first_indices;
		first_indices.resize(indices.size()); 
		if(DEBUG)std::cout << "strings: " << indices << "\n";fflush(stdout);
		
		for(unsigned int i = 0; i < indices.size(); i++){
			if(indices[i]=='y')
				first_indices.at(i) = 1;
			else
				first_indices.at(i) = 0;
		}
		for(unsigned  int i = 0; i < indices.size(); i++)
			if(DEBUG)std::cout << "first_indices: " << first_indices.at(i) << "\n";fflush(stdout);
		
		
		
		//Done processing the indices
		Eigen::Array<T, -1, -1> lambda = Eigen::Array<T, -1, -1 > :: Zero(Global.lambda.rows(),Global.lambda.cols());
		
		typedef typename extract_value_type<T>::value_type value_type;


		//make sure the local matrix is zeroed to be used again
		lambda.setZero();
    
    KPM_Vector<T,D> phi0(1, *this);
    KPM_Vector<T,D>  phi(2, *this);

    
    long average = 0;
    for(int disorder = 0; disorder < NDisorder; disorder++)
		{
			h.distribute_AndersonDisorder();
			for(int randV = 0; randV < NRandomV; randV++)
				{
					phi0.initiate_vector();
					phi.set_index(0);	
					phi.v.col(0) = phi0.v.col(0);
					phi.Exchange_Boundaries();
					
					
					
					// check which operator needs to be calculated.
					if(first_indices.size()==1){
						phi.Velocity(phi0.v.col(0).data(), phi.v.col(phi.get_index()).data(), first_indices.at(0));
						phi0.empty_ghosts(0);
					}
					
					if(first_indices.size()==2){
						phi.Velocity2(phi0.v.col(0).data(), phi.v.col(phi.get_index()).data(), first_indices.at(0), first_indices.at(1));
						phi0.empty_ghosts(0);
					}
					
					
					phi.template Multiply<0>();
					lambda.matrix().block(0,0,1,2) +=  (phi0.v.adjoint() * phi.v - lambda.matrix().block(0,0,1,2))/value_type(average + 1);
					
					for(int m = 2; m < lambda.cols(); m += 2)
					{	    
						phi.template Multiply<1>();
						phi.template Multiply<1>();
						lambda.matrix().block(0,m,1,2) +=  (phi0.v.adjoint() * phi.v - lambda.matrix().block(0,m,1,2))/value_type(average + 1);
					}
					average++;
				}
    }
    
     
#pragma omp critical

		Global.lambda += lambda;
    
#pragma omp barrier

    
#pragma omp master
    {
      H5::H5File * file = new H5::H5File(name, H5F_ACC_RDWR);
			write_hdf5(Global.lambda, file, filename_dataset);
      delete file;
      
			//make sure the global matrix is zeroed to be used again later
			Global.lambda.setZero();
    }
#pragma omp barrier
  }

  void Single_Shot(double EScale, int & NRandomV, int & NDisorder, int N_cheb_moments, Eigen::Array<double, -1, 1> energy_array, double finite_gamma, std::string indices, std::string name_dataset) {
		/*
		 * Calculate the longitudinal conductivity for a single value of the energy
		 */
		 
		//std::cout << "energies: " << energy_array << "\n";
		 
		// Process the indices
		if(DEBUG)std::cout << "entered singleshot\n";
		std::vector<int> first_indices, second_indices;
		
		int comma_location = indices.find_first_of(',');
	
		std::string first_string, second_string;
		first_string = indices.substr(0, comma_location);
		second_string = indices.substr(comma_location+1);		
		first_indices.resize(first_string.size()); 
		second_indices.resize(second_string.size());
		
		if(DEBUG)std::cout << "strings: " << first_string << " "<< second_string << "\n";
		
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
			/*
		for(unsigned  int i = 0; i < first_string.size(); i++)
			std::cout << "first_indices: " << first_indices.at(i) << "\n";
			
		for(unsigned  int i = 0; i < second_string.size(); i++)
			std::cout << "second_indices: " << second_indices.at(i) << "\n";
		
		*/
		
		//Done processing the indices
		
		
		
		
		 
		
    typedef typename extract_value_type<T>::value_type value_type;
    KPM_Vector<T,D> phi0(1, *this);
    KPM_Vector<T,D> phi (2, *this);
    
    KPM_Vector<T,D> phi1(2, *this);
    KPM_Vector<T,D> phi2(2, *this);
		
	
		//int a = 4;
		//int b = 4;
		
    Eigen::Array<T, -1, -1> cond_array;
    int N_energies = energy_array.cols()*energy_array.rows(); //one of them is one. I don't know if it's the columns or the rows, but it doesn't matter
    //~ std::cout << "before cond array\n";fflush(stdout);
    cond_array = Eigen::Array<T, -1, -1>::Zero(1, Global.singleshot_cond.cols());
    
    if(DEBUG)std::cout << "gonna start calculating\n";fflush(stdout);
    
    for(int ener = 0; ener < N_energies; ener++){
			std::complex<double> energy(energy_array(ener), finite_gamma);
			if(DEBUG)std::cout << "finished setting complex energy\n";fflush(stdout);
			long average = 0;
			for(int disorder = 0; disorder < NDisorder; disorder++){
				if(DEBUG)std::cout << "before disorder\n";fflush(stdout);
					h.generate_disorder();
					if(DEBUG)std::cout << "after disorder\n";fflush(stdout);
					for(int randV = 0; randV < NRandomV; randV++){
						
						if(DEBUG)std::cout << "started calculating the first vector\n";fflush(stdout);
																		
						phi0.initiate_vector();					
						phi0.Exchange_Boundaries(); 	
						
						// calculate the left vector
						phi.set_index(0);				
						
						
						// |phi> = v |phi_0>
						phi.Velocity(phi.v.col(0).data(), phi0.v.col(0).data(), first_indices.at(0));	
						if(DEBUG)std::cout << "multiplied by velocity\n";fflush(stdout);
						phi.Exchange_Boundaries();	
						//~ phi1.v.col(0) = phi.v.col(phi.get_index())*(a==0);
						phi1.v.col(0) = phi.v.col(phi.get_index())*green(0, 1, energy).imag()/2.0;
									
						phi.template Multiply<0>();		
						//~ phi1.v.col(0) += phi.v.col(1)*(a==1);
						phi1.v.col(0) += phi.v.col(1)*green(1, 1, energy).imag();
						
						if(DEBUG)std::cout << "finished first cheb\n";fflush(stdout);
						
						for(int n = 2; n < N_cheb_moments; n++){		
							phi.template Multiply<1>();
							//~ phi1.v.col(0) += phi.v.col(phi.get_index())*(a==n);
							phi1.v.col(0) += phi.v.col(phi.get_index())*green(n, 1, energy).imag();
						}
						
						// multiply phi1 by the velocity operator again. 
						// We need a temporary vector to mediate the operation, which will be |phi>
						phi.v.col(0) = phi1.v.col(0);
						phi.Velocity(phi1.v.col(0).data(), phi.v.col(0).data(), second_indices.at(0));
						phi1.empty_ghosts(0);
						
						if(DEBUG)std::cout << "Finished calculating the first vector\n";fflush(stdout);
						
						// calculate the right vector
						phi.set_index(0);			
						phi.v.col(0) = phi0.v.col(0);
						
						//~ phi2.v.col(0) = phi.v.col(phi.get_index())*(b==0);
						phi2.v.col(0) = phi.v.col(phi.get_index())*green(0, 1, energy).imag()/2.0;
									
						phi.template Multiply<0>();		
						//~ phi2.v.col(0) += phi.v.col(1)*(b==1);
						phi2.v.col(0) += phi.v.col(1)*green(1, 1, energy).imag();
						
						for(int n = 2; n < N_cheb_moments; n++){		
							phi.template Multiply<1>();
							//~ phi2.v.col(0) += phi.v.col(phi.get_index())*(b==n);
							phi2.v.col(0) += phi.v.col(phi.get_index())*green(n, 1, energy).imag();
						}
						
						cond_array(ener) += (T(phi1.v.col(0).adjoint()*phi2.v.col(0)) - cond_array(ener))/value_type(average+1);						
						average++;
						
						if(DEBUG)std::cout << "Finished calculating the second vector\n";fflush(stdout);
				}
			}
		}
		
		//std::cout << "gamma_00\n" << gamma.matrix().block(0,0,3,1) << "\n";fflush(stdout);
		
		// the gamma matrix has been calculated. Now we're going to use the 
		// property that gamma is hermitian: gamma_nm=gamma_mn*
	
		
#pragma omp critical

		//std::cout << "IMPORTANT ! ! !:\n V is not hermitian. Make sure you take this into account\n";
		// in this case there's no problem. both V are anti-hermitic, so the minus signs cancel
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
      if(DEBUG)std::cout << "left singleshot\n";fflush(stdout);
    }
#pragma omp barrier
  }
	
  
};


  
