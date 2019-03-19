
//template <typename T,unsigned D>
//void Simulation<T,D>::Measure_Gamma(measurement_queue queue) {
  //debug_message("Entered Measure_Gamma.\n");
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
  //int NRandomV = queue.NRandom;
  //int NDisorder = queue.NDisorder;
  //std::vector<int> N_moments = queue.NMoments;
  //std::string indices_string = queue.direction_string;
  //std::string name_dataset = queue.label;
 
  //debug_message("Indices: "); debug_message(indices_string); debug_message(".\n");
		


  //// First of all, we need to process the indices_string into something the program can use
//#pragma omp barrier
  //std::vector<std::vector<unsigned>> indices = process_string(indices_string);
  //int dim = indices.size();
		
		
		
  //// Check if the dimensions match
  //if(dim != int(N_moments.size())){
    //std::cout << "Dimension of the Gamma matrix does not match the number of chebyshev moments. Aborting.\n";
    //exit(1);
  //}
    
  //if(dim == 2 and MEMORY > 2){
    //if(N_moments.at(0)%MEMORY!=0 or N_moments.at(1)%MEMORY!=0){
      //std::cout << "The number of Chebyshev moments ("<< N_moments.at(0)<<","<< N_moments.at(1)<<")"; 
      //std::cout << "has to be a multiple of MEMORY ("<< MEMORY <<"). Exiting.\n";
      //exit(1);
    //}
    //Gamma2D(NRandomV, NDisorder, N_moments, indices, name_dataset);
  //} else {
    //if(dim == 3){
      //Gamma3D(NRandomV, NDisorder, N_moments, indices, name_dataset);
    //} else {
      //Gamma1D(NRandomV, NDisorder, N_moments.at(0), indices, name_dataset);
      ////GammaGeneral(NRandomV, NDisorder, N_moments, indices, name_dataset);
    //}
  //}
  
  
  //debug_message("Left Measure_Gamma\n");
//};


//template <typename T,unsigned D>
//void Simulation<T,D>::recursive_KPM(int depth, int max_depth, std::vector<int> N_moments, long *average, long *index_gamma, 
                   //std::vector<std::vector<unsigned>> indices, std::vector<KPM_Vector<T,D>*> *kpm_vector, Eigen::Array<T, -1, -1> *gamma){
  //debug_message("Entered recursive_KPM\n");
  //typedef typename extract_value_type<T>::value_type value_type;
  //Eigen::Matrix < T, 1, 2> tmp =  Eigen::Matrix < T, 1, 2> ::Zero();		
		
		
  //if(depth != max_depth){
    //KPM_Vector<T,D> *kpm1 = kpm_vector->at(depth);
    //KPM_Vector<T,D> *kpm2 = kpm_vector->at(depth + 1);
			
    //T * kpm1data;
    //T * kpm2data;
			
    ////std::cout << "first branch. Depth: " << depth << "\n" << std::flush;
			
			
    //for(int p = 0; p < N_moments.at(depth - 1); p++){
      //kpm2->set_index(0);
      //switch(indices.at(depth).size()){
      //case 0:
        //break;
      //default:
        //kpm1data = kpm1->v.col(kpm1->get_index()).data();
        //kpm2data = kpm2->v.col(kpm2->get_index()).data();
                  
        //kpm2->Velocity(kpm2data, kpm1data, max_depth - depth); 											
      //}
				
      //recursive_KPM(depth + 1, max_depth, N_moments, average, index_gamma, indices, kpm_vector, gamma);
      //if(p == 0){
        //kpm1->template Multiply<0>(); 
      //}
      //else if(p < N_moments.at(depth-1) - 1){
        //kpm1->template Multiply<1>(); 
      //}
    //}
			
  //} else {
    //KPM_Vector<T,D> *kpm0 = kpm_vector->at(0);
    //KPM_Vector<T,D> *kpm1 = kpm_vector->at(depth);
			
    //kpm1->template Multiply<0>();		
    //tmp.setZero();
    //for(std::size_t ii = 0; ii < r.Sized ; ii += r.Ld[0])
      //tmp += kpm0->v.block(ii,0, r.Ld[0], 1).adjoint() * kpm1->v.block(ii, 0, r.Ld[0], 2);
    //gamma->matrix().block(0,*index_gamma,1,2) += (tmp - gamma->matrix().block(0,*index_gamma,1,2))/value_type(*average + 1);			
    //*index_gamma += 2;
	
    //for(int m = 2; m < N_moments.at(depth - 1); m += 2){
      //kpm1->template Multiply<1>();
      //kpm1->template Multiply<1>();
      //tmp.setZero();
      //for(std::size_t ii = 0; ii < r.Sized ; ii += r.Ld[0])
        //tmp += kpm0->v.block(ii,0, r.Ld[0], 1).adjoint() * kpm1->v.block(ii, 0, r.Ld[0], 2);
      //gamma->matrix().block(0, *index_gamma,1,2) += (tmp - gamma->matrix().block(0,*index_gamma,1,2))/value_type(*average + 1);
            
      //*index_gamma += 2;
    //}
  //}
		
  //debug_message("Left recursive_KPM\n");
//};
//
//



//template <typename T,unsigned D>
//void Simulation<T,D>::GammaGeneral(int NRandomV, int NDisorder, std::vector<int> N_moments,
                                   //std::vector<std::vector<unsigned>> indices, std::string name_dataset){
    
  //int dim = indices.size();
		
  //// Check if the dimensions match
  //if(dim != int(N_moments.size())){
    //std::cout << "Dimension of the Gamma matrix does not match the number of chebyshev moments. Aborting.\n";
    //exit(1);
  //}
			
  //// Determine the size of the gamma matrix we want to calculate
  //int size_gamma = 1;
  //for(int i = 0; i < dim; i++){
    //if(N_moments.at(i) % 2 != 0){
      //std::cout << "The number of moments must be an even number, due to limitations of the program. Aborting\n";
      //exit(1);
    //}
    //size_gamma *= N_moments.at(i);
  //}
		
		
  //// Initialize the KPM vectors that will be needed to run the program 
  //std::vector<KPM_Vector<T,D>*> kpm_vector(dim+1);
  //kpm_vector.at(0) = new KPM_Vector<T,D> (1, *this);
  //for(int i = 0; i < dim; i++)
    //kpm_vector.at(i+1) = new KPM_Vector<T,D> (2, *this);
		
		
  //// Define some pointers to make the code easier to read
  //KPM_Vector<T,D> *kpm0 = kpm_vector.at(0);
  //KPM_Vector<T,D> *kpm1 = kpm_vector.at(1);
			
  //// Make sure the local gamma matrix is zeroed
  //Eigen::Array<T, -1, -1> gamma = Eigen::Array<T, -1, -1 >::Zero(1, size_gamma);

  //long average = 0;
  //for(int disorder = 0; disorder < NDisorder; disorder++){
    //h.generate_disorder();
    //for(unsigned it = 0; it < indices.size(); it++)
      //h.build_velocity(indices.at(it), it);

    //for(int randV = 0; randV < NRandomV; randV++){
        
      //kpm0->initiate_vector();			// original random vector
      //kpm1->set_index(0);
      //kpm1->v.col(0) = kpm0->v.col(0);
      //kpm1->Exchange_Boundaries();

      //// replace <0| by  <0|v. Note that v is not self-adjoint in this formulation
      //generalized_velocity(kpm0, kpm1, indices, 0);
      //int factor = 1 - (indices.at(0).size() % 2)*2;		
      //kpm0->v.col(0) = factor*kpm0->v.col(0); // This factor is due to the fact that this Velocity operator is not self-adjoint

      //kpm0->empty_ghosts(0);
      //long index_gamma = 0;
      //recursive_KPM(1, dim, N_moments, &average, &index_gamma, indices, &kpm_vector, &gamma);
      //average++;
    //}
  //} 
		
		
  //store_gamma(&gamma, N_moments, indices, name_dataset);
		
  //// delete the kpm_vector
  //delete kpm_vector.at(0);
  //for(int i = 0; i < dim; i++)
    //delete kpm_vector.at(i+1);
	
//}

//template <typename T,unsigned D>
//void Simulation<T,D>::GammaGeneral(int NRandomV, int NDisorder, std::vector<int> N_moments,
                                   //std::vector<std::vector<unsigned>> indices, std::string name_dataset){
    
  //int dim = indices.size();
		
  //// Check if the dimensions match
  //if(dim != int(N_moments.size())){
    //std::cout << "Dimension of the Gamma matrix does not match the number of chebyshev moments. Aborting.\n";
    //exit(1);
  //}
			
  //// Determine the size of the gamma matrix we want to calculate
  //int size_gamma = 1;
  //for(int i = 0; i < dim; i++){
    //if(N_moments.at(i) % 2 != 0){
      //std::cout << "The number of moments must be an even number, due to limitations of the program. Aborting\n";
      //exit(1);
    //}
    //size_gamma *= N_moments.at(i);
  //}
		
		
  //// Initialize the KPM vectors that will be needed to run the program 
  //std::vector<KPM_Vector<T,D>*> kpm_vector(dim+1);
  //kpm_vector.at(0) = new KPM_Vector<T,D> (1, *this);
  //for(int i = 0; i < dim; i++)
    //kpm_vector.at(i+1) = new KPM_Vector<T,D> (2, *this);
		
		
  //// Define some pointers to make the code easier to read
  //KPM_Vector<T,D> *kpm0 = kpm_vector.at(0);
  //KPM_Vector<T,D> *kpm1 = kpm_vector.at(1);
			
  //// Make sure the local gamma matrix is zeroed
  //Eigen::Array<T, -1, -1> gamma = Eigen::Array<T, -1, -1 >::Zero(1, size_gamma);

  //long average = 0;
  //for(int disorder = 0; disorder < NDisorder; disorder++){
    //h.generate_disorder();
    //for(unsigned it = 0; it < indices.size(); it++)
      //h.build_velocity(indices.at(it), it);

    //for(int randV = 0; randV < NRandomV; randV++){
        
      //kpm0->initiate_vector();			// original random vector
      //kpm1->set_index(0);
      //kpm1->v.col(0) = kpm0->v.col(0);
      //kpm1->Exchange_Boundaries();

      //// replace <0| by  <0|v. Note that v is not self-adjoint in this formulation
      //generalized_velocity(kpm0, kpm1, indices, 0);
      //int factor = 1 - (indices.at(0).size() % 2)*2;		
      //kpm0->v.col(0) = factor*kpm0->v.col(0); // This factor is due to the fact that this Velocity operator is not self-adjoint

      //kpm0->empty_ghosts(0);
      //long index_gamma = 0;
      //recursive_KPM(1, dim, N_moments, &average, &index_gamma, indices, &kpm_vector, &gamma);
      //average++;
    //}
  //} 
		
		
  //store_gamma(&gamma, N_moments, indices, name_dataset);
		
  //// delete the kpm_vector
  //delete kpm_vector.at(0);
  //for(int i = 0; i < dim; i++)
    //delete kpm_vector.at(i+1);
	
//}

