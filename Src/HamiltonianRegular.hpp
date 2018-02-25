#ifndef DEBUG
#define DEBUG 0
#endif

template <typename T, unsigned D>
struct Periodic_Operator {
  Simulation<T,D> & simul;
  // Non-diagonal component of the operator
  Eigen::Array<unsigned,    Eigen::Dynamic, 1 >            NHoppings;         // Number of elements different from Zero from each orbital
  Eigen::Array<std::ptrdiff_t, Eigen::Dynamic, Eigen::Dynamic> distance;          // Distance in the basis 
  Eigen::Array<   T, Eigen::Dynamic, Eigen::Dynamic> hopping;           // Hopping
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> dist;
  Eigen::Array<   T, Eigen::Dynamic, Eigen::Dynamic> hopping_magnetic; 	// Hoppings with magnetic field, periodic part
  Eigen::Array<   T, Eigen::Dynamic, Eigen::Dynamic> V[D];              // Velocity [r,h]
  Eigen::Array<   T, Eigen::Dynamic, Eigen::Dynamic> V2[D][D];              // Velocity [r,[r,h]]
  
  Periodic_Operator(Simulation<T,D> & sim) : simul(sim) {
    debug_message("Entered Periodic_Operator constructor.\n");
    
    NHoppings =  Eigen::Array<unsigned, Eigen::Dynamic, 1 > (sim.r.Orb);
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(sim.name, H5F_ACC_RDONLY);
      get_hdf5<unsigned>(NHoppings.data(), file, (char *) "/Hamiltonian/NHoppings");

      std::size_t max  	= NHoppings.maxCoeff();
      distance 			= Eigen::Matrix< std::ptrdiff_t, Eigen::Dynamic, Eigen::Dynamic>  (max, sim.r.Orb);
      dist 				= Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic>  (max, sim.r.Orb);
      hopping  			= Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(max, sim.r.Orb);
      hopping_magnetic  = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(max, sim.r.Orb);
      
      
      for(unsigned i = 0; i < D; i++)
        V[i]    = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(max, sim.r.Orb);
       
      
      for(unsigned i = 0; i < D; i++)
				for(unsigned j = 0; j < D; j++)
					V2[i][j] = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(max, sim.r.Orb);

      get_hdf5<T>(hopping.data(), file, (char *) "/Hamiltonian/Hoppings");          // Read Hoppings
      get_hdf5<int>(dist.data(), file, (char *) "/Hamiltonian/d");                  // Read the distances
      
      for(std::size_t i = 0; i < max; i++ )
	for(std::size_t j = 0; j <  sim.r.Orb; j++ )      
	  distance(i,j) = dist(i,j);
      delete file;
      Convert_Build(sim.r);
    }
    debug_message("Left Periodic_Operator constructor.\n");
  }







  void Convert_Build (  LatticeStructure <D> & r )
  {
	debug_message("Entered Convert_Build\n");
	  
    Eigen::Matrix<double, D, 1> dr; 
    Eigen::Matrix<double, D, 1> dr_a; //difference vector in coordinates of the lattice
    
    // These quantities will be needed to calculate the phase due to the magnetic field
    Eigen::Matrix<double, D, 1> orbital_difference_R; //vectors in real space
    Eigen::Matrix<double, D, 1> lattice_difference_R;
    Eigen::Matrix<double, D, 1> orbital_difference_a; //vectors in coordinates of the lattice
    Eigen::Matrix<double, D, 1> lattice_difference_a;
    
    
    // l = (3,3,Number_of_orbitals)
    unsigned l[D + 1];
    std::fill_n(l, D, 3);
    l[D]  = r.Orb;
    
    // rLat is organized per columns each column is a lattice vector
    
    Coordinates<std::ptrdiff_t, D + 1> b3(l), Ld(r.Ld);    
    Eigen::Map<Eigen::Matrix<std::ptrdiff_t,D, 1>> v(b3.coord); // Column vector
    
    for(unsigned io = 0; io < r.Orb; io++)
      for(unsigned  i = 0; i < NHoppings(io); i++)
        {
          b3.set_coord(distance(i,io));                              // Get Coordinates in Basis 3,  The last is the Final Orbital.
          b3.coord[D] -= io;                                         // Obtain the difference in orbital coordinates
          v.array() -= 1;                                            // Subtract to the first D elements of v to get v(i) in (-1, 0 , 1)
          distance(i,io) = Ld.set_index(b3.coord).index;             // Convert in distances in this lattice
          
          
          // difference vectors in real-space coordinates
          orbital_difference_R = r.rOrb.col(b3.coord[D] + io) - r.rOrb.col(io)  ;       // The D components of the vector difference in orbital positions in real space
          lattice_difference_R = r.rLat * v.template cast<double>();
          dr = orbital_difference_R + lattice_difference_R;
		  
          // difference vectors in lattice coordinates
          orbital_difference_a = r.rLat.inverse() * orbital_difference_R;       // vector difference in orbital positions, expressed in terms of the lattice vectors      
          lattice_difference_a = r.rLat.inverse() * lattice_difference_R;
          dr_a = orbital_difference_a + lattice_difference_a;
          
          
          // periodic part of the Peierls phase. 
          // If the kpm vector is not complex, peierls1(phase) will return 1.0, so it is effectively harmless.
          //double phase = ((0.5*dr_a.transpose() + (r.rLat.inverse()*r.rOrb.col(io)).transpose())*r.vect_pot*dr_a - lattice_difference_a.transpose()*r.vect_pot*r.rLat.inverse()*r.rOrb.col(b3.coord[D] + io))(0,0);
          double phase = ((-0.5*dr_a.transpose() + (r.rLat.inverse()*r.rOrb.col(b3.coord[D] + io)).transpose())*r.vect_pot*(-dr_a) + lattice_difference_a.transpose()*r.vect_pot*r.rLat.inverse()*r.rOrb.col(io))(0,0);
          
          
          hopping_magnetic(i,io) = hopping(i,io) * peierls1(phase);
          hopping(i,io) = hopping_magnetic(i,io);
          
          //if(DEBUG) std::cout << "Finished calculating the phase.\n" << std::flush;
          
          for(unsigned dim = 0; dim < D; dim++)
            V[dim](i,io) = hopping(i,io) * T( dr(dim) );
          
            
          for(unsigned dim1 = 0; dim1 < D; dim1++)
            for(unsigned dim2 = 0; dim2 < D; dim2++)
              V2[dim1][dim2](i,io) = hopping(i,io) * T( dr(dim1) )* T( dr(dim2) );
              
          
						
	  
        }
    debug_message("Left Convert_Build.\n");
  };
	
  template <typename U = T>
  typename std::enable_if<is_tt<std::complex, U>::value, U>::type peierls1(double phase) {
	std::complex<double> im(0,1.0);
    return U(exp(im*phase));
  };
  
  template <typename U = T>
  typename std::enable_if<!is_tt<std::complex, U>::value, U>::type peierls1(double phase) {
    return 1.0;
  };
};
