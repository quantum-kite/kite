template <typename T, unsigned D>
struct Periodic_Operator {
  Simulation<T,D> & simul;
  // Non-diagonal component of the operator
  Eigen::Array<unsigned, Eigen::Dynamic, 1 >         NHoppings;         // Number of elements different from Zero from each orbital
  Eigen::Array< int, Eigen::Dynamic, Eigen::Dynamic> distance;          // Distance in the basis 
  Eigen::Array<   T, Eigen::Dynamic, Eigen::Dynamic> hopping;           // Hopping
  Eigen::Array<   T, Eigen::Dynamic, Eigen::Dynamic> V[D];              // Velocity
  
  Periodic_Operator(Simulation<T,D> & sim) : simul(sim) {
    
    NHoppings =  Eigen::Array<unsigned, Eigen::Dynamic, 1 > (sim.r.Orb);
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(sim.name, H5F_ACC_RDONLY);
      get_hdf5<unsigned>(NHoppings.data(), file, (char *) "/Hamiltonian/NHoppings");

      int max  = NHoppings.maxCoeff();
      distance = Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic>  (max, sim.r.Orb);
      hopping    = Eigen::Matrix<   T, Eigen::Dynamic, Eigen::Dynamic>(max, sim.r.Orb);
      for(unsigned i = 0; i < D; i++)
        V[i]    = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(max, sim.r.Orb);
      
      get_hdf5<T>(hopping.data(), file, (char *) "/Hamiltonian/Hoppings");          // Read Hoppings
      get_hdf5<int>(distance.data(), file, (char *) "/Hamiltonian/d");              // Read the distances
      delete file;
      Convert_Build(sim.r);
    }
  }

  void Convert_Build (  LatticeStructure <D> & r )
  {
    Eigen::Matrix<double, D, 1> dr;
    unsigned l[D + 1];
    std::fill_n(l, D, 3);
    l[D]  = r.Orb;
    
    Coordinates<int, D + 1> b3(l), Ld(r.Ld);    
    Eigen::Map<Eigen::Matrix<int,D, 1>> v(b3.coord);
    
    for(unsigned io = 0; io < r.Orb; io++)
      for(unsigned  i = 0; i < NHoppings(io); i++)
        {
          b3.set_coord(distance(i,io));                              // Get Coordinates in Basis 3,  The last is the Final Orbital.
          b3.coord[D] -= io;                                         // Obtain the difference in orbital coordinates
          v.array() -= 1;                                            // Subtract to the first D elements of v to get v(i) in (-1, 0 , 1)
          
          distance(i,io) = Ld.set_index(b3.coord).index;             // Convert in distances in this lattice
          
          dr = r.rOrb.col(io) - r.rOrb.col(b3.coord[D] + io) ;       // The D components of the vector difference in orbital positions      
          dr += r.rLat * v.template cast<double>();
          
          for(unsigned dim = 0; dim < D; dim++)
            V[dim](i,io) = hopping(i,io) * T( dr(dim) );
	  
        }
  };
};
