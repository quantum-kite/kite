#ifndef _HAMILTONIAN_HPP
#define _HAMILTONIAN_HPP

template <typename T, unsigned D>
class Hamiltonian {

public:
  KPMRandom <T> & rng;  
  LatticeStructure <D> & r;
  
  Eigen::Array< int, Eigen::Dynamic, Eigen::Dynamic> d;          // Distance
  Eigen::Array<   T, Eigen::Dynamic, Eigen::Dynamic> t[D + 1];   // Hopping values
  Eigen::Array<   T, Eigen::Dynamic, Eigen::Dynamic> U ;   // Hopping values
  Eigen::Array<unsigned, Eigen::Dynamic, 1 >         NHoppings;
  
  Hamiltonian (KPMRandom <T> & rng1, LatticeStructure <D> & r1,  char *name) : rng(rng1), r(r1) {
    Eigen::Matrix<double, D, 1> dr;
    unsigned l[D + 1];
    NHoppings =  Eigen::Array<unsigned, Eigen::Dynamic, 1 > (r.Orb);
    std::fill_n(l, D, 3);
    l[D]  = r.Orb;
    Coordinates<int, D + 1> b3(l), Ld(r.Ld);    
    Eigen::Map<Eigen::Matrix<int,D, 1>> v(b3.coord);
	    
#pragma omp critical
    {        
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      get_hdf5<unsigned>(NHoppings.data(), file, (char *) "/Hamiltonian/NHoppings");
      
      // Allocate distances, Hamiltonian and velocities : t[0] --> Hopping, t[1] --> Velocity component 0 ...
      d = Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic>(NHoppings.maxCoeff(), r.Orb);
      for(int i = 0; i < int(D + 1); i++)
	t[i] = Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic>(NHoppings.maxCoeff(), r.Orb);
      
      // Read Hoppings
      get_hdf5<T>(t[0].data(), file, (char *) "/Hamiltonian/Hoppings");
      
      // Read the distances
      get_hdf5<int>(d.data(), file, (char *) "/Hamiltonian/d");

      // convert the distances
      
      for(unsigned io = 0; io < r.Orb; io++)
	for(unsigned  i = 0; i < NHoppings[io]; i++)
	  {
	    b3.set_coord(d(i,io));                              // Get Coordintes in Basis 3,  The last is the Final Orbital.
	    b3.coord[D] -= io;                                  // Obtain the difference in orbital coordinates
	    v.array() -= 1; 	                                // Subtract to the first D elements of v to get v(i) in (-1, 0 , 1)
	    
	    d(i,io) = Ld.set_index(b3.coord).index;             // Convert in distances in this lattice
	    dr = r.rOrb.col(io) - r.rOrb.col(b3.coord[D] + io) ; // The D components of the vector difference in orbital positions	    

	    dr += r.rLat * v.template cast<double>();
	    
	    for(unsigned dim = 0; dim < D; dim++)
	      t[1 + dim](i,io) = t[0](i,io) * T( dr(dim) );	      
	  }
      U = Eigen::Array<   T, Eigen::Dynamic, Eigen::Dynamic> (r.Sized, 1);
      for(unsigned i = 0; i < r.Sized; i++)
	U(i,0) = rng.uniform( double(0.), double(0.1) );
      file->close();
    }
  };  
};

#endif
