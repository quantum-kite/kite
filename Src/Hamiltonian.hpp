#ifndef _HAMILTONIAN_HPP
#define _HAMILTONIAN_HPP

template <typename T, unsigned D>
class Hamiltonian {

public:
  KPMRandom <T> & rng;  
  LatticeStructure <D> & r;
  
  Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic> d;  // Distance
  Eigen::Matrix<   T, Eigen::Dynamic, Eigen::Dynamic> t[D + 1];  // Hopping values
  unsigned * NHoppings;
  
  Hamiltonian (KPMRandom <T> & rng1, LatticeStructure <D> & r1,  char *name) : rng(rng1), r(r1) {
    unsigned l[D + 1], max = 0;
    for(unsigned i = 0; i < D; i++) l[i] = 3;
    l[D] = r.Orb;
    
    Coordinates<int, D + 1> b3(l), Ld(r.Ld);    

#pragma omp critical
    {        
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      NHoppings = new unsigned [r.Orb];
      get_hdf5<unsigned>(NHoppings, file, (char *) "/Hamiltonian/NHoppings");

      for(unsigned i = 0; i < r.Orb; i++)
	max = (max < NHoppings[i] ? NHoppings[i] : max );


      // Read The Hamiltonian 
      for(int i = 0; i < int(D + 1); i++)
	t[i] = Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic>(max, r.Orb);
      get_hdf5<T>(t[0].data(), file, (char *) "/Hamiltonian/Hoppings");
      
      // Read the distances
      d = Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic>(max, r.Orb);
      get_hdf5<int>(d.data(), file, (char *) "/Hamiltonian/d");
#pragma omp master
      {
	std::cout << " max: " << max << std::endl;
	std::cout << " d: " << std::endl << d << std::endl;
      }
      // convert the distances

      for(unsigned io = 0; io < r.Orb; io++)
	for(unsigned  i = 0; i < NHoppings[io]; i++)
	  {
	    
	    b3.set_coord(d(i,io));
	    b3.coord[D] -= int(io);
	    
	    for(unsigned k = 0; k < D; k++)
	      b3.coord[k] -= 1;

	    d(i,io) = Ld.set_index(b3.coord).index;
	    
	    // Build Velocity
	    for(unsigned dim = 0; dim < D; dim++)
	      {
		t[1 + dim](i,io) = t[0](i,io) * T( r.rOrb[b3.coord[D] + int(io)][dim] - r.rOrb[io][dim] );
		for(unsigned k = 0; k < D; k++)
		  t[1 + dim](i,io) += t[0](i,io) * T(b3.coord[k] * r.rLat[k][dim]);
	      }
	  }

      file->close();
    }
  };

};

#endif
