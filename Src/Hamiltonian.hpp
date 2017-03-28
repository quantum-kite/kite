#ifndef _HAMILTONIAN_HPP
#define _HAMILTONIAN_HPP

template <typename T, unsigned D>
class Hamiltonian {

public:
  KPMRandom <T> & rng;  
  LatticeStructure <D> & r;

  Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic> d;  // Distance
  Eigen::Matrix<   T, Eigen::Dynamic, Eigen::Dynamic> t;  // Hopping values
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
      
      d = Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic>(max, r.Orb);
      get_hdf5<int>(d.data(), file, (char *) "/Hamiltonian/d");

      t = Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic>(max, r.Orb);
      get_hdf5<T>(t.data(), file, (char *) "/Hamiltonian/d");
      
      for(unsigned  i = 0; i < max; i++)
	for(unsigned io = 0; io < r.Orb; io++)
	  {
	    b3.set_coord(d(i,io));	    
	    d(i,io) = (b3.coord[k] - io)* Ld.basis[D];
	    
	    for(unsigned k = 0; k < D; k++)
	      d(i,io) += (b3.coord[k] - 1) * Ld.basis[i];
	  }
      
      
      file->close();
    }
  };
};

#endif
