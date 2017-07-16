#ifndef _HAMILTONIAN_HPP
#define _HAMILTONIAN_HPP


template <typename T>
struct Periodic_Operator {
  Eigen::Array<unsigned, Eigen::Dynamic, 1 >         NHoppings;         // Number of elements different from Zero from each orbital
  Eigen::Array< int, Eigen::Dynamic, Eigen::Dynamic> distance;          // Distance in the basis 
  Eigen::Array<   T, Eigen::Dynamic, Eigen::Dynamic> value;             // operator values

  Periodic_Operator(int orb) {
    NHoppings =  Eigen::Array<unsigned, Eigen::Dynamic, 1 > (orb);    
  }
  
  void finishbuild() {
    // Once we know NHoppings we can allocate the full matrix
    int max  = NHoppings.maxCoeff();
    int orb = NHoppings.rows();
    distance = Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic>(max, orb);
    value    = Eigen::Matrix<   T, Eigen::Dynamic, Eigen::Dynamic>(max, orb);
  }
};

template <typename T>
struct Defect_Operator  {
  std::vector <int> position;                                           // vector with position in the lattice of the Orbital 0 of the defects
  int NumberNodes;                                                      // Number of nodes
  std::vector <int> node_position;                                      // Relative distances of the nodes to the defect position
  std::vector <T>   value;                                              // vector of the non-zero values of the operator
  std::vector <int> element1;                                           // vector with the nodes 
  std::vector <int> element2;                                           // vector with the nodes

  /* Broken patterns */
  
  std::vector <T>   border_value;                                      // value of the operator
  std::vector <int> border_element1_position;                          // positions of the memory elements involved
  std::vector <int> border_element2_position;                          // positions of the memory elements involved
};


template <typename T, unsigned D>
class Hamiltonian {
public:
  typedef typename extract_value_type<T>::value_type value_type;
  KPMRandom <T> & rng;  
  LatticeStructure <D> & r;
  Periodic_Operator<T> hr;

  std::vector <int> orb_num;
  std::vector <int> model;
  std::vector <double> mu;
  std::vector <double> sigma;
  
  std::vector<value_type> U_Orbital;
  std::vector<value_type> U_Anderson;      // Local disorder  
  std::vector<int> Anderson_orb_address;
  

  Hamiltonian (KPMRandom <T> & rng1, LatticeStructure <D> & r1,  char *name) : rng(rng1), r(r1), hr(r1.Orb) {
    Eigen::Matrix<double, D, 1> dr;
    unsigned l[D + 1];
    std::fill_n(l, D, 3);
    l[D]  = r.Orb;
    Coordinates<int, D + 1> b3(l), Ld(r.Ld);    
    Eigen::Map<Eigen::Matrix<int,D, 1>> v(b3.coord);
	    
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      get_hdf5<unsigned>(hr.NHoppings.data(), file, (char *) "/Hamiltonian/NHoppings");
      hr.finishbuild();
      
      get_hdf5<T>(hr.value.data(), file, (char *) "/Hamiltonian/Hoppings");          // Read Hoppings
      get_hdf5<int>(hr.distance.data(), file, (char *) "/Hamiltonian/d");            // Read the distances

      // convert the distances
      
      for(unsigned io = 0; io < r.Orb; io++)
	for(unsigned  i = 0; i < hr.NHoppings(io); i++)
	  {
	    b3.set_coord(hr.distance(i,io));                              // Get Coordinates in Basis 3,  The last is the Final Orbital.
	    b3.coord[D] -= io;                                            // Obtain the difference in orbital coordinates
	    v.array() -= 1; 	                                          // Subtract to the first D elements of v to get v(i) in (-1, 0 , 1)
	    hr.distance(i,io) = Ld.set_index(b3.coord).index;             // Convert in distances in this lattice
	    
	    //	    dr = r.rOrb.col(io) - r.rOrb.col(b3.coord[D] + io) ; // The D components of the vector difference in orbital positions	    
	    //	    dr += r.rLat * v.template cast<double>();
	    
	    //	    for(unsigned dim = 0; dim < D; dim++)
	    //	      t[1 + dim](i,io) = t[0](i,io) * T( dr(dim) );	      
	  }
      
      build_Anderson_disorder(file);
      file->close();
    }
    
    distribute_AndersonDisorder();    
  };

  void build_Anderson_disorder(H5::H5File * file)
  {
    H5::DataSet   dataset    = H5::DataSet(file->openDataSet("/Hamiltonian/Disorder/OrbitalNum"));
    H5::DataSpace dataspace  = dataset.getSpace();
    size_t        m          = dataspace.getSimpleExtentNpoints();
    orb_num.resize(m);
    model.resize(m);
    mu.resize(m);
    sigma.resize(m);

    get_hdf5<int>(orb_num.data(), file, (char *) "/Hamiltonian/Disorder/OrbitalNum");               // read the orbitals that have local disorder
    get_hdf5<int> (model.data(), file, (char *) "/Hamiltonian/Disorder/OnsiteDisorderModelType");   // read the the type        of local disorder
    get_hdf5<double> (mu.data(), file, (char *) "/Hamiltonian/Disorder/OnsiteDisorderMeanValue");   // read the the mean value
    get_hdf5<double> (sigma.data(), file, (char *) "/Hamiltonian/Disorder/OnsiteDisorderMeanStdv"); // read the the variance
    
    Anderson_orb_address.resize(r.Orb);
    U_Orbital.resize(r.Orb);
    
    std::fill_n ( Anderson_orb_address.begin(), r.Orb, -2 );
    std::fill_n ( U_Orbital.begin(), r.Orb,  0 );
    /*
     * Gaussian      : 1
     * Uniform       : 2
     * Deterministic : 3
     */
    int sum = 0;
    for (unsigned i = 0; i < model.size(); i++)
      {
	if(model.at(i) < 3)
	  {
	    Anderson_orb_address.at( orb_num.at(i) ) = sum;
	    sum++;
	  }
	
	if(model.at(i) == 3) // Deterministic
	  {
	    Anderson_orb_address.at( orb_num.at(i) ) = -1;
	    U_Orbital.at( orb_num.at(i) ) = mu.at(i);
	  }
	
      }
    
    if(sum > 0)
      U_Anderson.resize( sum * r.Nd);    
  }
  
  
  void distribute_AndersonDisorder()
  {
    int sum = 0;
    for (unsigned i = 0; i < model.size(); i++)
      if(model.at(i) == 1 )
	{
	  for(unsigned j = 0; j < r.Nd; j++)
	    U_Anderson.at(sum + j) = rng.gaussian(mu.at(i), sigma.at(i)) ;
	  sum += r.Nd;
	}
      else if ( model.at(i) == 2 )
	{
	  for(unsigned j = 0; j < r.Nd; j++)
	    U_Anderson.at(sum + j) = rng.uniform(mu.at(i), sigma.at(i)) ;
	  sum += r.Nd;
	}
  }
};




#endif
