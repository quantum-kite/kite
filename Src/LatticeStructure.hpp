/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/


template <typename T, unsigned D>
struct Coordinates {
  T index = 0;
  T coord[D] = {};
  unsigned (&L)[D];
  T basis[D];
  
  Coordinates(T x, unsigned (&b)[D] ) : L(b) {
    buildBasis();
    set_coord(x);
  };
  
  Coordinates(T (&coord)[D], unsigned (&b)[D] ) : L(b) {
    buildBasis();
    set_index(coord);
  };
  
  Coordinates(unsigned (&b)[D] ) : L(b) {
    buildBasis();
    index = 0;
  };
  
  void buildBasis() {
    T p = 1;
    
    for(unsigned i = 0; i < D; i++)
	{
		basis[i] = p;
		p *= T(L[i]);
	}
  }
  void print() {
    for(unsigned i = 0; i < D ; i++)
      std::cout << coord[i] << " ";
    std::cout << std::endl;
  }
  
  template <typename T1>
  Coordinates & set(std::initializer_list<T1> a_args)  {
    int k = 0;
    for (auto i: a_args) coord[k++] = T(i);
    set_index(coord);
    return *this;
  }

  template <typename T1>
  Coordinates & set_index(T1 (&c)[D] ) {
    index = 0;
    for(int i = D - 1; i >= 0; i--)
	{
		coord[i] = T(c[i]);
		index += T(c[i]) * basis[i];
	}
    return *this;
  }
  
  Coordinates & set_coord(T x ) {
    index = x;
    for(int i = D - 1; i >= 0; i--)
      {
	coord[i] = x / basis[i];
	x = x % basis[i];
      };
    return *this;
  }
  
  Coordinates & add( Coordinates<T,D> & x) {
    for(int i = 0; i < int(D); i++)
      coord[i] = (L[i] + coord[i] + x.coord[i]) % L[i];

    set_index(coord);
    return *this;
  }

  Coordinates & subtract( Coordinates<T,D> & x) {
    for(int i = 0; i < int(D); i++)
      coord[i] = coord[i] - x.coord[i];
    set_index(coord);
    return *this;
  }
};


template <unsigned D>
struct LatticeStructure {
private:
  
public:
  std::size_t Sized; // Size of vector for each subdomain (with ghosts)
  std::size_t Size; // Size of vector for each subdomain (without ghosts)
  std::size_t Sizet; // Size of full Hilbert Space
  std::size_t SizetVacancies;
  Eigen::Matrix<double, D, D> rLat;  // The vectors are organized by columns 
  Eigen::MatrixXd rOrb;              // The vectors of each orbital are organized by columns
  unsigned nd[D + 1]; // Number of domains in each dimension (the last dimension corresponding with Orbitals are not decomposed)
  unsigned n_threads; // Number of threads
  
  unsigned Lt[D+1]; // Dimensions of the global sample
  unsigned Ld[D+1]; // Dimensions of each sub-domain (domain  + ghosts) 
  unsigned ld[D+1]; // Dimensions of each sub-domain (domain) 
  unsigned Bd[D+1]; // Information about periodic or non-periodic boundary conditions
  unsigned lStr[D + 1];
  std::size_t Nt; // Number of lattice postions of the global sample
  std::size_t Nd; // Number of lattice postions of the sub-domain with ghosts
  std::size_t N; // Number of lattice postions of the sub-domain without ghosts
  std::size_t NStr; // Number of lattice postions of the sub-domain without ghosts
  unsigned Orb; // Number of orbitals
  unsigned thread_id; // thread identification
  int MagneticField;

  
  Eigen::Matrix<double, D, D> vect_pot; // vector potential
  
  LatticeStructure(char *name ) {
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      get_hdf5<unsigned>(&Orb, file, (char *) "/NOrbitals");
      get_hdf5<double>(rLat.data(), file, (char *) "/LattVectors");    
      rOrb = Eigen::MatrixXd::Zero(D, Orb);
      get_hdf5<double>(rOrb.data(), file, (char *) "/OrbPositions");
      
      get_hdf5<unsigned>(Lt, file, (char *) "/L");
      get_hdf5<unsigned>(Bd, file, (char *) "/Boundaries");
      get_hdf5<unsigned>(nd, file, (char *) "/Divisions");      
      
      // verify if the number of boundaries divides the length
      if(Lt[0]%nd[0] != 0){
        std::cout << "The number of divisions in the x direction ("<< nd[0] <<") must ";
        std::cout << "be a divisor of the length of that side ("<< Lt[0] <<"). Exiting.\n";
        exit(1);
      }
      if(Lt[1]%nd[1] != 0){
        std::cout << "The number of divisions in the y direction ("<< nd[1] <<") must ";
        std::cout << "be a divisor of the length of that side ("<< Lt[1] <<"). Exiting.\n";
        exit(1);
      }



      MagneticField = 0;
      try {
        H5::Exception::dontPrint();
        get_hdf5<int>(&MagneticField,  file, (char *)   "/Hamiltonian/MagneticField");
      }
      catch (H5::Exception& e){}
	        
      file->close();
    }


    
    // Set the vector potential
    vect_pot.setZero();
		vect_pot(0,1) = MagneticField*1.0/Lt[0]*2.0*M_PI;
    
    

    Nd = 1;
    N = 1;
    Nt = 1;
    NStr = 1;
    n_threads = 1;
    
    for(unsigned i = 0; i < D; i++)
      {
	ld[i] = Lt[i]/nd[i];
	Ld[i] = ld[i] + 2*NGHOSTS;
	lStr[i] = ld[i]/STRIDE;
	Nd *= Ld[i];
	N  *= ld[i];
	Nt *= Lt[i] ;
	NStr *= lStr[i] ;
	n_threads *= nd[i];
      }
    
    Lt[D] = Orb;
    Ld[D] = Orb;
    ld[D] = Orb;
    lStr[D] = Orb;
    nd[D] = 1;
    
    Size = N * Orb;
    Sized = Nd * Orb;
    Sizet = Nt * Orb;
    SizetVacancies = 0;
    thread_id = omp_get_thread_num();
  };
  
  unsigned get_BorderSize() {
    unsigned size;
    switch (D) {
    case 1 :
      size = 2 * Orb * n_threads * NGHOSTS;
      break;
    case 2:
      size = 2 * std::max(Ld[0],ld[1]) * Orb * n_threads * NGHOSTS;
      break;
    case 3:
      size = (2 * std::max(Ld[0] * Ld[1], std::max(Ld[0] * ld[2] , ld[1]*ld[2]) ) * Orb * n_threads) * NGHOSTS;
      break;
    default:
      std::cout << "Error in LatticeBuilding.hpp. Exiting.\n";
      exit(1);
    }
    return size;
  };
  
  template <typename T1>
  void  convertCoordinates(Coordinates<T1, D + 1> & dest, Coordinates<T1, D + 1> & source)
  {
    /*
     * Convert between the types basis defined in the LatticeStructure
     */
    Coordinates<T1, D + 1> xd(T1(thread_id), nd);
    // Convert from Ld to Lt
    if( std::equal(std::begin(source.L), std::end(source.L), std::begin(Ld)) && std::equal(std::begin(dest.L), std::end(dest.L), std::begin(Lt)))
      {
	for(unsigned i = 0; i < D; i++)
	  dest.coord[i] =  (source.coord[i] + xd.coord[i] * ld[i] - NGHOSTS + Lt[i])%Lt[i] ;
	dest.coord[D] = source.coord[D];
	dest.set_index(dest.coord);
      }
    // Convert from ld to Lt
    if( std::equal(std::begin(source.L), std::end(source.L), std::begin(ld)) && std::equal(std::begin(dest.L), std::end(dest.L), std::begin(Lt)))
      {
	for(unsigned i = 0; i < D; i++)
	  dest.coord[i] =  source.coord[i] + xd.coord[i] * ld[i];
	dest.coord[D] = source.coord[D];
	dest.set_index(dest.coord);
      }
    
    // Convert from ld to Ld
    if( std::equal(std::begin(source.L), std::end(source.L), std::begin(ld)) && std::equal(std::begin(dest.L), std::end(dest.L), std::begin(Ld)))
      {
	for(unsigned i = 0; i < D; i++)
	  dest.coord[i] =  source.coord[i] + NGHOSTS ;
	dest.coord[D] = source.coord[D];
	dest.set_index(dest.coord);
      }
    
    
    // Convert from Lt to ld
    if( std::equal(std::begin(source.L), std::end(source.L), std::begin(Lt)) && std::equal(std::begin(dest.L), std::end(dest.L), std::begin(ld)))
      {
	for(unsigned i = 0; i < D; i++)
	  dest.coord[i] =  source.coord[i] - xd.coord[i] * ld[i];
	dest.coord[D] = source.coord[D];
	dest.set_index(dest.coord);
      }

    // Convert from Lt to Ld
    if( std::equal(std::begin(source.L), std::end(source.L), std::begin(Lt)) && std::equal(std::begin(dest.L), std::end(dest.L), std::begin(Ld)))
      {
	for(unsigned i = 0; i < D; i++)
	  dest.coord[i] =  source.coord[i] - xd.coord[i] * ld[i] + NGHOSTS;
	dest.coord[D] = source.coord[D];
	dest.set_index(dest.coord);
      }
    
    // Convert from Ld to LStr
    if( std::equal(std::begin(source.L), std::end(source.L), std::begin(Ld)) && std::equal(std::begin(dest.L), std::end(dest.L), std::begin(lStr)))
      {
	for(unsigned i = 0; i < D; i++)
	  dest.coord[i] =  (source.coord[i] -NGHOSTS)/STRIDE;
	dest.coord[D] = 0;
	dest.set_index(dest.coord);
      }


    // Convert from ld to LStr
    if( std::equal(std::begin(source.L), std::end(source.L), std::begin(ld)) && std::equal(std::begin(dest.L), std::end(dest.L), std::begin(lStr)))
      {
	for(unsigned i = 0; i < D; i++)
	  dest.coord[i] =  source.coord[i]/STRIDE; 
	dest.coord[D] = 0;
	dest.set_index(dest.coord);
      }
    
  };
  
  unsigned domain_number (long index) {
    Coordinates<long, D + 1> LATT(Lt);
    Coordinates<long, D + 1> n(nd);
    LATT.set_coord(index);
    for (unsigned i = 0; i < D; i++)
      LATT.coord[i] /= ld[i];
    LATT.coord[D] = 0;
    return unsigned(n.set_index(LATT.coord).index);
  }

  void print_coordinates(std::size_t pos1, std::size_t pos2)
  {
    Coordinates<long, D + 1> Latt1(Ld);
    Coordinates<long, D + 1> Latt2(Ld);
    Eigen::Matrix<double, D, 1> r1, r2;
    Eigen::Map<Eigen::Matrix<std::ptrdiff_t, D, 1>> v1(Latt1.coord), v2(Latt2.coord); // Column vectors
    Latt1.set_coord(pos1);
    Latt2.set_coord(pos2);
    
    r1 = rOrb.col(Latt1.coord[D]) + rLat * v1.template cast<double>();
    r2 = rOrb.col(Latt2.coord[D]) + rLat * v2.template cast<double>();
    //std::cout << r1.transpose() << std::endl;
    //std::cout << r2.transpose() << std::endl << std::endl;
  }

  bool test_ghosts(  Coordinates<std::size_t, D + 1> & Latt)
  {
    // This function tests if the coordinates are in the ghosts
    // 0 is in the ghosts
    // 1 isn't in the ghosts
    
    bool teste = 1;
  
    for(int j = 0; j < int(D); j++)
      if(teste && (Latt.coord[j] < NGHOSTS || Latt.coord[j] >= std::ptrdiff_t(Ld[j] - NGHOSTS)) )
	teste = 0;                                        // node is in the ghosts!
      else  if(Latt.coord[j] < 0 || Latt.coord[j] > std::ptrdiff_t(Ld[j] - 1))
	{
	  //std::cout << "Big Mistake" << std::endl;
	  //std::cout.flush();
	}
    return teste;
  }
    

  
  
};

