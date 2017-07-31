
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
  unsigned Sized; // Size of vector for each subdomain (with ghosts)
  unsigned Size; // Size of vector for each subdomain (without ghosts)
  unsigned long Sizet; // Size of full Hilbert Space

  Eigen::Matrix<double, D, D> rLat;  // The vectors are organized by columns 
  Eigen::MatrixXd rOrb;
  
  unsigned nd[D + 1]; // Number of domains in each dimension (the last dimension corresponding with Orbitals are not decomposed)
  unsigned n_threads; // Number of threads
  
  unsigned Lt[D+1]; // Dimensions of the global sample
  unsigned Ld[D+1]; // Dimensions of each sub-domain (domain  + ghosts) 
  unsigned ld[D+1]; // Dimensions of each sub-domain (domain) 
  unsigned Bd[D+1]; // Information about periodic or non-periodic boundary conditions

  unsigned long Nt; // Number of lattice postions of the global sample
  unsigned Nd; // Number of lattice postions of the sub-domain with ghosts
  unsigned N; // Number of lattice postions of the sub-domain without ghosts
  
  unsigned Orb; // Number of orbitals
  unsigned thread_id; // thread identification

  
  LatticeStructure(char *name )
  {
    
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      get_hdf5<unsigned>(&Orb, file, (char *) "/NOrbitals");
      get_hdf5<double>(rLat.data(), file, (char *) "/LattVectors");    
      rOrb = Eigen::MatrixXd::Zero(Orb,D);
      get_hdf5<double>(rOrb.data(), file, (char *) "/OrbPositions");
      
      get_hdf5<unsigned>(Lt, file, (char *) "/L");
      get_hdf5<unsigned>(Bd, file, (char *) "/Boundaries");
      get_hdf5<unsigned>(nd, file, (char *) "/Divisions");
      file->close();
    }
    Nd = 1;
    N = 1;
    Nt = 1;
    n_threads = 1;

    for(unsigned i = 0; i < D; i++)
      {
	
	ld[i] = Lt[i]/nd[i];
	Ld[i] = ld[i] + 2;
	Nd *= Ld[i];
	N  *= ld[i];
	Nt *= ((unsigned long) Lt[i] );
	n_threads *= nd[i];
      }
    
    Lt[D] = Orb;
    Ld[D] = Orb;
    ld[D] = Orb;
    nd[D] = 1;
    
    Size = N * Orb;
    Sized = Nd * Orb;
    Sizet = Nt * Orb;
    thread_id = omp_get_thread_num();
  };

  unsigned get_BorderSize() {
    unsigned size;
    switch (D) {
    case 1 :
      size = 2 * Orb * n_threads;
      break;
    case 2:
      size = 2 * std::max(Ld[0],ld[1]) * Orb * n_threads;
      break;
    case 3:
      size = (2 * std::max(Ld[0] * Ld[1], std::max(Ld[0] * ld[2] , ld[1]*ld[2]) ) * Orb * n_threads);
      break;
    default:
      exit(0);
    }
    return size;
  };

  template <typename T1>
  void  buildGlobalCoordinates(Coordinates<T1, D + 1> & global) {
    // This function receive a Coordinate in some basis  
    Coordinates<T1, D + 1> xd(nd);
    xd.set_coord(T1(thread_id));
    for(unsigned i = 0; i < D; i++)
      global.coord[i] =  (global.L[i] + global.coord[i] + xd.coord[i] * ld[i] - 1  ) % global.L[i];
    global.set_index(global.coord);
  }
};

