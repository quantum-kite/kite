
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
  unsigned Sized;
  unsigned Size;

  double rLat[D][D];   // first index is the vector the second the component in cartesian coordinates 
  double (*rOrb)[D];
  unsigned nd[D + 1];
  unsigned n_threads;
  
  unsigned Lt[D+1];
  unsigned Ld[D+1];
  unsigned ld[D+1];
  unsigned Bd[D+1];

  unsigned long Nt;
  unsigned Nd;
  unsigned N;
  
  unsigned Orb;
  unsigned thread_id;

  
  LatticeStructure(char *name )
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    get_hdf5<unsigned>(&Orb, file, (char *) "/NOrbitals");
    get_hdf5<double>(&rLat[0][0], file, (char *) "/LattVector");
    
    rOrb =  (double(*)[D]) malloc(Orb * D * sizeof(double) );
    get_hdf5<double>(&rOrb[0][0], file, (char *) "/OrbPositions");
    
    get_hdf5<unsigned>(Lt, file, (char *) "/L");
    get_hdf5<unsigned>(Bd, file, (char *) "/Boundaries");
    get_hdf5<unsigned>(nd, file, (char *) "/Divisions");
    file->close();
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
    Coordinates<T1, D + 1> xd(nd);
    xd.set_coord(T1(thread_id));
    for(unsigned i = 0; i < D; i++)
      global.coord[i] =  (global.L[i] + global.coord[i] + xd.coord[i] * ld[i] - 1  ) % global.L[i];
    global.set_index(global.coord);
  }
};

