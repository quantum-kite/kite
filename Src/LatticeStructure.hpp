/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/


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
  unsigned lStr[D+1];
  unsigned lB3[D+1];
  
  std::size_t Nt; // Number of lattice postions of the global sample
  std::size_t Nd; // Number of lattice postions of the sub-domain with ghosts
  std::size_t N; // Number of lattice postions of the sub-domain without ghosts
  std::size_t NStr; // Number of lattice postions of the sub-domain without ghosts
  unsigned Orb; // Number of orbitals
  unsigned thread_id; // thread identification
  int MagneticField = 0;
  bool boundary[D][2]; // Information about the Global border in the subdomain 
  Eigen::Matrix<double, D, D> ghost_pot; // ghosts_correlation potential
  
  LatticeStructure(char *);
  unsigned get_BorderSize();
  template <typename T>
  void     convertCoordinates(Coordinates<T, D + 1> & dest, Coordinates<T, D + 1> & source);
  unsigned domain_number (long index);
  void     print_coordinates(std::size_t pos1, std::size_t pos2);
  bool     test_ghosts(  Coordinates<std::size_t, D + 1> & Latt);
  void     test_divisibility();
  
};

