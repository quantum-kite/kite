#include "Generic.hpp"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"

template <unsigned D>
LatticeStructure<D>::LatticeStructure(char *name )
{
  
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
    

    try {
      H5::Exception::dontPrint();
      get_hdf5<int>(&MagneticField, file, (char *) "/Hamiltonian/MagneticFieldMul");
    }
    catch (H5::Exception& e){}
    file->close();
  }

  // Set the ghost_correlation potential and normalize it to the size of the system
  ghost_pot.setZero();
  ghost_pot(0,1) = MagneticField * 1.0 /Lt[1]*2.0*M_PI;

  test_divisibility();
    
  Nd = 1;
  N = 1;
  Nt = 1;
  NStr = 1;
  n_threads = 1;
    
  for(unsigned i = 0; i < D; i++)
    {
      ld[i] = Lt[i]/nd[i];
      Ld[i] = ld[i] + 2*NGHOSTS;
      lStr[i] = ld[i] / TILE;
      Nd *= Ld[i];
      N  *= ld[i];
      Nt *= Lt[i] ;
      NStr *= lStr[i] ;
      n_threads *= nd[i];
    }

  std::fill_n(lB3, D, 3); 
  lB3[D]  = Orb;  
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
    
  Coordinates<unsigned, D + 1> dist(nd);
  dist.set_coord(unsigned(thread_id));
  
  // Test if subdomain is in the Global border and if it has open boundaries set to FALSE
  for(unsigned i = 0; i < D; i++)
    {
      boundary[i][0] = (dist.coord[i] == 0         && Bd[i] == 0 ? false : true); 
      boundary[i][1] = (dist.coord[i] == nd[i] - 1 && Bd[i] == 0 ? false : true);	
    }
    
}

template <unsigned D>
unsigned LatticeStructure<D>::get_BorderSize() {
  unsigned size;
  switch (D) {
  case 1 :
    size = 2 * Orb * n_threads * NGHOSTS;
    break;
  case 2:
    size = 2 * std::max(Ld[0],ld[1]) * Orb * n_threads * NGHOSTS;
    break;
  case 3:
    size = 2 * std::max(Ld[0] * Ld[1], std::max(Ld[0] * ld[2] , ld[1]*ld[2]) ) * Orb * n_threads * NGHOSTS;
    break;
  default:
    std::cout << "Error in LatticeBuilding.hpp. Exiting.\n";
    exit(1);
  }
  return size;
}


template <unsigned D>
void LatticeStructure<D>::test_divisibility() {
  debug_message("Entered LatticeStructure::test_divisibility.\n");
  // Test if TILE x nd divides the length

  for(unsigned i = 0; i < D; i++){
    if(Lt[i]%(nd[i] * TILE) != 0){
      std::cout << "The system size in direction " << i << " (" << Lt[i] <<  ") ";
      std::cout << "must be a multiple of the number of divisions in that ";
      std::cout << "direction (" << nd[i] << ") times TILE (" << TILE << "). ";
      std::cout << "Exiting.\n";
      exit(1);
    }
  } 

  debug_message("Left LatticeStructure::test_divisibility.\n");
}

template <unsigned D>
template <typename T1>
void  LatticeStructure<D>::convertCoordinates(Coordinates<T1, D + 1> & dest, Coordinates<T1, D + 1> & source)
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
        dest.coord[i] = (source.coord[i] -NGHOSTS) / TILE;
      dest.coord[D] = 0;
      dest.set_index(dest.coord);
    }


  // Convert from ld to LStr
  if( std::equal(std::begin(source.L), std::end(source.L), std::begin(ld)) && std::equal(std::begin(dest.L), std::end(dest.L), std::begin(lStr)))
    {
      for(unsigned i = 0; i < D; i++)
        dest.coord[i] = source.coord[i] / TILE; 
      dest.coord[D] = 0;
      dest.set_index(dest.coord);
    }

  // Convert from Ld to nd
  if( std::equal(std::begin(source.L), std::end(source.L), std::begin(Lt)) && std::equal(std::begin(dest.L), std::end(dest.L), std::begin(nd)))
    {
      for(unsigned i = 0; i < D; i++)
        dest.coord[i] =  source.coord[i]/ld[i]; 
      dest.coord[D] = 0;
      dest.set_index(dest.coord);
    }

  
}

template <unsigned D>
unsigned LatticeStructure<D>::domain_number (long index) {
  Coordinates<long, D + 1> LATT(Lt);
  Coordinates<long, D + 1> n(nd);
  LATT.set_coord(index);
  for (unsigned i = 0; i < D; i++)
    LATT.coord[i] /= ld[i];
  LATT.coord[D] = 0;
  return unsigned(n.set_index(LATT.coord).index);
}

template <unsigned D>
void LatticeStructure<D>::print_coordinates(std::size_t pos1, std::size_t pos2)
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

template <unsigned D>
bool LatticeStructure<D>::test_ghosts(  Coordinates<std::size_t, D + 1> & Latt)
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

template struct LatticeStructure<1u>;
template struct LatticeStructure<2u>;
template struct LatticeStructure<3u>;
template void LatticeStructure<1u>::convertCoordinates<std::size_t>(Coordinates<std::size_t, 2u> &, Coordinates<std::size_t, 2u> &);
template void LatticeStructure<2u>::convertCoordinates<std::size_t>(Coordinates<std::size_t, 3u> &, Coordinates<std::size_t, 3u> &);
template void LatticeStructure<3u>::convertCoordinates<std::size_t>(Coordinates<std::size_t, 4u> &, Coordinates<std::size_t, 4u> &);
template void LatticeStructure<1u>::convertCoordinates<long>(Coordinates<long, 2u> &, Coordinates<long, 2u> &);
template void LatticeStructure<2u>::convertCoordinates<long>(Coordinates<long, 3u> &, Coordinates<long, 3u> &);
template void LatticeStructure<3u>::convertCoordinates<long>(Coordinates<long, 4u> &, Coordinates<long, 4u> &);

