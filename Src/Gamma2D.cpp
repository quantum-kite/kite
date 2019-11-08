#include "Generic.hpp"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"
#include "Global.hpp"
#include "Random.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
template <typename T, unsigned D>
class Hamiltonian;
template <typename T, unsigned D>
class KPM_Vector;
#include "Simulation.hpp"
#include "Hamiltonian.hpp"
#include "KPM_VectorBasis.hpp"
#include "KPM_Vector.hpp"


template <typename T,unsigned D>
void Simulation<T,D>::Gamma2D(int NRandomV, int NDisorder, std::vector<int> N_moments, 
                              std::vector<std::vector<unsigned>> indices, std::string name_dataset){
  Eigen::Matrix<T, MEMORY, MEMORY> tmp;
  // This function calculates all kinds of two-dimensional gamma matrices such
  // as Tr[V^a Tn v^b Tm] = G_nm
  //
  // The matrices are stored as
  //
  // | G_00   G_01   G_02   ...   G_0M | 
  // | G_10   G_11   G_12   ...   G_1M | 
  // | G_20   G_21   G_22   ...   G_2M | 
  // | ...    ...    ...    ...   ...  |
  // | G_N0   G_N1   G_N2   ...   G_NM | 
  //
  // For example, a 3x3 matrix would be represented as
  //
  // | G_00   G_01   G_02 | 
  // | G_10   G_11   G_12 | 
  // | G_20   G_21   G_22 | 
  //
  // This function calculates all the kinds of one-dimensional Gamma matrices
  // such as Tr[Tn]    Tr[v^xx Tn]     etc

  typedef typename extract_value_type<T>::value_type value_type;

  int num_velocities = 0;
  for(int i = 0; i < int(indices.size()); i++)
    num_velocities += indices.at(i).size();
  int factor = 1 - (num_velocities % 2)*2;

  //  --------- INITIALIZATIONS --------------
    
  KPM_Vector<T,D> kpm0(1, *this);      // initial random vector
  KPM_Vector<T,D> kpm1(2, *this); // left vector that will be Chebyshev-iterated on
  KPM_Vector<T,D> kpm2(MEMORY, *this); // right vector that will be Chebyshev-iterated on
  KPM_Vector<T,D> kpm3(MEMORY, *this); // kpm1 multiplied by the velocity

  // initialize the local gamma matrix and set it to 0
  int size_gamma = 1;
  for(int i = 0; i < 2; i++){
    if(N_moments.at(i) % 2 != 0){
      std::cout << "The number of moments must be an even number, due to limitations of the program. Aborting\n";
      exit(1);
    }
    size_gamma *= N_moments.at(i);
  }

  Eigen::Array<T, -1, -1> gamma = Eigen::Array<T, -1, -1 >::Zero(1, size_gamma);
 
  // finished initializations


    
    
  // start the kpm iteration
  long average = 0;
  for(int disorder = 0; disorder < NDisorder; disorder++){
    h.generate_disorder();
    for(unsigned it = 0; it < indices.size(); it++)
      h.build_velocity(indices.at(it), it);
    for(int randV = 0; randV < NRandomV; randV++){
        

      kpm0.initiate_vector();			// original random vector. This sets the index to zero
      kpm0.Exchange_Boundaries();
      kpm1.set_index(0);

      generalized_velocity(&kpm1, &kpm0, indices, 0);
        
      // run through the left loop MEMORY iterations at a time
      for(int n = 0; n < N_moments.at(0); n+=MEMORY){
          
        // Iterate MEMORY times. The first time this occurs, we must exclude the zeroth
        // case, because it is already calculated, it's the identity
        for(int i = n; i < n + MEMORY; i++){
          if(i!=0){
            cheb_iteration(&kpm1, i-1);
          }

          kpm3.set_index(i%MEMORY);
          generalized_velocity(&kpm3, &kpm1, indices, 1);
          kpm3.empty_ghosts(i%MEMORY);
        }
          
        // copy the |0> vector to |kpm2>
        kpm2.set_index(0);
        kpm2.v.col(0) = kpm0.v.col(0);
        for(int m = 0; m < N_moments.at(1); m+=MEMORY){

          // iterate MEMORY times, just like before. No need to multiply by v here
          for(int i = m; i < m + MEMORY; i++){
            if(i!=0){
              cheb_iteration(&kpm2, i-1);
            }
          }
          //std::cout << "index2: " << kpm2.get_index() << "\n";
          // Finally, do the matrix product and store the result in the Gamma matrix
          tmp.setZero();
          for(std::size_t ii = 0; ii < r.Sized ; ii += r.Ld[0])
            tmp += kpm3.v.block(ii,0, r.Ld[0], MEMORY).adjoint() * kpm2.v.block(ii, 0, r.Ld[0], MEMORY);
          T flatten;
          long int ind;
          for(int j = 0; j < MEMORY; j++)
            for(int i = 0; i < MEMORY; i++){
              flatten = tmp(i,j);
              ind = (m+j)*N_moments.at(0) + n+i;
              gamma(ind) += (flatten - gamma(ind))/value_type(average + 1);			
            }
        }
      }
      average++;
    }
  } 
  gamma = gamma*factor;
            
  store_gamma(&gamma, N_moments, indices, name_dataset);
}




template <typename T,unsigned D>
void Simulation<T,D>::store_gamma(Eigen::Array<T, -1, -1> *gamma, std::vector<int> N_moments, 
                                  std::vector<std::vector<unsigned>> indices, std::string name_dataset){
  debug_message("Entered store_gamma\n");
  // The whole purpose of this function is to take the Gamma matrix calculated by



  long int size_gamma = gamma->cols();
  int dim = indices.size();

		
  // Number of commutators inside the Gamma matrix. 
  // V^{x}  = [x,H]		-> one commutator
  // V^{xy} = [x,[y,H]]	-> two commutators
  // This is important because the commutator is anti-hermitian. So, an odd number of commutators
  int num_velocities = 0;
  for(int i = 0; i < int(indices.size()); i++)
    num_velocities += indices.at(i).size();
  int factor = 1 - (num_velocities % 2)*2;
		
  switch(dim){
  case 2: {
    Eigen::Array<T,-1,-1> general_gamma = Eigen::Map<Eigen::Array<T,-1,-1>>(gamma->data(), N_moments.at(0), N_moments.at(1));
#pragma omp master
    Global.general_gamma = Eigen::Array<T, -1, -1 > :: Zero(N_moments.at(0), N_moments.at(1));
#pragma omp barrier
#pragma omp critical
    Global.general_gamma.matrix() += (general_gamma.matrix() + factor*general_gamma.matrix().adjoint())/2.0;
#pragma omp barrier
    break;
  }
  case 1: {
    Eigen::Array<T,-1,-1> general_gamma = Eigen::Map<Eigen::Array<T,-1,-1>>(gamma->data(), 1, size_gamma);
#pragma omp master
    Global.general_gamma = Eigen::Array<T, -1, -1 > :: Zero(1, size_gamma);
#pragma omp barrier
#pragma omp critical
    Global.general_gamma += general_gamma;
#pragma omp barrier
    break;
  }
  default:
    std::cout << "You're trying to store a matrix that is not expected by the program. Exiting.\n";
    exit(1);
  }
    

    
#pragma omp master
  {
    H5::H5File * file = new H5::H5File(name, H5F_ACC_RDWR);
    write_hdf5(Global.general_gamma, file, name_dataset);
    delete file;
  }
#pragma omp barrier    

    
  debug_message("Left store_gamma\n");
}


template void Simulation<float ,1u>::Gamma2D(int, int, std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<double ,1u>::Gamma2D(int, int, std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<long double ,1u>::Gamma2D(int, int, std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<std::complex<float> ,1u>::Gamma2D(int, int, std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<std::complex<double> ,1u>::Gamma2D(int, int, std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<std::complex<long double> ,1u>::Gamma2D(int, int, std::vector<int>, std::vector<std::vector<unsigned>>, std::string);

template void Simulation<float ,3u>::Gamma2D(int, int, std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<double ,3u>::Gamma2D(int, int, std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<long double ,3u>::Gamma2D(int, int, std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<std::complex<float> ,3u>::Gamma2D(int, int, std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<std::complex<double> ,3u>::Gamma2D(int, int, std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<std::complex<long double> ,3u>::Gamma2D(int, int, std::vector<int>, std::vector<std::vector<unsigned>>, std::string);

template void Simulation<float ,2u>::Gamma2D(int, int, std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<double ,2u>::Gamma2D(int, int, std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<long double ,2u>::Gamma2D(int, int, std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<std::complex<float> ,2u>::Gamma2D(int, int, std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<std::complex<double> ,2u>::Gamma2D(int, int, std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<std::complex<long double> ,2u>::Gamma2D(int, int, std::vector<int>, std::vector<std::vector<unsigned>>, std::string);



//void Simulation<T,D>::store_gamma3D(Eigen::Array<T, -1, -1> *gamma, std::vector<int> N_moments, 
                                    //std::vector<std::vector<unsigned>> indices, std::string name_dataset){

template void Simulation<float ,1u>::store_gamma(Eigen::Array<float, -1, -1>* , std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<double ,1u>::store_gamma(Eigen::Array<double, -1, -1>* , std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<long double ,1u>::store_gamma(Eigen::Array<long double, -1, -1>* , std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<std::complex<float> ,1u>::store_gamma(Eigen::Array<std::complex<float>, -1, -1>* , std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<std::complex<double> ,1u>::store_gamma(Eigen::Array<std::complex<double>, -1, -1>* , std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<std::complex<long double> ,1u>::store_gamma(Eigen::Array<std::complex<long double>, -1, -1>* , std::vector<int>, std::vector<std::vector<unsigned>>, std::string);

template void Simulation<float ,2u>::store_gamma(Eigen::Array<float, -1, -1>* , std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<double ,2u>::store_gamma(Eigen::Array<double, -1, -1>* , std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<long double ,2u>::store_gamma(Eigen::Array<long double, -1, -1>* , std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<std::complex<float> ,2u>::store_gamma(Eigen::Array<std::complex<float>, -1, -1>* , std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<std::complex<double> ,2u>::store_gamma(Eigen::Array<std::complex<double>, -1, -1>* , std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<std::complex<long double> ,2u>::store_gamma(Eigen::Array<std::complex<long double>, -1, -1>* , std::vector<int>, std::vector<std::vector<unsigned>>, std::string);

template void Simulation<float ,3u>::store_gamma(Eigen::Array<float, -1, -1>* , std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<double ,3u>::store_gamma(Eigen::Array<double, -1, -1>* , std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<long double ,3u>::store_gamma(Eigen::Array<long double, -1, -1>* , std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<std::complex<float> ,3u>::store_gamma(Eigen::Array<std::complex<float>, -1, -1>* , std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<std::complex<double> ,3u>::store_gamma(Eigen::Array<std::complex<double>, -1, -1>* , std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
template void Simulation<std::complex<long double> ,3u>::store_gamma(Eigen::Array<std::complex<long double>, -1, -1>* , std::vector<int>, std::vector<std::vector<unsigned>>, std::string);
