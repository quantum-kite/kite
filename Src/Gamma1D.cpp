/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2020, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/


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

void Simulation<T,D>::Gamma1D(int NRandomV, int NDisorder, int N_moments,
    std::vector<std::vector<unsigned>> indices, const std::string & name_dataset){

  int num_velocities = 0;
  for(int i = 0; i < int(indices.size()); i++)
    num_velocities += indices.at(i).size();
  int factor = 1 - (num_velocities % 2)*2;
    
  // Initialize the KPM vectors that will be needed to run the 1D Gamma matrix
  KPM_Vector<T,D> kpm0(1, *this);
  KPM_Vector<T,D> kpm1(2, *this);
		
  // Make sure the local gamma matrix is zeroed
  Eigen::Array<T, -1, -1> gamma = Eigen::Array<T, -1, -1 >::Zero(1, N_moments);
  Eigen::Matrix<T, 1, 2> tmp =  Eigen::Matrix < T, 1, 2> ::Zero();		

  long average = 0;
  for(int disorder = 0; disorder < NDisorder; disorder++){
    h.generate_disorder();

    for(unsigned it = 0; it < indices.size(); it++){
      h.build_velocity(indices.at(it), it);
    }

    for(int randV = 0; randV < NRandomV; randV++){
        
      kpm0.initiate_vector();			// original random vector
      kpm1.set_index(0);
      kpm1.v.col(0) = kpm0.v.col(0);
      kpm1.Exchange_Boundaries();

      if(indices.size() != 0)
        generalized_velocity(&kpm1, &kpm0, indices, 0);

      kpm0.v.col(0) = factor*kpm0.v.col(0); // This factor is due to the fact that this Velocity operator is not self-adjoint
      kpm0.empty_ghosts(0);

			
      kpm1.template Multiply<0>();		
      tmp.setZero();
      for(std::size_t ii = 0; ii < r.Sized ; ii += r.Ld[0])
        tmp += kpm0.v.block(ii,0, r.Ld[0], 1).adjoint() * kpm1.v.block(ii, 0, r.Ld[0], 2);

      gamma.matrix().block(0,0,1,2) += (tmp - gamma.matrix().block(0,0,1,2))/value_type(average + 1);			
	
      for(int m = 2; m < N_moments; m += 2){
        kpm1.template Multiply<1>();
        kpm1.template Multiply<1>();
        tmp.setZero();
        for(std::size_t ii = 0; ii < r.Sized ; ii += r.Ld[0])
          tmp += kpm0.v.block(ii,0, r.Ld[0], 1).adjoint() * kpm1.v.block(ii, 0, r.Ld[0], 2);

        gamma.matrix().block(0, m,1,2) += (tmp - gamma.matrix().block(0,m,1,2))/value_type(average + 1);

      }
  //std::cout << "got to line " << __LINE__ << " in file " << __FILE__ << "\n" << std::flush;

      average++;
    }
  } 

  store_gamma1D(&gamma, name_dataset);
}



template <typename T,unsigned D>
void Simulation<T,D>::store_gamma1D(Eigen::Array<T, -1, -1> *gamma, 
                                  const std::string & name_dataset){
  debug_message("Entered store_gamma\n");
  // The whole purpose of this function is to take the Gamma matrix calculated by



  long int size_gamma = gamma->cols();
#pragma omp master
  Global.general_gamma = Eigen::Array<T, -1, -1 > :: Zero(1, size_gamma);
#pragma omp barrier
#pragma omp critical
  Global.general_gamma += *gamma;
#pragma omp barrier

    
#pragma omp master
  {
    H5::H5File * file = new H5::H5File(name, H5F_ACC_RDWR);
    write_hdf5(Global.general_gamma, file, name_dataset);
    delete file;
  }
#pragma omp barrier    

    
  debug_message("Left store_gamma\n");
}


template void Simulation<float ,1u>::Gamma1D(int, int, int, std::vector<std::vector<unsigned>>, const std::string &);
template void Simulation<double ,1u>::Gamma1D(int, int, int, std::vector<std::vector<unsigned>>, const std::string &);
template void Simulation<long double ,1u>::Gamma1D(int, int, int, std::vector<std::vector<unsigned>>, const std::string &);
template void Simulation<std::complex<float> ,1u>::Gamma1D(int, int, int, std::vector<std::vector<unsigned>>, const std::string &);
template void Simulation<std::complex<double> ,1u>::Gamma1D(int, int, int, std::vector<std::vector<unsigned>>, const std::string &);
template void Simulation<std::complex<long double> ,1u>::Gamma1D(int, int, int, std::vector<std::vector<unsigned>>, const std::string &);

template void Simulation<float ,3u>::Gamma1D(int, int, int, std::vector<std::vector<unsigned>>, const std::string &);
template void Simulation<double ,3u>::Gamma1D(int, int, int, std::vector<std::vector<unsigned>>, const std::string &);
template void Simulation<long double ,3u>::Gamma1D(int, int, int, std::vector<std::vector<unsigned>>, const std::string &);
template void Simulation<std::complex<float> ,3u>::Gamma1D(int, int, int, std::vector<std::vector<unsigned>>, const std::string &);
template void Simulation<std::complex<double> ,3u>::Gamma1D(int, int, int, std::vector<std::vector<unsigned>>, const std::string &);
template void Simulation<std::complex<long double> ,3u>::Gamma1D(int, int, int, std::vector<std::vector<unsigned>>, const std::string &);

template void Simulation<float ,2u>::Gamma1D(int, int, int, std::vector<std::vector<unsigned>>, const std::string &);
template void Simulation<double ,2u>::Gamma1D(int, int, int, std::vector<std::vector<unsigned>>, const std::string &);
template void Simulation<long double ,2u>::Gamma1D(int, int, int, std::vector<std::vector<unsigned>>, const std::string &);
template void Simulation<std::complex<float> ,2u>::Gamma1D(int, int, int, std::vector<std::vector<unsigned>>, const std::string &);
template void Simulation<std::complex<double> ,2u>::Gamma1D(int, int, int, std::vector<std::vector<unsigned>>, const std::string &);
template void Simulation<std::complex<long double> ,2u>::Gamma1D(int, int, int, std::vector<std::vector<unsigned>>, const std::string &);


template void Simulation<float ,1u>::store_gamma1D(Eigen::Array<float, -1, -1>* , const std::string &);
template void Simulation<double ,1u>::store_gamma1D(Eigen::Array<double, -1, -1>* , const std::string &);
template void Simulation<long double ,1u>::store_gamma1D(Eigen::Array<long double, -1, -1>* , const std::string &);
template void Simulation<std::complex<float> ,1u>::store_gamma1D(Eigen::Array<std::complex<float>, -1, -1>* , const std::string &);
template void Simulation<std::complex<double> ,1u>::store_gamma1D(Eigen::Array<std::complex<double>, -1, -1>* , const std::string &);
template void Simulation<std::complex<long double> ,1u>::store_gamma1D(Eigen::Array<std::complex<long double>, -1, -1>* , const std::string &);

template void Simulation<float ,2u>::store_gamma1D(Eigen::Array<float, -1, -1>* , const std::string &);
template void Simulation<double ,2u>::store_gamma1D(Eigen::Array<double, -1, -1>* , const std::string &);
template void Simulation<long double ,2u>::store_gamma1D(Eigen::Array<long double, -1, -1>* , const std::string &);
template void Simulation<std::complex<float> ,2u>::store_gamma1D(Eigen::Array<std::complex<float>, -1, -1>* , const std::string &);
template void Simulation<std::complex<double> ,2u>::store_gamma1D(Eigen::Array<std::complex<double>, -1, -1>* , const std::string &);
template void Simulation<std::complex<long double> ,2u>::store_gamma1D(Eigen::Array<std::complex<long double>, -1, -1>* , const std::string &);

template void Simulation<float ,3u>::store_gamma1D(Eigen::Array<float, -1, -1>* , const std::string &);
template void Simulation<double ,3u>::store_gamma1D(Eigen::Array<double, -1, -1>* , const std::string &);
template void Simulation<long double ,3u>::store_gamma1D(Eigen::Array<long double, -1, -1>* , const std::string &);
template void Simulation<std::complex<float> ,3u>::store_gamma1D(Eigen::Array<std::complex<float>, -1, -1>* , const std::string &);
template void Simulation<std::complex<double> ,3u>::store_gamma1D(Eigen::Array<std::complex<double>, -1, -1>* , const std::string &);
template void Simulation<std::complex<long double> ,3u>::store_gamma1D(Eigen::Array<std::complex<long double>, -1, -1>* , const std::string &);




