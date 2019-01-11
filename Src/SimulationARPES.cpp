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
#include "queue.hpp"
#include "Simulation.hpp"
#include "Hamiltonian.hpp"
#include "KPM_VectorBasis.hpp"
#include "KPM_Vector.hpp"

template <typename T,unsigned D>
void Simulation<T,D>::store_ARPES(Eigen::Array<T, -1, -1> *gamma){
    debug_message("Entered store_ARPES\n");

    long int nMoments   = gamma->rows();
    long int nPositions = gamma->cols();

#pragma omp master
	Global.general_gamma = Eigen::Array<T, -1, -1 > :: Zero(nMoments, nPositions);
#pragma omp barrier
#pragma omp critical
	Global.general_gamma += *gamma;
#pragma omp barrier
    
    
#pragma omp master
{
    //std::cout << "Printing huge matrix, brb\n";
    //std::cout << Global.general_gamma << "\n";
    H5::H5File * file = new H5::H5File(name, H5F_ACC_RDWR);
    write_hdf5(Global.general_gamma, file, "/Calculation/ARPES/kMU");
    file->close();
    delete file;
}
#pragma omp barrier    
    debug_message("Left store_lmu\n");
  };

template <typename T,unsigned D>
void Simulation<T,D>::ARPES(int NDisorder, int NMoments, Eigen::Array<T, -1, -1> k_vectors){
    typedef typename extract_value_type<T>::value_type value_type;

    Eigen::Matrix<T, 1, 2> tmp;
    int Nk_vectors = k_vectors.size();
    Eigen::Matrix<T, -1, 1> k;

    KPM_Vector<T,D> kpm0(1, *this); // initial random vector
    KPM_Vector<T,D> kpm1(2, *this); // left vector that will be Chebyshev-iterated on

    // initialize the local gamma matrix and set it to 0
    Eigen::Array<T, -1, -1> gamma = Eigen::Array<T, -1, -1 >::Zero(NMoments, Nk_vectors);

    // average for each k value
    Eigen::Array<long, -1, 1> average;
    average = Eigen::Array<long, -1, 1>::Zero(Nk_vectors,1);

    // start the kpm iteration
    for(int disorder = 0; disorder < NDisorder; disorder++){
        h.generate_disorder();

        for(int k_index = 0; k_index < Nk_vectors; k_index++){

            // Iterate over the list of k vectors
            k = k_vectors(k_index);

            kpm0.v.setZero();
            kpm0.build_wavepacket(k, weight); // already sets index=0
            kpm0.Exchange_Boundaries();

            kpm1.set_index(0);
            kpm1.v.col(0) = kpm0.v.col(0);
            kpm0.empty_ghosts(0);

            for(int n = 0; n < NMoments; n+=2){
                for(int i = n; i < n+2; i++)
                    if(i!=0) cheb_iteration(&kpm1, i-1);

              
                tmp.setZero();
                for(std::size_t ii = 0; ii < r.Sized ; ii += r.Ld[0])
                    tmp += kpm0.v.block(ii,0, r.Ld[0], 1).adjoint() * kpm1.v.block(ii, 0, r.Ld[0], 2);

                gamma(n, k_index) += (tmp(0,0) - gamma(n, k_index))/value_type(average(k_index) + 1);			
                gamma(n+1, k_index) += (tmp(0,1) - gamma(n+1, k_index))/value_type(average(k_index) + 1);			
            }
            average(k_index)++;
        } 
    }
    store_ARPES(&gamma);
}

template <typename T, unsigned DIM>
void Simulation<T, DIM>::calc_ARPES(){
    // Checks if ARPES needs to be calculated. If it does, it will search
    // the input .h5 file for the needed parameters for that calculation
    // and then performs the calculation


    int NDisorder = 1;
    int NMoments = 128;
    Eigen::k_vectors


    ARPES(NDisorder, NMoments, k_vectors);

}

template class Simulation<float ,1u>;
template class Simulation<double ,1u>;
template class Simulation<long double ,1u>;
template class Simulation<std::complex<float> ,1u>;
template class Simulation<std::complex<double> ,1u>;
template class Simulation<std::complex<long double> ,1u>;

template class Simulation<float ,3u>;
template class Simulation<double ,3u>;
template class Simulation<long double ,3u>;
template class Simulation<std::complex<float> ,3u>;
template class Simulation<std::complex<double> ,3u>;
template class Simulation<std::complex<long double> ,3u>;

template class Simulation<float ,2u>;
template class Simulation<double ,2u>;
template class Simulation<long double ,2u>;
template class Simulation<std::complex<float> ,2u>;
template class Simulation<std::complex<double> ,2u>;
template class Simulation<std::complex<long double> ,2u>;
