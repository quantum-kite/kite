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
void Simulation<T,D>::store_LMU(Eigen::Array<T, -1, -1> *gamma){
    debug_message("Entered store_lmu\n");

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
    write_hdf5(Global.general_gamma, file, "/Calculation/ldos/lMU");
    file->close();
    delete file;
}
#pragma omp barrier    
    debug_message("Left store_lmu\n");
  };

template <typename T,unsigned D>
void Simulation<T,D>::LMU(int NDisorder, int NMoments, Eigen::Array<unsigned long, -1, 1> positions){
    typedef typename extract_value_type<T>::value_type value_type;
    Eigen::Matrix<T, 1, 2> tmp;
    int NPositions = positions.size();
    unsigned long pos;

    unsigned long T_thread[D+1]; // index of the thread
    unsigned long x_thread[D+1]; // position within the thread

    Coordinates<unsigned long, D+1> thread_coords(r.ld);
    Coordinates<unsigned long, D+1> thread_coords_gh(r.Ld);
    Coordinates<unsigned long, D+1> thread(r.nd);
    Coordinates<unsigned long, D+1> total_coords(r.Lt);
    bool correct_thread;

    
    KPM_Vector<T,D> kpm0(1, *this);      // initial random vector
    KPM_Vector<T,D> kpm1(2, *this); // left vector that will be Chebyshev-iterated on

    // initialize the local gamma matrix and set it to 0
    Eigen::Array<T, -1, -1> gamma = Eigen::Array<T, -1, -1 >::Zero(NMoments, NPositions);

    // start the kpm iteration
    Eigen::Array<long, -1, 1> average;
    average = Eigen::Array<long, -1, 1>::Zero(NPositions,1);

    for(int disorder = 0; disorder < NDisorder; disorder++){
        h.generate_disorder();

        for(int pos_index = 0; pos_index < NPositions; pos_index++){

            // Iterate over the list of positions 
            // each position is a global position. We need to find which thread contains
            // this position and the coordinates in that thread. Only THAT thread should
            // have a starting vector different from zero
            pos = positions(pos_index);
#pragma omp critical
{
            //std::cout << "Thread id: " << r.thread_id << "\n";
            total_coords.set_coord(pos);
            //std::cout << "position: " << pos << "\n" << std::flush;
            for(unsigned d = 0; d < D; d++){
                //std::cout << "total_coord: " << total_coords.coord[d] << "\n";
                T_thread[d] = total_coords.coord[d]/r.ld[d];
                x_thread[d] = total_coords.coord[d]%r.ld[d];
                //std::cout << "index: " << d << "\n";
                //std::cout << "   thread in direction: " << T_thread[d] << "\n";
                //std::cout << "   x in thread: " << x_thread[d] << "\n";
            }
            T_thread[D] = 0;
            x_thread[D] = total_coords.coord[D];
            thread_coords.set_index(x_thread);
            thread.set_index(T_thread);

            // convert to coordinates with ghosts
            r.convertCoordinates(thread_coords_gh, thread_coords); 
            //for(unsigned d = 0; d < D+1; d++){
                //std::cout << thread_coords_gh.coord[d] << "\n"; 
            //}
            
            // check if the site is in the current thread
            correct_thread = thread.index == r.thread_id;


            //std::cout << "thread_id: " << r.thread_id << "\n";
            //std::cout << "index: " << thread.index << "\n";
            //std::cout << "Is correct thread? " << T(correct_thread) << "\n";

}
#pragma omp barrier
            kpm0.v.setZero();
            kpm0.v(thread_coords_gh.index,0) = T(correct_thread);
            kpm0.set_index(0);
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

                gamma(n, pos_index) += (tmp(0,0) - gamma(n, pos_index))/value_type(average(pos_index) + 1);			
                gamma(n+1, pos_index) += (tmp(0,1) - gamma(n+1, pos_index))/value_type(average(pos_index) + 1);			
            }
            average(pos_index)++;
        } 
    }
    store_LMU(&gamma);
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
