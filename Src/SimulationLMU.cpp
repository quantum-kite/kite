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
  }

template <typename T,unsigned D>
void Simulation<T,D>::LMU(int NDisorder, int NMoments, Eigen::Array<unsigned long, -1, 1> positions){
    debug_message("Entered Simulation::MU\n");

    typedef typename extract_value_type<T>::value_type value_type;
    Eigen::Matrix<T, 1, 2> tmp;
    int NPositions = positions.size();
    unsigned long pos;

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
            kpm0.build_site(pos);
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
    debug_message("Left Simulation::MU\n");
}

template <typename T,unsigned D>
void Simulation<T,D>::calc_LDOS(){
    debug_message("Entered Simulation::calc_LDOS\n");

    // Make sure that all the threads are ready before opening any files
    // Some threads could still be inside the Simulation constructor
    // This barrier is essential
#pragma omp barrier

    //Check if the local density of states needs to be calculated
  bool local_calculate_ldos = false;
#pragma omp master
  {
        Global.calculate_ldos = false;
        H5::H5File * file = new H5::H5File(name, H5F_ACC_RDONLY);
        try{
          int dummy_var;
          get_hdf5<int>(&dummy_var, file, (char *) "/Calculation/ldos/NumDisorder");
          Global.calculate_ldos = true;
        } catch(H5::Exception& e) {
          debug_message("ldos: no need to calculate.\n");
        }
          file->close();  
          delete file;
  }
#pragma omp barrier
        
        // Now calculate it
          unsigned ldos_NumMoments;
          unsigned ldos_NumDisorder;
          Eigen::Array<unsigned long, -1, 1> ldos_Orbitals;
          Eigen::Array<unsigned long, -1, 1> ldos_Positions;

        local_calculate_ldos = Global.calculate_ldos;
        if(local_calculate_ldos){
#pragma omp master
      {
        std::cout << "Calculating LDoS.\n";
      }
#pragma omp barrier

#pragma omp critical
  {

          H5::DataSet * dataset;
          H5::DataSpace * dataspace;
          hsize_t dim[1];
          H5::H5File * file = new H5::H5File(name, H5F_ACC_RDONLY);
          dataset            = new H5::DataSet(file->openDataSet("/Calculation/ldos/Orbitals")  );
          dataspace          = new H5::DataSpace(dataset->getSpace());
          dataspace -> getSimpleExtentDims(dim, NULL);
          dataspace->close(); delete dataspace;
          dataset->close();   delete dataset;

          ldos_Orbitals = Eigen::Array<unsigned long, -1, 1>::Zero(dim[0],1);
          ldos_Positions = Eigen::Array<unsigned long, -1, 1>::Zero(dim[0],1);

          get_hdf5<unsigned>(&ldos_NumMoments, file, (char *) "/Calculation/ldos/NumMoments");
          get_hdf5<unsigned>(&ldos_NumDisorder, file, (char *) "/Calculation/ldos/NumDisorder");
          get_hdf5<unsigned long>(ldos_Orbitals.data(), file, (char *) "/Calculation/ldos/Orbitals");
          get_hdf5<unsigned long>(ldos_Positions.data(), file, (char *) "/Calculation/ldos/FixPosition");
          file->close();  
          delete file;
  }
#pragma omp barrier



          
          Eigen::Array<unsigned long, -1, 1> total_positions;
          total_positions = ldos_Positions + ldos_Orbitals*r.Lt[0]*r.Lt[1];
          LMU(ldos_NumDisorder, ldos_NumMoments, total_positions);
        }
    

    debug_message("Left Simulation::calc_LDOS\n");
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
