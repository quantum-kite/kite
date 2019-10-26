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

    // Make sure that all the threads are ready before opening any files
    // Some threads could still be inside the Simulation constructor
    // This barrier is essential
#pragma omp barrier

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
    write_hdf5(Global.general_gamma, file, "/Calculation/arpes/kMU");
    file->close();
    delete file;
}
#pragma omp barrier    
    debug_message("Left store_lmu\n");
  }

template <typename T,unsigned D>
void Simulation<T,D>::ARPES(int NDisorder, int NMoments, Eigen::Array<double, -1, -1> & k_vectors, Eigen::Matrix<T, -1, 1> & weight){
    typedef typename extract_value_type<T>::value_type value_type;

    Eigen::Matrix<T, 1, 2> tmp;
    int Nk_vectors = k_vectors.rows();
    Eigen::Matrix<double, -1, 1> k;

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
            k = k_vectors.row(k_index);

            kpm0.v.setZero();
            kpm0.build_planewave(k, weight); // already sets index=0
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



    bool local_calculate_arpes = false;
#pragma omp master
{
    H5::H5File * file = new H5::H5File(name, H5F_ACC_RDONLY);
        Global.calculate_arpes = false;
    try{
        int dummy_var;
        get_hdf5<int>(&dummy_var, file, (char *) "/Calculation/arpes/NumDisorder");
        Global.calculate_arpes = true;
    } catch(H5::Exception& e) {debug_message("ARPES: no need to calculate.\n");}
        file->close();  
        delete file;
}
#pragma omp barrier


      int NumDisorder;
      int NumMoments;
      Eigen::Array<double, -1, -1> k_vectors;
      Eigen::Matrix<T, -1, 1> weight;
    // Fetch the data from the hdf file
    local_calculate_arpes = Global.calculate_arpes;
    if(local_calculate_arpes){
#pragma omp critical
{

      H5::DataSet * dataset;
      H5::DataSpace * dataspace;
      hsize_t dim_k[2], dim_w[2];
      H5::H5File * file  = new H5::H5File(name, H5F_ACC_RDONLY);
      dataset            = new H5::DataSet(file->openDataSet("/Calculation/arpes/k_vector")  );
      dataspace          = new H5::DataSpace(dataset->getSpace());
      dataspace -> getSimpleExtentDims(dim_k, NULL);
      dataspace->close(); delete dataspace;
      dataset->close();   delete dataset;
      
      dataset            = new H5::DataSet(file->openDataSet("/Calculation/arpes/OrbitalWeights")  );
      dataspace          = new H5::DataSpace(dataset->getSpace());
      dataspace -> getSimpleExtentDims(dim_w, NULL);
      dataspace->close(); delete dataspace;
      dataset->close();   delete dataset;

      // Make sure the number of entries in the weight vector is consistent with the
      // number of orbitals
      if(dim_w[1]*dim_w[0] != r.Orb){
        std::cout << "Error in Simulation::calc_ARPES. The number of entries in the orbital "
          "weight vector (" << dim_w[1]*dim_w[0] << ") has to be the same as the number of "
          "orbitals (" << r.Orb << "). Exiting.\n";
        exit(1);
      }

      
      k_vectors = Eigen::Array<double,-1, -1>::Zero(dim_k[1], dim_k[0]);
      weight    = Eigen::Matrix<    T,-1,  1>::Zero(dim_w[1], dim_w[0]);

      //std::cout << "dim_k: " << dim_k[0] << " " << dim_k[1] << "\n";
      //std::cout << "dim_w: " << dim_w[0] << " " << dim_w[1] << "\n";

      // The weights have to be read in doubles before being cast into type T
      Eigen::Matrix<double, -1, 1> weight_test;
      weight_test = Eigen::Matrix<double, -1, 1>::Zero(r.Orb, 1);
      
      get_hdf5    <int>(&NumDisorder,    file, (char *) "/Calculation/arpes/NumDisorder");
      get_hdf5    <int>(&NumMoments,     file, (char *) "/Calculation/arpes/NumMoments" );
      get_hdf5 <double>(weight_test.data(),   file, (char *) "/Calculation/arpes/OrbitalWeights");
      get_hdf5 <double>(k_vectors.data(), file, (char *) "/Calculation/arpes/k_vector");

      file->close();  
      delete file;

      for(unsigned i = 0; i < r.Orb; i++)
        weight(i) = T(weight_test(i));
      
      //std::cout << "weights: " << weight << "\n";

}
#pragma omp barrier




     Eigen::Array<double, -1, -1> k_transposed;
     k_transposed = k_vectors.transpose();
     ARPES(NumDisorder, NumMoments, k_transposed, weight);
    }

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
