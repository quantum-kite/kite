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

#if !(COMPILE_WAVEPACKET)
#warning "Cannot compile SimulationGaussianWavepacket.cpp. This error is not fatal, but KITE will not be able to run GaussianWavepacket(). A more recent version of gcc (8.0) is required."
#endif

template <typename T, unsigned DIM>
void Simulation<T, DIM>::calc_wavepacket(){
    // Make sure that all the threads are ready before opening any files
    // Some threads could still be inside the Simulation constructor
    // This barrier is essential
    //

#pragma omp barrier

    // Check if the Gaussian_Wave_Packet needs to be calculated
    bool local_calculate_wavepacket;
#pragma omp master
    {
        Global.calculate_wavepacket = 0;
        H5::H5File * file = new H5::H5File(name, H5F_ACC_RDONLY);
      try{
        int dummy_var;
        get_hdf5<int>(&dummy_var, file, (char *) "/Calculation/gaussian_wave_packet/NumDisorder");
        Global.calculate_wavepacket = 1;
      } catch(H5::Exception& e) {debug_message("Wavepacket: no need to calculate.\n");}
        file->close();  
        delete file;
    }
#pragma omp barrier
#pragma omp critical
    local_calculate_wavepacket = Global.calculate_wavepacket;
      
    // Now calculate it
#pragma omp barrier
    if(local_calculate_wavepacket){
      Gaussian_Wave_Packet();
    }
}

template <typename T, unsigned D>
void Simulation<T,D>::Gaussian_Wave_Packet(){
#if COMPILE_WAVEPACKET
#pragma omp master
      {
        std::cout << "Calculating Wavepacket.\n";
      }
#pragma omp barrier
  ComplexTraits<T> CT;
  KPM_Vector<T,D> phi (2, *this), sum_ket(1u,*this);
  int NumDisorder, NumMoments, NumPoints;
    
  T II   = CT.assign_value(double(0), double(1));
  T zero = CT.assign_value(double(0),  double(0));
  T one  = CT.assign_value(double(1),  double(0));
  std::vector<double> times;
  
  Eigen::Array<T,-1,-1> avg_x, avg_y, avg_z, avg_ident;
  Eigen::Matrix<T, 2, 2> ident, spin_x, spin_y, spin_z;

  Eigen::Map<Eigen::Matrix<T,-1,-1>> vket (sum_ket.v.data(), r.Sized/2, 2);
  Eigen::Map<Eigen::Matrix<T,-1,-1>> vtmp (phi.v.data()    , r.Sized/2, 2);
  Eigen::Map<Eigen::Matrix<T,-1,-1>> vtmp1(phi.v.data()    , r.Sized  , 1);
  float timestep;
  double width;
  Eigen::Array<T, -1, -1> avg_results;
  Eigen::Array<T, -1, -1> results(2*D, 1);
  H5::DataSet * dataset;
  H5::DataSpace * dataspace;
  hsize_t dim[2];
  Eigen::Matrix <double,-1, -1> k_vector;
  Eigen::Matrix <double ,1, 2> vb;
  Eigen::Matrix <T,-1, -1>        spinor;

  ident <<  one, zero,
    zero, one;
  
  spin_x << zero, one,
    one, zero;
    
  spin_y << zero, -II,
    II, zero;
    
  spin_z << one, zero,
    zero, -one;

  //Load bra and ket
#pragma omp critical
  {
    H5::H5File * file  = new H5::H5File(name, H5F_ACC_RDONLY);
    dataset            = new H5::DataSet(file->openDataSet("/Calculation/gaussian_wave_packet/k_vector")  );
    dataspace          = new H5::DataSpace(dataset->getSpace());
    dataspace -> getSimpleExtentDims(dim, NULL);
    dataspace->close(); delete dataspace;
    dataset->close();   delete dataset;

    k_vector  = Eigen::Matrix<double,-1, -1>::Zero(dim[1],dim[0]);
    spinor    = Eigen::Matrix<     T,-1, -1>::Zero(r.Orb,dim[0]);
    vb = Eigen::Matrix<     double,1, 2>::Zero(2);
    get_hdf5    <int>(&NumDisorder,    file, (char *) "/Calculation/gaussian_wave_packet/NumDisorder");
    get_hdf5    <int>(&NumMoments,     file, (char *) "/Calculation/gaussian_wave_packet/NumMoments" );
    get_hdf5    <int>(&NumPoints,      file, (char *) "/Calculation/gaussian_wave_packet/NumPoints"  );
    get_hdf5  <float>(&timestep,       file, (char *) "/Calculation/gaussian_wave_packet/timestep"   );
    get_hdf5 <double>(&width,          file, (char *) "/Calculation/gaussian_wave_packet/width"      );
    get_hdf5      <T>(spinor.data(),   file, (char *) "/Calculation/gaussian_wave_packet/spinor");
    get_hdf5     <double>(vb.data(),   file, (char *) "/Calculation/gaussian_wave_packet/mean_value");
    get_hdf5 <double>(k_vector.data(), file, (char *) "/Calculation/gaussian_wave_packet/k_vector");

    file->close();  delete file;
  }
#pragma omp barrier
  avg_x       = Eigen::Matrix<T,-1,-1>::Zero(NumPoints,1);
  avg_y       = Eigen::Matrix<T,-1,-1>::Zero(NumPoints,1);
  avg_z       = Eigen::Matrix<T,-1,-1>::Zero(NumPoints,1);
  avg_ident   = Eigen::Matrix<T,-1,-1>::Zero(NumPoints,1);
  avg_results = Eigen::Array<T, -1,-1>::Zero(2*D, NumPoints);
    
#pragma omp master
  {
    Global.avg_x       = Eigen::Matrix<T,-1,1>::Zero(NumPoints,1);
    Global.avg_y       = Eigen::Matrix<T,-1,1>::Zero(NumPoints,1);
    Global.avg_z       = Eigen::Matrix<T,-1,1>::Zero(NumPoints,1);
    Global.avg_ident   = Eigen::Matrix<T,-1,1>::Zero(NumPoints,1);
    Global.avg_results = Eigen::Array<T,-1,-1>::Zero(2*D,NumPoints);
  }
    
    
  NumMoments = (NumMoments/2)*2;
  Eigen::Matrix<T,-1,1> m(NumMoments);
  for(unsigned n = 0; n < unsigned(NumMoments); n++)
    m(n) = value_type((n == 0 ? 1 : 2 )*std::cyl_bessel_j(n, timestep )) * T(pow(-II,n));
    
  for(int id = 0; id < NumDisorder; id++)
    {
      sum_ket.set_index(0);
      sum_ket.v.setZero();
      sum_ket.build_wave_packet(k_vector, spinor, width, vb);
      h.generate_disorder();	
      sum_ket.empty_ghosts(0);
      for(unsigned t = 0; t < unsigned(NumPoints); t++)
        {
          if(t > 0)
            {
              phi.v.setZero();
              phi.set_index(0);
              phi.v.col(0) = sum_ket.v.col(0);
              phi.Exchange_Boundaries();
              cheb_iteration(&phi, 0); // multiply by H
              sum_ket.v.col(0) = phi.v * m.segment(0, 2);
              for(unsigned n = 2; n < unsigned(NumMoments); n += 2)
                {
                  cheb_iteration(&phi, n - 1);
                  cheb_iteration(&phi, n);
                  sum_ket.v.col(0) += phi.v * m.segment(n,2);
                }
            }
          sum_ket.empty_ghosts(0);
	    
          // In the multiplication of a matrix the number of columns of the first should be equal to
          // the number of rows of the second, because the spin is organized by columns we have:

          //	    vtmp = (vket * ident.transpose()) - vket;
          //	    auto x3 = sum_ket2.v.adjoint() * vtmp1;
          auto x3 = sum_ket.get_point();
          avg_ident(t) += (x3 - avg_ident(t) ) /T(id + 1);


          vtmp = vket * spin_x.transpose();
          auto x0 = sum_ket.v.adjoint() * vtmp1;
          avg_x(t) += (x0(0,0) - avg_x(t) ) /T(id + 1);
	    
          vtmp = vket * spin_y.transpose();
          auto x1 = sum_ket.v.adjoint() * vtmp1;
          avg_y(t) += (x1(0,0) - avg_y(t) ) /T(id + 1);
	    
          vtmp = vket * spin_z.transpose();
          auto x2 = sum_ket.v.adjoint() * vtmp1;

          avg_z(t) += (x2(0,0) - avg_z(t) ) /T(id + 1);
          phi.measure_wave_packet(sum_ket.v.data(), sum_ket.v.data(), results.data());
          avg_results.col(t) += (results - avg_results.col(t) ) /T(id + 1); 
        }
	
    }

#pragma omp critical
  {      
    Global.avg_x += avg_x;
    Global.avg_y += avg_y;
    Global.avg_z += avg_z;
    Global.avg_ident += avg_ident;
    Global.avg_results += avg_results;
  }
#pragma omp barrier


    
#pragma omp master
  {
    std::cout << name << std::endl;
    H5::H5File * file = new H5::H5File(name, H5F_ACC_RDWR);
    write_hdf5(Global.avg_x, file, (char *) "/Calculation/gaussian_wave_packet/Sx");
    write_hdf5(Global.avg_y, file, (char *) "/Calculation/gaussian_wave_packet/Sy");
    write_hdf5(Global.avg_z, file, (char *) "/Calculation/gaussian_wave_packet/Sz");
    write_hdf5(Global.avg_ident, file, (char *) "/Calculation/gaussian_wave_packet/Id");


    for(unsigned i = 0; i < D; i++)
      {
        std::string orient = "xyz";
        char name[200];
        avg_z.col(0) = Global.avg_results.row(2*i) ;
        sprintf(name,"/Calculation/gaussian_wave_packet/mean_value%c", orient.at(i));
        write_hdf5(avg_z, file, name);
      }
      
    for(unsigned i = 0; i < D; i++)
      {
        std::string orient = "xyz";
        char name[200];
        avg_z.col(0) = Global.avg_results.row(2*i + 1) - Global.avg_results.row(2*i)*Global.avg_results.row(2*i);
        sprintf(name,"/Calculation/gaussian_wave_packet/Var%c", orient.at(i));
			  
        write_hdf5(avg_z, file, name);
      }
      
    file->close();      
    delete file;
  }
#pragma omp barrier
#endif
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
