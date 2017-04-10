#ifndef _SIMULATION_HPP
#define _SIMULATION_HPP


template <typename T,unsigned D>
class GlobalSimulation {
private:
  std::vector<T> Border;
  LatticeStructure <D> rglobal;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mu;
  const int MUSIZE;
public:
  GlobalSimulation( char *name ) : rglobal(name), MUSIZE(100) {    
    Border.resize( rglobal.get_BorderSize() );
    std::fill(Border.begin(), Border.end(), 0);
    mu = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > :: Zero(1, MUSIZE); 
    omp_set_num_threads(rglobal.n_threads);

#pragma omp parallel default(shared)
    {
      Simulation<T,D> simul(rglobal, name, Border, mu);
      simul.Measure_Dos(10000);
    }
  };
  
};
#endif


template <typename T,unsigned D>
class Simulation  {
  KPMRandom <T> rnd;
  LatticeStructure <D>  r;      
  Hamiltonian<T,D> h;
  std::vector<T> Border;
  std::vector<T> & GlobalBorder;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mu;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & mu_global;
  char *filename;
public:
  Simulation(LatticeStructure <D> rglobal,
	     char *name ,
	     std::vector<T> & borders,
	     Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & m ) :
    r(rglobal), h(rnd, r, name), GlobalBorder(borders), mu_global(m), filename(name)  {
    r.thread_id = omp_get_thread_num();    
    Border.resize(GlobalBorder.size()/r.n_threads);
    mu = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > :: Zero(mu_global.rows(), mu_global.cols()); 
  };
  
  
  void Measure_Dos(int NMedias ) {
    typedef typename extract_value_type<T>::value_type value_type;
    KPM_Vector<T,D> phi0(1, h, r, Border, GlobalBorder);
    KPM_Vector<T,D>  phi(2, h, r, Border, GlobalBorder);
    
    //    phi0.test_boundaries_system();

    for(int medias = 0; medias < NMedias; medias++)
      {
	phi0.initiate_vector();
	phi.set_index(0);
	phi.v.col(0) = phi0.v.col(0);
	phi.Exchange_Boundaries();

	phi.template Multiply<0>(0);
	mu.block(0,0,1,2) +=  (phi0.v.adjoint() * phi.v - mu.block(0,0,1,2))/value_type(medias + 1);
	
	for(int m = 2; m < mu.cols(); m += 2)
	  {	    
	    phi.template Multiply<1>(0);
	    phi.template Multiply<1>(0);
	    mu.block(0,m,1,2) +=  (phi0.v.adjoint() * phi.v - mu.block(0,m,1,2))/value_type(medias + 1);
	  }		
      }
#pragma omp critical
    mu_global += mu;
	
#pragma omp master
    {
      H5::H5File file = H5::H5File(filename, H5F_ACC_RDWR);
      Eigen::MatrixXd  dMUR = Eigen::MatrixXd::Zero(1,mu_global.cols()); //mu_global.real();
      Eigen::MatrixXd  dMUI = Eigen::MatrixXd::Zero(1,mu_global.cols()); //mu_global.imag();
      H5::DataSet dataset;
      H5::DataSpace dataspace;
      H5::IntType datatype(H5::PredType::NATIVE_DOUBLE );
      H5::DSetCreatPropList plist;
      const int RANK = 1;
      hsize_t    dims[RANK], chunk_dims[RANK]; // dataset dimensions

      for(int i = 0; i < mu_global.cols(); i++)
	{
	  dMUR(0,i) = std::real(mu_global(0,i));
	  dMUI(0,i) = std::imag(mu_global(0,i));
	}
      
      dims[0] = chunk_dims[0] = mu_global.cols();

      
      dataspace = H5::DataSpace( RANK, dims );
      plist.setChunk(RANK, chunk_dims);
    
     
      dataspace = H5::DataSpace( 1, dims );
      plist.setChunk(1, chunk_dims);
      plist.setDeflate(6);
      //dataset = file.createDataSet( "MUR", datatype, dataspace );
      dataset = file.openDataSet( "MUR");
      dataset.write(dMUR.data(), H5::PredType::NATIVE_DOUBLE );  
      //dataset = file.createDataSet( "MUI", datatype, dataspace );
      dataset = file.openDataSet( "MUI");
      dataset.write(dMUI.data(), H5::PredType::NATIVE_DOUBLE );  
  

      file.close();
    }
#pragma omp barrier
  }
};


  
