#ifndef _SIMULATION_HPP
#define _SIMULATION_HPP


template <typename T,unsigned D>
class GlobalSimulation {
private:
  std::vector<T> Border;
  LatticeStructure <D> rglobal;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> mu;
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
  Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> mu;
  Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> & mu_global;
  char *filename;
public:
  Simulation(LatticeStructure <D> rglobal,
	     char *name ,
	     std::vector<T> & borders,
	     Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> & m ) :
    r(rglobal), h(rnd, r, name), GlobalBorder(borders), mu_global(m), filename(name)  {
    r.thread_id = omp_get_thread_num();    
    Border.resize(GlobalBorder.size()/r.n_threads);
    mu = Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic > :: Zero(mu_global.rows(), mu_global.cols()); 
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
	mu.matrix().block(0,0,1,2) +=  (phi0.v.adjoint() * phi.v - mu.matrix().block(0,0,1,2))/value_type(medias + 1);
	
	for(int m = 2; m < mu.cols(); m += 2)
	  {	    
	    phi.template Multiply<1>(0);
	    phi.template Multiply<1>(0);
	    mu.matrix().block(0,m,1,2) +=  (phi0.v.adjoint() * phi.v - mu.matrix().block(0,m,1,2))/value_type(medias + 1);
	  }		
      }
#pragma omp critical
    mu_global += mu;
	
#pragma omp master
    {
      H5::H5File file = H5::H5File(filename, H5F_ACC_RDWR);
      write_hdf5(mu_global, file, "MU");
      file.close();
    }
#pragma omp barrier
  }
};


  
