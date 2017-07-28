#ifndef _SIMULATION_HPP
#define _SIMULATION_HPP

template <typename T,unsigned D>
class GlobalSimulation {
private:
  std::vector<T> Border;
  LatticeStructure <D> rglobal;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> mu;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> gammaxx;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> gammaxy;
  std::vector <int> Quantities, NMoments, NRandomV, NDisorder; 
public:
  GlobalSimulation( char *name ) : rglobal(name) {
    Border.resize( rglobal.get_BorderSize() );
    std::fill(Border.begin(), Border.end(), 0);
    
    /*
     *   /Calculation/FunctionNum:
     *
     *   DOS - denstity of states == 1,
     *   CondXX - conductivity in xx direction == 2,
     *   CondXY - conductivity in xy direction == 3,
     *   OptCond - optical conductivity == 4
     *   SpinCond - spin conductivity == 5

     * /Calculation/NumRandoms   : number of different random vector realisations,
     * /Calculation/NumMoments   : number of moments for the calculation,
     * /Calculation/NumDisorder  : number of disorder realisations.
     */
    
    H5::H5File * file        = new H5::H5File(name, H5F_ACC_RDONLY);
    H5::DataSet * dataset     = new H5::DataSet(file->openDataSet("/Calculation/FunctionNum"));
    H5::DataSpace * dataspace = new H5::DataSpace(dataset->getSpace());
    size_t NQuantities        = dataspace->getSimpleExtentNpoints();
    Quantities.resize(NQuantities);
    NMoments.resize(NQuantities);
    NRandomV.resize(NQuantities);
    NDisorder.resize(NQuantities);

    get_hdf5<int>(Quantities.data(), file, (char *) "/Calculation/FunctionNum");
    get_hdf5<int>(NRandomV.data(), file, (char *) "/Calculation/NumRandoms");
    get_hdf5<int>(NMoments.data(), file, (char *) "/Calculation/NumMoments");
    get_hdf5<int>(NDisorder.data(), file, (char *) "/Calculation/NumDisorder");
    delete dataspace;
    delete dataset;
    delete file;

    bool DOS = any_of(Quantities.begin(), Quantities.end(), std::bind2nd(std::equal_to<int>(), 1));
    int index_dos = 0;
    if(DOS)
      {
	while(Quantities.at(0) != 1)
	  index_dos++;
	mu = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > :: Zero(1, NMoments.at(index_dos));
      }
    omp_set_num_threads(rglobal.n_threads);
#pragma omp parallel default(shared)
    {
      LatticeStructure <D> r(rglobal);                     // Create a local copy of the lattice Structure 
      r.thread_id = omp_get_thread_num();        
      Simulation<T,D> simul(r, name, Border, mu);      
      if(DOS) simul.Measure_Dos(NRandomV.at(index_dos), NDisorder.at(index_dos) );
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
    rnd.init_random();
  };
  
  
  void Measure_Dos(int & NRandomV, int & NDisorder) {
    typedef typename extract_value_type<T>::value_type value_type;
    KPM_Vector<T,D> phi0(1, h, r, Border, GlobalBorder);
    KPM_Vector<T,D>  phi(2, h, r, Border, GlobalBorder);

    
    
    //    phi0.test_boundaries_system();
    long average = 0;
    for(int disorder = 0; disorder < NDisorder; disorder++)
      {
	h.distribute_AndersonDisorder();
	for(int randV = 0; randV < NRandomV; randV++)
	  {
	    phi0.initiate_vector();
	    phi.set_index(0);
	    phi.v.col(0) = phi0.v.col(0);
	    phi.Exchange_Boundaries();
	    
	    phi.template Multiply<0>();
	    mu.matrix().block(0,0,1,2) +=  (phi0.v.adjoint() * phi.v - mu.matrix().block(0,0,1,2))/value_type(average + 1);
	    
	    for(int m = 2; m < mu.cols(); m += 2)
	      {	    
		phi.template Multiply<1>();
		phi.template Multiply<1>();
		mu.matrix().block(0,m,1,2) +=  (phi0.v.adjoint() * phi.v - mu.matrix().block(0,m,1,2))/value_type(average + 1);
	      }
	    average++;
	  }
      }
#pragma omp critical
    mu_global += mu;
#pragma omp barrier

    
#pragma omp master
    {
      H5::H5File * file = new H5::H5File(filename, H5F_ACC_RDWR);
      write_hdf5(mu_global, file, "MU");
      delete file;
    }
#pragma omp barrier
  }
};


  
