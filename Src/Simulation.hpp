#ifndef _SIMULATION_HPP
#define _SIMULATION_HPP

template <typename T,unsigned D>
class GlobalSimulation {
private:
  GLOBAL_VARIABLES <T> Global;
  LatticeStructure <D> rglobal;
  std::vector <int> Quantities, NMoments, NRandomV, NDisorder; 
public:
  GlobalSimulation( char *name ) : rglobal(name) {
    Global.ghosts.resize( rglobal.get_BorderSize() );
    std::fill(Global.ghosts.begin(), Global.ghosts.end(), 0);
    
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
    
    H5::H5File * file         = new H5::H5File(name, H5F_ACC_RDONLY);
    H5::DataSet * dataset     = new H5::DataSet(file->openDataSet("/Calculation/FunctionNum"));
    H5::DataSpace * dataspace = new H5::DataSpace(dataset->getSpace());
    size_t NQuantities        = dataspace->getSimpleExtentNpoints();
    Quantities.resize (NQuantities);
    NMoments.resize   (NQuantities);
    NRandomV.resize   (NQuantities);
    NDisorder.resize  (NQuantities);

    get_hdf5<int>(Quantities.data(), file, (char *)   "/Calculation/FunctionNum");
    get_hdf5<int>(NRandomV.data(),   file, (char *)   "/Calculation/NumRandoms");
    get_hdf5<int>(NMoments.data(),   file, (char *)   "/Calculation/NumMoments");
    get_hdf5<int>(NDisorder.data(),  file, (char *)   "/Calculation/NumDisorder");
    for(unsigned i = 0; i < NMoments.size(); i++)
      NMoments.at(i) = 2*(NMoments.at(i)/2);

    delete dataspace;
    delete dataset;
    delete file;

    bool DOS = any_of(Quantities.begin(), Quantities.end(), std::bind2nd(std::equal_to<int>(), 1));
    int index_dos = 0;
    if(DOS) {
      while(Quantities.at(index_dos) != 1) index_dos++;
      Global.mu = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > :: Zero(1, NMoments.at(index_dos));
    }

    omp_set_num_threads(rglobal.n_threads);
#pragma omp parallel default(shared)
    {
      Simulation<T,D> simul(name, Global);	    
      if(DOS) simul.Measure_Dos(NRandomV.at(index_dos), NDisorder.at(index_dos) );
    }
    
  };

  
};
#endif


template <typename T,unsigned D>
class Simulation  {
public:
  KPMRandom <T>          rnd;
  std::vector<T>         ghosts;
  LatticeStructure <D>   r;      
  GLOBAL_VARIABLES <T> & Global;
  char                 * name;
  Hamiltonian<T,D>       h;
  Simulation(char *filename, GLOBAL_VARIABLES <T> & Global1): r(filename),  Global(Global1), name(filename), h(*this)  {
    rnd.init_random();
    ghosts.resize(Global.ghosts.size()/r.n_threads);
  };
  
  
  void Measure_Dos(int & NRandomV, int & NDisorder) {
    typedef typename extract_value_type<T>::value_type value_type;
    Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> mu = Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic > :: Zero(Global.mu.rows(), Global.mu.cols());
    KPM_Vector<T,D>  phi0(1, *this);
    KPM_Vector<T,D>  phi (2, *this);

    //    phi0.test_boundaries_system();

#pragma omp barrier
    long average = 0;

    for(int disorder = 0; disorder <  NDisorder; disorder++)
      {
	h.generate_disorder();
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

#pragma omp critical
	Global.mu += mu;
#pragma omp barrier
	
	
#pragma omp master
	{
	  H5::H5File * file = new H5::H5File(name, H5F_ACC_RDWR);
	  write_hdf5(Global.mu, file, "MU");
	  delete file;
	  Global.mu.setZero();
	}
#pragma omp barrier	
      }
#pragma omp barrier	
  }
};


  
