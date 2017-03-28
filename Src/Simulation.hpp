#ifndef _SIMULATION_HPP
#define _SIMULATION_HPP


template <typename T,unsigned D>
class GlobalSimulation {
private:
  std::vector<T> Border;
  LatticeStructure <D> rglobal;  
public:
  GlobalSimulation( char *name ) : rglobal(name) {    
    Border.resize( rglobal.get_BorderSize() );
    std::fill(Border.begin(), Border.end(), 0);
    omp_set_num_threads(rglobal.n_threads);
    std::cout << rglobal.thread_id << std::endl;
    std::cout.flush();
#pragma omp parallel default(shared) 
    {
      Simulation<T,D> simul(rglobal, name, Border); 
      simul.Measure_Dos();
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
public:
  Simulation(LatticeStructure <D> rglobal,
	     char *name ,
	     std::vector<T> & borders) : r(rglobal), h(rnd, r, name), GlobalBorder(borders)  {
    
    r.thread_id = omp_get_thread_num();
    Border.resize(GlobalBorder.size()/r.n_threads);
  };
  
  void Measure_Dos() {
#pragma omp barrier     
    KPM_Vector<T,D> phi0( 1, h, r, Border, GlobalBorder);

    phi0.test_boundaries_system();
    //  phi0.initiate_vector();
    // phi0.template Multiply<0>(0);
  }
};


  
