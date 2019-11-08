#include "Generic.hpp"
#include "ComplexTraits.hpp"
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

template<typename T,unsigned D>
Simulation<T,D>::Simulation(char *filename, GLOBAL_VARIABLES <T> & Global1): r(filename),  Global(Global1), name(filename), h(name, r, Global1)  {
  // Initializes the Hamiltonian h, an instance of Lattice Structure r, 
  // and an instance of GLOBAL_VARIABLES Global1
  ghosts.resize(Global.ghosts.size()/r.n_threads);
}


template <typename T,unsigned D>
void Simulation<T,D>::generalized_velocity(KPM_Vector<T,D>* kpm0, KPM_Vector<T,D>* kpm1, std::vector<std::vector<unsigned>> indices, int pos){
  // Check which generalized velocity operator needs to be calculated. 
  // reads from kpm1 and writes on kpm0
  
  T * kpm1data = kpm1->v.col(kpm1->get_index()).data();
  T * kpm0data = kpm0->v.col(kpm0->get_index()).data();
  
  switch(indices.at(pos).size()){
  case 0:
    break;
  default:
    kpm0->Velocity(kpm0data, kpm1data, pos); 											
    break;
  }
}


template <typename T,unsigned D>
void Simulation<T,D>::cheb_iteration(KPM_Vector<T,D>* kpm, long int current_iteration){
  // Performs a chebyshev iteration
  if(current_iteration == 0){
    kpm->template Multiply<0>(); 
  } else {
    kpm->template Multiply<1>(); 
  }
}



template <typename T,unsigned D>	
std::vector<std::vector<unsigned>> Simulation<T,D>::process_string(std::string indices_string){
  // First of all, split the indices string by commas ','
  std::vector<std::string> strings;
  int end_pos = 0;
  while(end_pos != -1){
    end_pos = indices_string.find(',');
    strings.push_back(indices_string.substr(0, end_pos));
    indices_string = indices_string.substr(end_pos + 1, indices_string.size() - end_pos - 1);
  }
  
  int dim = strings.size();
  /*
    for(int i = 0; i < dim; i++)
    std::cout << strings.at(i) << "\n";*/
  
  
  std::vector<std::vector<unsigned>> indices;
  
  for(int i = 0; i < dim; i++){
			
    int len_str = strings.at(i).size();
    int single_digit;
    std::vector<unsigned> temp;
			
    for(int j = 0; j < len_str; j++){
      char single_char = strings.at(i)[j];
				
      if(single_char == 'x'){
        single_digit = 0;
      } else {
        if(single_char == 'y'){
          single_digit = 1;
        } else {
          if(single_char == 'z'){
            single_digit = 2;
          } else {
            // This block should never run
            std::cout << "Please enter a valid expression.\n";
            exit(1);
          }
        } 
      }
      temp.push_back(single_digit);
				
    }
    	
    indices.push_back(temp);
			   
  }
			
			
  return indices;
}


template <typename T,unsigned D>
double Simulation<T,D>::time_kpm(int N_average){
  debug_message("Entered time_kpm");
  //This function serves to provide an estimate of the time it takes for each kpm iteration      
#pragma omp barrier
  KPM_Vector<T,D> kpm0(1, *this);
  KPM_Vector<T,D> kpm1(2, *this);

  kpm0.initiate_vector();
  kpm1.set_index(0);
  kpm1.v.col(0) = kpm0.v.col(0);
  kpm1.template Multiply<0>();
  
  auto t0 =  std::chrono::high_resolution_clock::now();
  for(int i = 0; i < N_average; i++)
    kpm1.template Multiply<1>(); 
  auto t1 =  std::chrono::high_resolution_clock::now();
  
  std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);

		
#pragma omp barrier
  debug_message("Left time_kpm");
  return double(time_span.count())/N_average;
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
