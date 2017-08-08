#ifndef _HAMILTONIAN_HPP
#define _HAMILTONIAN_HPP

#include "HamiltonianDefects.hpp"
#include "HamiltonianRegular.hpp"

template <typename T, unsigned D>
class Hamiltonian {
private:
  LatticeStructure<D> & r;
public:
  typedef typename extract_value_type<T>::value_type value_type;
  Simulation<T,D> & simul;  
  Periodic_Operator<T,D> hr;

  /* Anderson disorder */
  std::vector <int> orb_num;
  std::vector <int> model;
  std::vector <double> mu;
  std::vector <double> sigma;
  
  std::vector<value_type> U_Orbital;
  std::vector<value_type> U_Anderson;      // Local disorder  
  std::vector<int> Anderson_orb_address;
  
  /*   Structural disorder    */
  Defect_Operator<T,D> hd;
  
  Hamiltonian (Simulation<T,D> & sim) : r(sim.r), simul(sim) , hr(sim), hd(sim)
  {  
    /* Anderson disorder */
    build_Anderson_disorder();
    distribute_AndersonDisorder();
  };

  void generate_disorder()
  {
    distribute_AndersonDisorder();
    hd.generate_disorder();
  };
  
  void build_Anderson_disorder() {
    /*
     * Gaussian      : 1
     * Uniform       : 2
     * Deterministic : 3
     */
#pragma omp critical
    {
      H5::H5File    *file      = new H5::H5File(simul.name, H5F_ACC_RDONLY);
      H5::DataSet   dataset    = H5::DataSet(file->openDataSet("/Hamiltonian/Disorder/OrbitalNum"));
      H5::DataSpace dataspace  = dataset.getSpace();
      size_t        m          = dataspace.getSimpleExtentNpoints();

      orb_num.resize(m);
      model.resize(m);
      mu.resize(m);
      sigma.resize(m);

      get_hdf5<int>(orb_num.data(), file, (char *) "/Hamiltonian/Disorder/OrbitalNum");               // read the orbitals that have local disorder
      get_hdf5<int> (model.data(), file, (char *) "/Hamiltonian/Disorder/OnsiteDisorderModelType");   // read the the type  of local disorder
      get_hdf5<double> (mu.data(), file, (char *) "/Hamiltonian/Disorder/OnsiteDisorderMeanValue");   // read the the mean value
      get_hdf5<double> (sigma.data(), file, (char *) "/Hamiltonian/Disorder/OnsiteDisorderMeanStdv"); // read the the variance
      delete file;
    }
    
    Anderson_orb_address.resize(r.Orb);
    U_Orbital.resize(r.Orb);
    
    std::fill_n ( Anderson_orb_address.begin(), r.Orb, -2 );
    std::fill_n ( U_Orbital.begin(), r.Orb,  0 );
    
    int sum = 0;
    for (unsigned i = 0; i < model.size(); i++)
      {
	if(model.at(i) < 3)
	  {
	    Anderson_orb_address.at( orb_num.at(i) ) = sum;
	    sum++;
	  }
	
	if(model.at(i) == 3) // Deterministic
	  {
	    Anderson_orb_address.at( orb_num.at(i) ) = -1;
	    U_Orbital.at( orb_num.at(i) ) = mu.at(i);
	  }
      }
    
    if(sum > 0)
      U_Anderson.resize( sum * r.Nd);    
  }
  
  
  void distribute_AndersonDisorder()
  {
    
    /*
     * Gaussian      : 1
     * Uniform       : 2
     * Deterministic : 3
     */
    
    int sum = 0;
    for (unsigned i = 0; i < model.size(); i++)
      if(model.at(i) == 1 )
	{
	  for(unsigned j = 0; j < r.Nd; j++)
            U_Anderson.at(sum + j) = simul.rnd.gaussian(mu.at(i), sigma.at(i)) ;
	  sum += r.Nd;
	}
      else if ( model.at(i) == 2 )
	{
	  for(unsigned j = 0; j < r.Nd; j++)
            U_Anderson.at(sum + j) = simul.rnd.uniform(mu.at(i), sigma.at(i)) ;
	  sum += r.Nd;
	}
  }
};




#endif
