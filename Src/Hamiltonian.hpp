/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

#ifndef _HAMILTONIAN_HPP
#define _HAMILTONIAN_HPP

#include "HamiltonianDefects.hpp"
#include "HamiltonianRegular.hpp"
#include "HamiltonianVacancies.hpp"
extern "C" herr_t getMembers(hid_t loc_id, const char *name, void *opdata);

template <typename T, unsigned D>
class Hamiltonian {
private:
  LatticeStructure<D> & r;
public:
  typedef typename extract_value_type<T>::value_type value_type;
  Simulation<T,D> & simul;  
  Periodic_Operator<T,D>  hr;

  /* Anderson disorder */
  std::vector <int> orb_num;
  std::vector <int> model;
  std::vector <double> mu;
  std::vector <double> sigma;
  
  std::vector<value_type> U_Orbital;
  std::vector<value_type> U_Anderson;      // Local disorder  
  std::vector<int> Anderson_orb_address;
  
  /*   Structural disorder    */
  std::vector <bool>                   cross_mozaic;
  std::vector <std::size_t>            cross_mozaic_indexes;
  std::vector < Defect_Operator<T,D>>  hd;
  Vacancy_Operator<T,D>                hV;
  
  Hamiltonian (Simulation<T,D> & sim) : r(sim.r), simul(sim) , hr(sim), cross_mozaic(sim.r.NStr), hV(sim)
  {  
    /* Anderson disorder */
    build_Anderson_disorder();
    build_vacancies_disorder();    
    build_structural_disorder();
  };
  
  void generate_disorder()
  {
    distribute_AndersonDisorder();
    for(std::size_t istr = 0; istr < r.NStr; istr++)
      cross_mozaic[istr] = true;
    cross_mozaic_indexes.clear();

    hV.generate_disorder();
    for(auto id = hd.begin(); id != hd.end(); id++)
      id->generate_disorder();
    

  };
  
  void build_structural_disorder()
  {
    
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(simul.name, H5F_ACC_RDONLY);
      // Test if there is a strutural disorder to build
      H5::Group  grp;
      std::vector<std::string> defects;
      try {
	H5::Exception::dontPrint();
	grp = file->openGroup("/Hamiltonian/StructuralDisorder");
	grp.iterateElems(grp.getObjName(), NULL, getMembers, static_cast<void*>(&defects) );
	
	for(auto id = defects.begin(); id != defects.end(); id++)
	  hd.push_back(Defect_Operator<T,D> ( simul, *id, file) );
      }
      catch(H5::Exception& e) {
	// Do nothing
      }
      delete file;
    }
  }

  void build_vacancies_disorder()
  {
    r.SizetVacancies = 0;
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(simul.name, H5F_ACC_RDONLY);
      // Test if there is vacancies to build
      H5::Group  grp;
      double p;
      std::vector<int> orbit;
      std::vector<std::string> vacancies;
      int n;
      
      try {
	H5::Exception::dontPrint();
	grp = file->openGroup("/Hamiltonian/Vacancy");
	// Get the names of the Vacancies Types
	grp.iterateElems(grp.getObjName(), NULL, getMembers, static_cast<void*>(&vacancies)); 
	for(auto id = vacancies.begin(); id != vacancies.end(); id++)
	  {
	    std::string field = *id + std::string("/Concentration");
	    get_hdf5<double> ( &p, file, field );
	    field = *id + std::string("/NumOrbitals");
	    get_hdf5<int> ( &n, file, field );
	    orbit.resize(n);
	    field = *id + std::string("/Orbitals");
	    get_hdf5<int> ( orbit.data(), file, field );
	    hV.add_model(p, orbit);
	  }
      }
      catch(H5::Exception& e) {
	// Do nothing
      }
      delete file;
    }
    
  }
  
  
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
      try {
	get_hdf5<int>(orb_num.data(), file, (char *) "/Hamiltonian/Disorder/OrbitalNum");               // read the orbitals that have local disorder
      	get_hdf5<int> (model.data(), file, (char *) "/Hamiltonian/Disorder/OnsiteDisorderModelType");   // read the the type  of local disorder
      	get_hdf5<double> (mu.data(), file, (char *) "/Hamiltonian/Disorder/OnsiteDisorderMeanValue");   // read the the mean value
      	get_hdf5<double> (sigma.data(), file, (char *) "/Hamiltonian/Disorder/OnsiteDisorderMeanStdv"); // read the the variance
      }
	catch (...){}
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
  
  void build_velocity(std::vector<unsigned> & components, unsigned n)
  {
    hr.build_velocity(components, n);
    for(auto i = hd.begin(); i != hd.end(); i++)
      i->build_velocity(components, n);
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
  template <typename U = T>
  typename std::enable_if<is_tt<std::complex, U>::value, U>::type ghosts_correlation(double phase) {
    std::complex<double> im(0,phase);
    return U(exp(im));
  };
  
  template <typename U = T>
  typename std::enable_if<!is_tt<std::complex, U>::value, U>::type ghosts_correlation(double phase) {
    //std::cout << "Shouldn't be here..\n";
    return 1.0;
  };
};


herr_t getMembers(hid_t loc_id, const char *name, void *opdata)
{
  H5::Group  grp(loc_id);
  std::string Disorder = grp.getObjName();
  std::string  sep = "/";
  std::string Defect = name;
  std::string group = Disorder+sep+name; 
  std::vector<std::string> * v = static_cast<std::vector<std::string> *> (opdata);

  try {
    H5::Exception::dontPrint();
    grp.openGroup(group);
    v->push_back(group);
  }
  catch(H5::Exception& e) {
    // Don't do nothing 
  }
  return 0;
}

#endif
