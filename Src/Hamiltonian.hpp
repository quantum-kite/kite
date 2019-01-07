/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

template <typename T, unsigned D>
class Hamiltonian;

#include "HamiltonianRegular.hpp"
#include "HamiltonianVacancies.hpp"
#include "HamiltonianDefects.hpp"


template <typename T, unsigned D>
class Hamiltonian {
private:
public:
  KPMRandom <T>          rnd;
  char                 *name;
  LatticeStructure<D>    & r;
  typedef typename extract_value_type<T>::value_type value_type;
  GLOBAL_VARIABLES <T>  & Global;
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
  Hamiltonian(char *name,  LatticeStructure<D> & rr, GLOBAL_VARIABLES <T> &  );
  void generate_disorder();
  void build_structural_disorder();
  void build_vacancies_disorder();
  void build_Anderson_disorder();
  void build_velocity(std::vector<unsigned> & components, unsigned n);
  void distribute_AndersonDisorder();
};





