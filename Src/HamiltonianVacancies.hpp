/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

template <typename T,unsigned D>
struct Vacancy_Operator {
  LatticeStructure <D>                   & r;
  Simulation <T,D>                       & simul;
  // vector of vectors with positions in the lattice of Vacancies for each stride block
  std::vector <std::vector<std::size_t>> position;                   
  std::vector <double>                   concentration;
  std::vector <std::vector<int>>         orbitals;
  std::vector <std::size_t>              vacancies_with_defects; 
  
  Vacancy_Operator(Simulation <T,D> & sim) : r(sim.r), simul(sim), position(sim.r.NStr)
  {
  };
  
  
  void generate_disorder()
  {
    Coordinates<std::size_t,D + 1> latt(r.ld), Latt(r.Ld), latStr(r.lStr);
    
    // Clean former vacancy distribution
    for(unsigned i = 0; i < r.NStr ; i++)
      position.at(i).clear();
    vacancies_with_defects.clear();
    // Distribute Vacancies
    
    for( std::size_t i = 0; i < r.N ; i++)
      for(unsigned k = 0; k < concentration.size(); k++)
	if(simul.rnd.get() < concentration[k])
	  for(auto o = orbitals.at(k).begin(); o != orbitals.at(k).end(); o++)
	    {
	      latt.set_coord(i + std::size_t(*o) * r.N );         // Set Vacancy position inside the domain
	      r.convertCoordinates(latStr,latt);                  // Get stride position
	      r.convertCoordinates(Latt,latt);                    // Get Domain coordinates
	      position.at(latStr.index).push_back(Latt.index);    // store the position by stride 
	    }
  };
  
  void add_model(double p, std::vector <int> & orb)
  {
    orbitals.push_back(orb);
    concentration.push_back(p);
  }
  
  void add_conflict_with_defect(std::size_t element, unsigned istride)
  {
    std::vector<std::size_t> & v = position.at(istride);  
    for(unsigned i = 0; i < v.size(); i++)
      if(element == v.at(i))
	vacancies_with_defects.push_back(element);
  }
  
  bool test_vacancy(Coordinates<std::size_t,D + 1> & Latt)
  {
    /*
      1 if is a vacancy
      0 if not
    */
    Coordinates<std::size_t,D + 1> latStr(r.lStr);
    r.convertCoordinates(latStr, Latt);
    auto & vc = position.at(latStr.index);
    if( !any_of(vc.begin(), vc.end(), std::bind2nd(std::equal_to<std::size_t>(), Latt.index)))
      return 0;
    else
      return 1;
  }

  void test_field( T * phi0 )
  {
    // The field should be zero in the vacancies
    for(unsigned i = 0; i  < r.NStr; i++)
      for(auto vc = position.at(i).begin(); vc != position.at(i).end(); vc++)
	if( abs(phi0[*vc]) > __DBL_EPSILON__)
	  {
	    std::cout << "Disparate" << std::endl;
	    exit(1);
	  }
  }
  
};
