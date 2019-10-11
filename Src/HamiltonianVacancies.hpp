/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

template <typename T,unsigned D>
struct Vacancy_Operator {
  LatticeStructure <D>                   & r;
  KPMRandom <T>                        & rnd;
  char*                                 name;
  // vector of vectors with positions in the lattice of Vacancies for each tile block
  std::vector <std::vector<std::size_t>> position;                   
  std::vector <double>                   concentration;
  std::vector <std::vector<int>>         orbitals;
  std::vector <std::size_t>              vacancies_with_defects;
  std::vector <std::vector<int>>         positions_fixed;
  Vacancy_Operator(char *, LatticeStructure <D> & , KPMRandom <T> &);
  void generate_disorder();
  void add_model(double p, std::vector <int> & orb, std::vector<int> & positions);
  void add_conflict_with_defect(std::size_t element, unsigned istride);
  bool test_vacancy(Coordinates<std::size_t,D + 1> & Latt);
  void test_field( T * phi0 );
};
