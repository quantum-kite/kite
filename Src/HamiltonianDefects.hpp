/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

template <typename T,unsigned D>
struct Defect_Operator: public ComplexTraits<T> {
  typedef typename extract_value_type<T>::value_type value_type;
  using ComplexTraits<T>::multEiphase;
  double                                   p;                        // Concentration of defects
  unsigned                       NumberNodes;                        // Number of nodes in the deffect
  std::vector <T>                          U;                        // local energies
  std::vector <unsigned>             element;                        // nodes with local energies
  std::vector <T>                    hopping;                        // vector of the non-zero values of the hopping operator operator                // 
  std::vector<std::vector <value_type>>    v;                        // temporary generalized velocities
  std::vector <unsigned>            element1;                        // vector with the nodes 
  std::vector <unsigned>            element2;                        // vector with the nodes
  std::vector <std::ptrdiff_t> node_position;                        // Relative distances of the nodes to the defect position
  
  // Contributions from the borders
  std::vector <std::size_t>  border_element1;                        // Position of broken deffects in hopping terms                     
  std::vector <std::size_t>  border_element2;                        // Position of broken deffects in hopping terms                     
  std::vector <T>             border_hopping;
  std::vector<std::vector <value_type>>  border_v;                        // temporary generalized velocities
  std::vector <T>                   border_U;                                                 
  std::vector <std::size_t>   border_element;                         // Position of broken deffects site energy
  Hamiltonian<T,D>                       & h;
  LatticeStructure <D>                   & r;
  GLOBAL_VARIABLES <T>              & Global;
  std::vector <std::vector<std::size_t>> position;                   // vector of vectors with positions in the lattice of the Orbital 0 of the defects for each stride block
  Eigen::Array<T, -1, -1>        new_hopping;
  KPMRandom <T>                        & rnd; 
  
  Defect_Operator(Hamiltonian<T,D> &,  std::string & defect, H5::H5File *file);
  void generate_disorder();
  template <unsigned MULT, bool VELOCITY>
  void multiply_defect(std::size_t istr, T* & phi0, T* & phiM1, unsigned axis);
  template <unsigned MULT, bool VELOCITY>
  void multiply_broken_defect(T* & phi0, T* & phiM1, unsigned axis);
  void build_velocity(std::vector<unsigned> & components, unsigned n);

  void interface_multiply_defect(std::size_t istr, T* & phi0, T* & phiM1, unsigned axis);
  void interface_multiply_broken_defect(T* & phi0, T* & phiM1, unsigned axis);  
};
