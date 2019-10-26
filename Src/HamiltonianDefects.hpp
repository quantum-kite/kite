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
  std::vector <std::vector<std::size_t>> position;                   // vector of vectors with positions in the lattice of the Orbital 0 of the defects for each tile block
  Eigen::Array<T, -1, -1>        new_hopping;
  KPMRandom <T>                        & rnd;
  std::vector <int>          positions_fixed;
  
  Defect_Operator(Hamiltonian<T,D> &,  std::string & defect, H5::H5File *file);
  void generate_disorder();
  void build_velocity(std::vector<unsigned> & components, unsigned n);


  template <unsigned MULT, bool VELOCITY>  
  void multiply_defect(std::size_t istr, T* & phi0, T* & phiM1, unsigned axis) {
    Coordinates<std::ptrdiff_t, D + 1>  local1(r.Ld);
    for(std::size_t i = 0; i <  position.at(istr).size(); i++)
      {
        std::size_t ip = position.at(istr)[i];
        std::size_t iv = local1.set_coord(ip).coord[D - 1];
        for(unsigned k = 0; k < hopping.size(); k++)
          {
            std::size_t k1 = ip + node_position[element1[k]];
            std::size_t k2 = ip + node_position[element2[k]];
            
            if(VELOCITY)
              phi0[k1] += value_type(MULT + 1) * v.at(axis).at(k)*new_hopping(k, iv) * phiM1[k2] ;
            else
              phi0[k1] += value_type(MULT + 1) * new_hopping(k, iv) * phiM1[k2] ;
          }
        
        if(!VELOCITY)
          for(std::size_t k = 0; k < U.size(); k++)
            {
              std::size_t k1 = ip + node_position[element[k]];
              phi0[k1] += value_type(MULT + 1) * U[k] * phiM1[k1];
            }
      }
  }

  
  template <unsigned MULT, bool VELOCITY>
  void multiply_broken_defect(T* & phi0, T* & phiM1, unsigned axis) {
    Coordinates<std::ptrdiff_t, D + 1> global1(r.Lt), global2(r.Lt), local1(r.Ld) ;
    Eigen::Map<Eigen::Matrix<std::ptrdiff_t,2,1>> v_global1(global1.coord), v_global2(global2.coord);
    double phase;
    
    Eigen::Matrix<double, 1, 2> temp_vect;
    for(std::size_t i = 0; i < border_element1.size(); i++)
      {
        std::size_t i1 = border_element1[i];
        std::size_t i2 = border_element2[i];
        
        // These four lines pertrain only to the ghost_correlation
        r.convertCoordinates(global1, local1.set_coord(i1));
        r.convertCoordinates(global2, local1.set_coord(i2));
        temp_vect  = (v_global2 - v_global1).template cast<double>().matrix().transpose();
        phase = temp_vect(0)*r.ghost_pot(0,1)*v_global1(1); //.template cast<double>().matrix();
        
        if(VELOCITY)
          phi0[i1] += value_type(MULT + 1) * border_v.at(axis).at(i) * border_hopping[i] * phiM1[i2] * multEiphase(phase);
        else
          phi0[i1] += value_type(MULT + 1) * border_hopping[i] * phiM1[i2] * multEiphase(phase);
      }
    
    if(!VELOCITY)
      for(std::size_t i = 0; i < border_element.size(); i++)
        {
          std::size_t i1 = border_element[i];
          phi0[i1] += value_type(MULT + 1) * border_U[i] * phiM1[i1];
        }
  }
};
