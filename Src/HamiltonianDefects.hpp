
template <typename T,unsigned D>
struct Defect_Operator  {
  double                                   p;                        // Concentration of defects
  unsigned                       NumberNodes;                        // Number of nodes in the deffect
  std::vector <T>                          U;                        // local energies
  std::vector <unsigned>             element;                        // nodes with local energies
  std::vector <T>                    hopping;                        // vector of the non-zero values of the operator
  std::vector <T>           hopping_magnetic;                        
  std::vector <T>                       V[D];
  std::vector <T>                   V2[D][D];
  std::vector <unsigned>            element1;                        // vector with the nodes 
  std::vector <unsigned>            element2;                        // vector with the nodes
  std::vector <std::ptrdiff_t> node_position;                        // Relative distances of the nodes to the defect position
  
  // Contributions from the borders
  std::vector <std::size_t>  border_element1;                        // Position of broken deffects in hopping terms                     
  std::vector <std::size_t>  border_element2;                        // Position of broken deffects in hopping terms                     
  std::vector <T>             border_hopping;
  std::vector <T>    border_hopping_magnetic;
  std::vector <T>                border_V[D];                                           
  std::vector <T>            border_V2[D][D];   
  std::vector <T>                   border_U;                                                 
  std::vector <std::size_t>   border_element;                        // Position of broken deffects site energy                                   
  LatticeStructure <D>                   & r;
  Simulation <T,D>                   & simul;
  std::vector <std::vector<std::size_t>> position;                   // vector of vectors with positions in the lattice of the Orbital 0 of the defects for each stride block
  
  Defect_Operator(Simulation <T,D> & sim, std::string & defect, H5::H5File *file) : r(sim.r), simul(sim), position(sim.r.NStr) {
	debug_message("Entered Defect_Operator\n");

    std::string field = defect + std::string("/Concentration");
    get_hdf5<double> ( &p, file, field );

    /* Read Number of nodes and their relative  positions */
    field = defect + std::string("/NumNodes");
    get_hdf5<unsigned> ( &NumberNodes, file, field );
    {
      std::vector<unsigned> tmp(NumberNodes);
      field = defect + std::string("/NodePosition");
      get_hdf5<unsigned> ( tmp.data(), file, field );
      node_position.resize(NumberNodes);
      for(unsigned i = 0; i < NumberNodes; i++)
	node_position.at(i) = tmp.at(i);
    }
    
    /* Read Hoppings */

    {
      int n;
      std::vector<int> tmp;
      field = defect + std::string("/NumBondDisorder");
      get_hdf5<int> ( &n, file, field );
      
      tmp.resize(n);
      hopping.resize(n);
      hopping_magnetic.resize(n);
      element1.resize(n);
      element2.resize(n);
      
      field = defect + std::string("/NodeTo");
      get_hdf5<int> (tmp.data(), file, field );
      for(int i = 0; i < n; i++)
	element1.at(i) = tmp.at(i);
      
      field = defect + std::string("/NodeFrom");
      get_hdf5<int> (tmp.data(), file, field );
      for(int i = 0; i < n; i++)
	element2.at(i) = tmp.at(i);
      
      field = defect + std::string("/Hopping");
      get_hdf5<T> (hopping.data(), file, field );
    }
    
    /* Read local Disorder */
    {
      int n;
      std::vector<int> tmp;
      field = defect + std::string("/NumOnsiteDisorder");
      get_hdf5<int> ( &n, file, field );
      tmp.resize(n);     
      U.resize(n);
      element.resize(n);
      
      field = defect + std::string("/NodeOnsite");
      get_hdf5<int> (tmp.data(), file, field );
      for(int i = 0; i < n; i++)
	element.at(i) = tmp.at(i);
      
      field = defect + std::string("/U0");
      get_hdf5<T> (U.data(), file, field );
    }
    
    /* Translate Positions */
    unsigned l[D + 1];
    std::fill_n(l, D, 3);
    l[D]  = r.Orb;
    
    Coordinates<std::ptrdiff_t, D + 1> b3(l), Ld(r.Ld);    
    Eigen::Map<Eigen::Matrix<std::ptrdiff_t,D, 1>> v(b3.coord);
    
    for(unsigned node = 0; node < node_position.size(); node++)
      {
	b3.set_coord(node_position[node]);                             // Get Coordinates in Basis 3,  The last is the Final Orbital.
	v.array() -= 1; 	                                       // Subtract to the first D elements of v to get v(i) in (-1, 0 , 1)
	node_position[node] = Ld.set_index(b3.coord).index;            // Convert in distances in this lattice
      }
    
    /* Build Velocity */ 

    Coordinates<std::ptrdiff_t, D + 1> Lda(r.Ld), Ldb(r.Ld);
    Eigen::Map<Eigen::Matrix<std::ptrdiff_t,D, 1>> va(Lda.coord), vb(Ldb.coord); // Column vector
    Eigen::Matrix<double, D, 1> dr_R;
    Eigen::Matrix<double, D, 1> dr_a;
    Eigen::Matrix<double, D, 1> orbital_difference_R;
    Eigen::Matrix<double, D, 1> lattice_difference_R;
    Eigen::Matrix<double, D, 1> orbital_difference_a;
    Eigen::Matrix<double, D, 1> lattice_difference_a;


    for(unsigned ih = 0; ih < hopping.size(); ih++)
      {
	
	Lda.set_coord(static_cast<std::ptrdiff_t>(r.Nd/2 + node_position[element1.at(ih)]));
	Ldb.set_coord(static_cast<std::ptrdiff_t>(r.Nd/2 + node_position[element2.at(ih)]));

	// difference vectors in real-space coordinates
	orbital_difference_R = r.rOrb.col(Ldb.coord[D]) - r.rOrb.col(Lda.coord[D]);
	lattice_difference_R = r.rLat * (vb - va).template cast<double>();
	dr_R = orbital_difference_R + lattice_difference_R;
	
	// difference vectors in lattice coordinates
	orbital_difference_a = r.rLat.inverse() * orbital_difference_R;
	lattice_difference_a = r.rLat.inverse() * lattice_difference_R;
	dr_a = orbital_difference_a + lattice_difference_a;
	
	
	// periodic part of the Peierls phase. 
	// If the kpm vector is not complex, peierls1(phase) will return 1.0, so it is effectively harmless.
	//double phase = ((0.5*dr_a.transpose() + (r.rLat.inverse()*r.rOrb.col(Lda.coord[D])).transpose())*r.vect_pot*dr_a - lattice_difference_a.transpose()*r.vect_pot*r.rLat.inverse()*r.rOrb.col(Ldb.coord[D]))(0,0);
	double phase = (-(-0.5*dr_a.transpose() + (r.rLat.inverse()*r.rOrb.col(Ldb.coord[D])).transpose())*r.vect_pot*dr_a + lattice_difference_a.transpose()*r.vect_pot*r.rLat.inverse()*r.rOrb.col(Lda.coord[D]))(0,0);
	hopping_magnetic.at(ih) = hopping.at(ih) * peierls1(phase);
	hopping.at(ih) = hopping_magnetic.at(ih);
	
	
	for(unsigned dim = 0; dim < D; dim++)
	  V[dim].push_back( hopping.at(ih) * T(dr_R(dim)) );

	for(unsigned dim1 = 0; dim1 < D; dim1++)
	  for(unsigned dim2 = 0; dim2 < D; dim2++)
	    V2[dim1][dim2].push_back( hopping.at(ih) * T(dr_R(dim1)) * T(dr_R(dim2)));
      }
      
    debug_message("Left Defect_Operator constructor.\n");
  };



  
  void generate_disorder()  {
	  debug_message("Entered generate_disorder\n");
    /* Structural disorder*/

    /*
     *   Empty the positions borders of the defects 
     *  
     */

    border_element1.clear();                                           
    border_element2.clear();                                           
    border_hopping.clear();                                           
    border_U.clear();                                                 
    border_element.clear();
    for(std::size_t istr = 0; istr < r.NStr; istr++)
      position.at(istr).clear();
    
#pragma omp master
    {
      simul.Global.element1.clear();
      simul.Global.element2_diff.clear();
      simul.Global.hopping.clear();
      simul.Global.element.clear();
      simul.Global.U.clear();
    }
#pragma omp barrier

    Coordinates<std::size_t,D + 1> latt(r.ld), LATT(r.Lt), Latt(r.Ld), latStr(r.lStr);
    // Distribute the local disorder

    std::size_t ndefects= p * r.N , count = 0;

    while(count < ndefects)
      {
	std::size_t pos = r.N * simul.rnd.get();
	latt.set_coord(pos);
	r.convertCoordinates(Latt,latt);
	r.convertCoordinates(latStr,latt);
	
	if( !any_of(position.at(latStr.index).begin(), position.at(latStr.index).end(), std::bind2nd(std::equal_to<std::size_t>(), Latt.index)))
	  {	    
	    position.at(latStr.index).push_back(Latt.index);
	    count++;
	  }
      }    
    // Test if any of the defect cross the borders
    for(std::size_t istr = 0; istr < r.NStr; istr++)
      for(auto it = position.at(istr).begin(); it != position.at(istr).end(); it++)
	for(unsigned node = 0; node < NumberNodes; node++ )
	  {
	    auto node_pos =  *it + node_position.at(node);
	    Latt.set_coord(node_pos);                             // coordinates of the node in Ld Lattice
	    
	    if(r.test_ghosts(Latt) == 0)
	      {
		r.convertCoordinates(LATT, Latt);
#pragma omp critical
		{
		  for(unsigned i = 0; i < element1.size(); i++)
		    if(node == element1[i])
		      simul.Global.addbond(LATT.index,  node_position.at(element2[i]) - node_position.at(element1[i]), hopping[i]);
		  
		  for(unsigned i = 0; i < element.size(); i++)
		    if(node == element[i])
		      simul.Global.addlocal(LATT.index, U[i]);
		}
	      }
	  }
#pragma omp barrier
    /* Look for the extra bonds in this domain */


#pragma omp critical
    {
      Coordinates<std::ptrdiff_t, D + 1> latt(r.ld), LATT(r.Lt), Latt(r.Ld), Ldb(r.Ld);
      Eigen::Map<Eigen::Matrix<std::ptrdiff_t,D, 1>> va(Latt.coord), vb(Ldb.coord); // Column vector
      
      
		  
	  Eigen::Matrix<double, D, 1> dr_R;
	  Eigen::Matrix<double, D, 1> dr_a;
	  Eigen::Matrix<double, D, 1> orbital_difference_R;
	  Eigen::Matrix<double, D, 1> lattice_difference_R;
	  Eigen::Matrix<double, D, 1> orbital_difference_a;
	  Eigen::Matrix<double, D, 1> lattice_difference_a;
      
      for(unsigned i = 0; i < simul.Global.element1.size(); i++ )
	if(r.domain_number( long(simul.Global.element1[i]) ) == std::ptrdiff_t(r.thread_id))
	  {
	    LATT.set_coord(simul.Global.element1[i]);
	    r.convertCoordinates(Latt, LATT);
	    border_element1.push_back( Latt.index );	    
	    border_element2.push_back(Latt.index + simul.Global.element2_diff[i]);
	    border_hopping.push_back(simul.Global.hopping[i]);
	    Ldb.set_coord(static_cast<std::ptrdiff_t>(Latt.index + simul.Global.element2_diff[i] ));
		
		// difference vectors in real-space coordinates
		orbital_difference_R = r.rOrb.col(Ldb.coord[D]) - r.rOrb.col(Latt.coord[D]);
		lattice_difference_R = r.rLat * (vb - va).template cast<double>();
		dr_R = orbital_difference_R + lattice_difference_R;
		
        
	    
	    for(unsigned dim = 0; dim < D; dim++)
	      border_V[dim].push_back(simul.Global.hopping[i] * T(dr_R(dim)));

	    for(unsigned dim1 = 0; dim1 < D; dim1++)
	      for(unsigned dim2 = 0; dim2 < D; dim2++)
	        border_V2[dim1][dim2].push_back(simul.Global.hopping[i] * T(dr_R(dim1))* T(dr_R(dim2)));
	  }
      
      for(unsigned i = 0; i < simul.Global.element.size(); i++ )
	if(r.domain_number (std::ptrdiff_t(simul.Global.element[i])) == std::ptrdiff_t(r.thread_id))
	  {
	    LATT.set_coord(simul.Global.element[i] );
	    r.convertCoordinates(Latt, LATT );
	    border_element.push_back(Latt.index );	    
	    border_U.push_back(simul.Global.U[i] );
	  }
    }
#pragma omp barrier



    
    /*
      The multiplication is done through a mozaic structure where each tile is initialized
      in the beginning before the hoppings of the regular hamiltonian.
      The defects are added after the regular Hamiltonian and could involve hoppings between elements 
      of different tiles. Because the deffect structure involve the application of the disorder hamiltonian
      to all the nodes belonging to the deffect, for the deffects that cross the borders to a higher tile, we need to garantee that 
      that tile  was already initialized before the lattice run and will not be initialized again in the future.
    */
    
    /* Test Mozaic to implement in tile that have to be set to zero  */
    
    for(std::size_t istr = 0; istr < r.NStr; istr++)
      for(auto it = position.at(istr).begin(); it != position.at(istr).end(); it++)
	for(unsigned node = 0; node < NumberNodes; node++ )
	  {
	    std::size_t node_pos =  *it + node_position.at(node);
	    Latt.set_coord(node_pos);
	    r.convertCoordinates(latStr, Latt ); // Get stride index
	    // Tests if the node is in a higher stride,
	    
	    if(r.test_ghosts(Latt) == 1 && simul.h.cross_mozaic[latStr.index] &&  latStr.index > istr )
	      {	    
		simul.h.cross_mozaic[latStr.index] = false;
		simul.h.cross_mozaic_indexes.push_back(latStr.index); // Add  because
	      }
	  }
    
    for(std::size_t istr = 0; istr < r.NStr; istr++)
      std::sort(position[istr].begin(), position[istr].end() );
    
    
   debug_message("Left generate_disorder\n");
  }
  


    template <typename U = T>
  typename std::enable_if<is_tt<std::complex, U>::value, U>::type peierls1(double phase) {
	std::complex<double> im(0,1.0);
    return U(exp(im*phase));
  };
  
  template <typename U = T>
  typename std::enable_if<!is_tt<std::complex, U>::value, U>::type peierls1(double phase) {
    return 1.0;
  };

  
};
