
template <typename T,unsigned D>
struct Defect_Operator  {
  typedef typename extract_value_type<T>::value_type value_type;
  double                                   p;                        // Concentration of defects
  unsigned                       NumberNodes;                        // Number of nodes in the deffect
  std::vector <T>                          U;                        // local energies
  std::vector <unsigned>             element;                        // nodes with local energies
  std::vector <T>                    hopping;                        // vector of the non-zero values of the hopping operator operator                // 
  std::vector<std::vector <value_type>>         v;                        // temporary generalized velocities
  std::vector <unsigned>            element1;                        // vector with the nodes 
  std::vector <unsigned>            element2;                        // vector with the nodes
  std::vector <std::ptrdiff_t> node_position;                        // Relative distances of the nodes to the defect position
  
  // Contributions from the borders
  std::vector <std::size_t>  border_element1;                        // Position of broken deffects in hopping terms                     
  std::vector <std::size_t>  border_element2;                        // Position of broken deffects in hopping terms                     
  std::vector <T>             border_hopping;
  std::vector<std::vector <value_type>>  border_v;                        // temporary generalized velocities
  std::vector <T>                   border_U;                                                 
  std::vector <std::size_t>   border_element;                        // Position of broken deffects site energy                                   
  LatticeStructure <D>                   & r;
  Simulation <T,D>                   & simul;
  std::vector <std::vector<std::size_t>> position;                   // vector of vectors with positions in the lattice of the Orbital 0 of the defects for each stride block
  Eigen::Array<T, -1, -1>        new_hopping;

  
  Defect_Operator(Simulation <T,D> & sim, std::string & defect, H5::H5File *file) : r(sim.r), simul(sim), position(sim.r.NStr)
  {
    debug_message("Entered Defect_Operator\n");
    
    std::string field = defect + std::string("/Concentration");
    get_hdf5<double> ( &p, file, field );
    std::vector<unsigned> tmpU;
    std::vector<int> tmpI;
    int n;
    
    /* Read Number of nodes and their relative  positions */
    field = defect + std::string("/NumNodes");
    get_hdf5<unsigned> ( &NumberNodes, file, field );
    tmpU.resize(NumberNodes);
    field = defect + std::string("/NodePosition");
    get_hdf5<unsigned> ( tmpU.data(), file, field );
    node_position.resize(NumberNodes);
    for(unsigned i = 0; i < NumberNodes; i++)
      node_position.at(i) = tmpU.at(i);
    
    
    /* Read Hoppings */
    
    field = defect + std::string("/NumBondDisorder");
    get_hdf5<int> ( &n, file, field );
    
    tmpI.resize(n);
    hopping.resize(n);
    element1.resize(n);
    element2.resize(n);
    
    field = defect + std::string("/NodeTo");
    get_hdf5<int> (tmpI.data(), file, field );
    for(int i = 0; i < n; i++)
      element1.at(i) = tmpI.at(i);
    
    field = defect + std::string("/NodeFrom");
    get_hdf5<int> (tmpI.data(), file, field );
    for(int i = 0; i < n; i++)
      element2.at(i) = tmpI.at(i);
    
    field = defect + std::string("/Hopping");
    get_hdf5<T> (hopping.data(), file, field );
    
    
    /* Read local Disorder */
    
    field = defect + std::string("/NumOnsiteDisorder");
    get_hdf5<int> ( &n, file, field );
    tmpI.resize(n);     
    U.resize(n);
    element.resize(n);
    
    field = defect + std::string("/NodeOnsite");
    get_hdf5<int> (tmpI.data(), file, field );
    for(int i = 0; i < n; i++)
      element.at(i) = tmpI.at(i);
    
    field = defect + std::string("/U0");
    get_hdf5<T> (U.data(), file, field );
    
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
    
    Coordinates<std::ptrdiff_t, D + 1> Lda(r.Ld), Ldb(r.Ld),  Lta(r.Lt);
    Eigen::Map<Eigen::Matrix<std::ptrdiff_t,D, 1>> va(Lda.coord), vb(Ldb.coord), Vta(Lta.coord) ; // Column vector
    Eigen::Matrix<double, D, 1> dif_R, dif_O, sum_O, ra;
    Eigen::Matrix<double, D, D> matA = r.rLat.inverse().transpose() * r.vect_pot * r.rLat.inverse();
    new_hopping = Eigen::Matrix<T, Eigen::Dynamic,Eigen::Dynamic>::Zero(hopping.size(), r.Ld[D - 1]);
    
    
    for(unsigned ih = 0; ih < hopping.size(); ih++)
      {
	double phase1 = 0., phase2 = 0., phase3 = 0.;
	for(std::size_t iv = NGHOSTS; iv < r.Ld[D - 1] - NGHOSTS; iv++)
	  {
	    std::size_t ip = NGHOSTS + iv *Lda.basis[D - 1];
	    Lda.set_coord(static_cast<std::ptrdiff_t>( ip + node_position[element1.at(ih)]));
	    Ldb.set_coord(static_cast<std::ptrdiff_t>( ip + node_position[element2.at(ih)]));
	    dif_R = r.rLat * (va - vb).template cast<double>();
	    dif_O = r.rOrb.col(Lda.coord[D]) - r.rOrb.col(Ldb.coord[D]);
	    sum_O = r.rOrb.col(Ldb.coord[D]) + r.rOrb.col(Lda.coord[D]);
	    r.convertCoordinates(Lta, Lda);
	    ra = r.rLat * Vta.template cast<double>();
	
	    phase1 =  0.5 * (dif_R + sum_O).transpose() * (matA * (dif_R + dif_O)).matrix();
	    phase2 =  - dif_R.transpose() * (matA * r.rOrb.col(Lda.coord[D])).matrix();
	    phase3 =  - dif_R.transpose() * (matA * ra).matrix();	   
	    new_hopping(ih, iv) = hopping.at(ih) * simul.h.peierls(phase1 + phase2 + phase3);
	  }
	hopping.at(ih) *= simul.h.peierls(phase1 + phase2);
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

    Coordinates<std::size_t,D + 1> latt(r.ld), LATT(r.Lt), Latt(r.Ld), Latt2(r.Ld), latStr(r.lStr);
    // Distribute the local disorder

    std::size_t ndefects= p * r.N , count = 0;
    
    while(count < ndefects) 
      {
	std::size_t  pos = r.N * simul.rnd.get();
	latt.set_coord(pos);
	r.convertCoordinates(Latt,latt);
	r.convertCoordinates(latStr,latt);
	auto & st = position.at(latStr.index);

	if( !any_of(st.begin(), st.end(), std::bind2nd(std::equal_to<std::size_t>(), Latt.index)))
	  {	    
	    st.push_back(Latt.index);
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

	    if(r.test_ghosts(Latt) == 1)                               // test if it is in the ghosts
	      {
		/*
		  For the nodes of the defects inside the sample
		  I need to test if they will add amplitudes for 
		  nodes of already visited and emptied tiles
		 */

		r.convertCoordinates(latStr, Latt);
		if(latStr.index < istr)                           
		  simul.h.hV.add_conflict_with_defect(std::size_t(node_pos), latStr.index);
	      }
	    else
	      {
		// node is in the ghosts
		r.convertCoordinates(LATT, Latt);
#pragma omp critical
		{
		  /*
		    For the nodes of the defects outside the sample
		    I need to add the local terms and bonds of the Hamiltonian 
		    do the global memory to be passed to the neighbour domain

		    In the case of bonds, I will test if the pair corresponds to 
		    a vacancy in this domain. In this case, the bond will not be added.
		    
		  */
		  
		  for(unsigned i = 0; i < element1.size(); i++)
		    if(node == element1[i])
		      {
			std::ptrdiff_t  dd =  node_position.at(element2[i]) - node_position.at(element1[i]);
			auto node2_pos =  *it + node_position.at(element2[i]);
			Latt2.set_coord(node2_pos);
			
			if(r.test_ghosts(Latt2) == 1)
			  {
			    // Not in the ghosts
			    if( simul.h.hV.test_vacancy(Latt2) == 0)  // Not a Vacancy
			      simul.Global.addbond(LATT.index,  dd, hopping[i]);
			  }
			else
			  simul.Global.addbond(LATT.index,  dd, hopping[i]); // If both nodes are the ghosts I add it
		      }
		  
		  
		  for(unsigned i = 0; i < element.size(); i++)
		    if(node == element[i])
		      simul.Global.addlocal(LATT.index, U[i]);
		}
	      }
	  }
#pragma omp barrier
    /* 
       Look for the extra bonds in this domain 
       That come from broken defects in the boundaries 
    */
    
    
#pragma omp critical
    {
      Coordinates<std::ptrdiff_t, D + 1> latt(r.ld), LATT(r.Lt), Latt(r.Ld), Ldb(r.Ld), latStr(r.lStr);
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
	    r.convertCoordinates(latStr, Latt);
	    
	    auto & st = simul.h.hV.position.at(latStr.index); 
	    if( !any_of(st.begin(), st.end(), std::bind2nd(std::equal_to<std::size_t>(), Latt.index)))
	      {     
		border_element1.push_back( Latt.index );	    
		border_element2.push_back(Latt.index + simul.Global.element2_diff[i]);
		border_hopping.push_back(simul.Global.hopping[i]);
		
	      }
	  }
      
      for(unsigned i = 0; i < simul.Global.element.size(); i++ )
	if(r.domain_number (std::ptrdiff_t(simul.Global.element[i])) == std::ptrdiff_t(r.thread_id))
	  {
	    LATT.set_coord(simul.Global.element[i] );
	    r.convertCoordinates(Latt, LATT );
	    r.convertCoordinates(latStr, Latt);
	    auto & st = simul.h.hV.position.at(latStr.index); 
	    if( !any_of(st.begin(), st.end(), std::bind2nd(std::equal_to<std::size_t>(), Latt.index)))
	      {
		border_element.push_back(Latt.index );	    
		border_U.push_back(simul.Global.U[i] );
	      }
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
  
  template <unsigned MULT, bool VELOCITY>
  void multiply_defect(std::size_t istr, T* & phi0, T* & phiM1, unsigned axis)
  {
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
  void multiply_broken_defect(T* & phi0, T* & phiM1, unsigned axis)
  {
    Coordinates<std::ptrdiff_t, D + 1> global1(r.Lt), global2(r.Lt), local1(r.Ld) ;
    Eigen::Map<Eigen::Matrix<std::ptrdiff_t,2,1>> v_global1(global1.coord), v_global2(global2.coord);
    double phase;

    Eigen::Matrix<double, 1, 2> temp_vect;
    for(std::size_t i = 0; i < border_element1.size(); i++)
      {
	std::size_t i1 = border_element1[i];
	std::size_t i2 = border_element2[i];
	
	// These four lines pertrain only to the magnetic field
	r.convertCoordinates(global1, local1.set_coord(i1));
	r.convertCoordinates(global2, local1.set_coord(i2));
	temp_vect  = (v_global2 - v_global1).template cast<double>().matrix().transpose();
	phase = temp_vect(0)*r.vect_pot(0,1)*v_global1(1); //.template cast<double>().matrix();
	
	if(VELOCITY)
	  phi0[i1] += value_type(MULT + 1) * border_v.at(axis).at(i) * border_hopping[i] * phiM1[i2] * simul.h.peierls(phase);
	else
	phi0[i1] += value_type(MULT + 1) * border_hopping[i] * phiM1[i2] * simul.h.peierls(phase);
      }
    
    if(!VELOCITY)
    for(std::size_t i = 0; i < border_element.size(); i++)
      {
	std::size_t i1 = border_element[i];
	phi0[i1] += value_type(MULT + 1) * border_U[i] * phiM1[i1];
      }
  }
  void build_velocity(std::vector<unsigned> & components, unsigned n)
  {
    Coordinates<std::ptrdiff_t, D + 1> Lda(r.Ld), Ldb(r.Ld);
    Eigen::Map<Eigen::Matrix<std::ptrdiff_t,D, 1>> va(Lda.coord), vb(Ldb.coord); // Column vector
    Eigen::Matrix<double, D, 1> orbital_difference_R;
    Eigen::Matrix<double, D, 1> lattice_difference_R;
    Eigen::Matrix<double, D, 1> dr_R;
    if ( n == v.size())
      v.push_back(std::vector<value_type>(hopping.size()));
    if ( n == border_v.size())
      border_v.push_back(std::vector<value_type>(border_hopping.size()));
    if(n > v.size())
      std::cout << "Simao esta a fazer asneira" << std::endl;
    
    std::vector<value_type> & v1 = v.at(n);
    std::vector<value_type> & border_v1 = border_v.at(n);  

    std::ptrdiff_t ip = 0;
    for(unsigned i = 0; i < D; i++)
      ip += simul.r.Ld[i]/2 * Lda.basis[i];
    Lda.set_coord(ip);
    
    for(unsigned ih = 0; ih < hopping.size(); ih++)
      {
	
	Lda.set_coord(static_cast<std::ptrdiff_t>(ip + node_position[element1.at(ih)]));
	Ldb.set_coord(static_cast<std::ptrdiff_t>(ip + node_position[element2.at(ih)]));
	lattice_difference_R = r.rLat * (vb - va).template cast<double>();
	orbital_difference_R = r.rOrb.col(Ldb.coord[D]) - r.rOrb.col(Lda.coord[D]);
	dr_R = orbital_difference_R + lattice_difference_R;
	v1.at(ih) = value_type(1);
	for(unsigned i = 0; i < components.size(); i++)
	  v1.at(ih) *= value_type(dr_R(components.at(i)));
      }
	

  
    
    for(unsigned ih = 0; ih < border_hopping.size(); ih++)
      {
	Lda.set_coord(static_cast<std::ptrdiff_t>(border_element1.at(ih)));
	Ldb.set_coord(static_cast<std::ptrdiff_t>(border_element2.at(ih)));
	lattice_difference_R = r.rLat * (vb - va).template cast<double>();
	orbital_difference_R = r.rOrb.col(Ldb.coord[D]) - r.rOrb.col(Lda.coord[D]);
	dr_R = orbital_difference_R + lattice_difference_R;
	border_v1.at(ih) = value_type(1);
	for(unsigned i = 0; i < components.size(); i++)
	  border_v1.at(ih) *= value_type(dr_R(components.at(i)));
      }
  }
  
};
