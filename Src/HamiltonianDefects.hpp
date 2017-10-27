
template <typename T,unsigned D>
struct Defect_Operator  {
  double                                   p;                        // Concentration of defects
  unsigned                       NumberNodes;                        // Number of nodes in the deffect
  std::vector <T>                          U;                        // local energies
  std::vector <unsigned>             element;                        // nodes with local energies
  std::vector <T>                    hopping;                        // vector of the non-zero values of the operator
  std::vector <unsigned>            element1;                        // vector with the nodes 
  std::vector <unsigned>            element2;                        // vector with the nodes
  std::vector <std::ptrdiff_t> node_position;                        // Relative distances of the nodes to the defect position
  
  // Contributions from the borders
  std::vector <std::size_t>  border_element1;                        // Position of broken deffects in hopping terms                     
  std::vector <std::size_t>  border_element2;                        // Position of broken deffects in hopping terms                     
  std::vector <T>             border_hopping;                                           
  std::vector <T>                   border_U;                                                 
  std::vector <std::size_t>   border_element;                        // Position of broken deffects site energy                                   
  LatticeStructure <D>                   & r;
  Simulation <T,D>                   & simul;
  std::vector <std::vector<std::size_t>> position;                   // vector of vectors with positions in the lattice of the Orbital 0 of the defects for each stride block
  std::vector <bool> cross_mozaic;
  std::vector <std::size_t> cross_mozaic_indexes;
  
  Defect_Operator(Simulation <T,D> & sim) : r(sim.r), simul(sim), position(sim.r.NStr), cross_mozaic(sim.r.NStr) {
    T timp, nu, rho, u0;
    
#pragma omp critical
    {
      H5::H5File    *file      = new H5::H5File(simul.name, H5F_ACC_RDONLY);
      get_hdf5<double> ( &p, file, (char *) "/Concentration");
      get_hdf5<T>      (&nu, file, (char *) "/NU"  );
      get_hdf5<T>    (&timp, file, (char *) "/TIMP");
      get_hdf5<T>    ( &rho, file, (char *) "/RHO" );
      get_hdf5<T>    (  &u0, file, (char *) "/U0"  );
      file->close();	    
    }

    NumberNodes = 6;
    node_position.resize(6);
    node_position.at(0) = (1 + 0) + (1 + 0)*3 + 0 * 9 ;                                            //
    node_position.at(1) = (1 + 0) + (1 + 0)*3 + 1 * 9;                  //
    node_position.at(2) = (1 + 1) + (1 + 0)*3 + 0 * 9;                  //      4  3
    node_position.at(3) = (1 + 0) + (1 + 1)*3 + 1 * 9;                  //    5      2
    node_position.at(4) = (1 + 0) + (1 + 1)*3 + 0 * 9;                  //      0  1
    node_position.at(5) = (1 - 1) + (1 + 1)*3 + 1 * 9;                  //

    hopping.resize(30);
    element1.resize(30);
    element2.resize(30);
    
    std::size_t index = 0;

    //   timp  
    element1[index] = 0; element2[index] = 1; hopping[index] = timp; index++;
    element1[index] = 1; element2[index] = 2; hopping[index] = timp; index++;
    element1[index] = 2; element2[index] = 3; hopping[index] = timp; index++;
    element1[index] = 3; element2[index] = 4; hopping[index] = timp; index++;
    element1[index] = 4; element2[index] = 5; hopping[index] = timp; index++;
    element1[index] = 5; element2[index] = 0; hopping[index] = timp; index++;
    
    element1[index] = 1; element2[index] = 0; hopping[index] = std::conj(timp); index++;
    element1[index] = 2; element2[index] = 1; hopping[index] = std::conj(timp); index++;
    element1[index] = 3; element2[index] = 2; hopping[index] = std::conj(timp); index++;
    element1[index] = 4; element2[index] = 3; hopping[index] = std::conj(timp); index++;
    element1[index] = 5; element2[index] = 4; hopping[index] = std::conj(timp); index++;
    element1[index] = 0; element2[index] = 5; hopping[index] = std::conj(timp); index++;

    
    //    nu    
    element1[index] = 0; element2[index] = 2; hopping[index] = nu; index++;
    element1[index] = 2; element2[index] = 4; hopping[index] = nu; index++;
    element1[index] = 4; element2[index] = 0; hopping[index] = nu; index++;
    element1[index] = 1; element2[index] = 3; hopping[index] = nu; index++;
    element1[index] = 3; element2[index] = 5; hopping[index] = nu; index++;
    element1[index] = 5; element2[index] = 1; hopping[index] = nu; index++;

    element1[index] = 2; element2[index] = 0; hopping[index] = std::conj(nu); index++;
    element1[index] = 4; element2[index] = 2; hopping[index] = std::conj(nu); index++;
    element1[index] = 0; element2[index] = 4; hopping[index] = std::conj(nu); index++;
    element1[index] = 3; element2[index] = 1; hopping[index] = std::conj(nu); index++;
    element1[index] = 5; element2[index] = 3; hopping[index] = std::conj(nu); index++;
    element1[index] = 1; element2[index] = 5; hopping[index] = std::conj(nu); index++;
    
    //  rho 

    element1[index] = 0; element2[index] = 3; hopping[index] = rho; index++;
    element1[index] = 1; element2[index] = 4; hopping[index] = rho; index++;
    element1[index] = 2; element2[index] = 5; hopping[index] = rho; index++;
    element1[index] = 3; element2[index] = 0; hopping[index] = std::conj(rho); index++;
    element1[index] = 4; element2[index] = 1; hopping[index] = std::conj(rho); index++;
    element1[index] = 5; element2[index] = 2; hopping[index] = std::conj(rho); index++;
    
    /* U */ 
    
    element.resize(6);
    U.resize(6);
    for(int i = 0; i < 6; i++)
      {
	element[i] = i;
	U[i] = u0;
      }
    
    
    
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


  };



  
  void generate_disorder()  {
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
      {
	position.at(istr).clear();
	cross_mozaic[istr] = true;
      }
    cross_mozaic_indexes.clear();
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
      Coordinates<std::ptrdiff_t,D + 1> latt(r.ld), LATT(r.Lt), Latt(r.Ld);
      for(unsigned i = 0; i < simul.Global.element1.size(); i++ )
	if(r.domain_number( long(simul.Global.element1[i]) ) == std::ptrdiff_t(r.thread_id))
	  {
	    LATT.set_coord(simul.Global.element1[i]);
	    r.convertCoordinates(Latt, LATT);
	    border_element1.push_back( Latt.index );	    
	    border_element2.push_back( Latt.index + simul.Global.element2_diff[i]);
	    border_hopping.push_back(simul.Global.hopping[i]);  
	  }
      
      for(unsigned i = 0; i < simul.Global.element.size(); i++ )
	if(r.domain_number (std::ptrdiff_t(simul.Global.element[i])) == std::ptrdiff_t(r.thread_id))
	  {
	    LATT.set_coord(simul.Global.element[i] );
	    r.convertCoordinates( Latt, LATT );
	    border_element.push_back( Latt.index );	    
	    border_U.push_back( simul.Global.U[i] );
	  }
    }
#pragma omp barrier



    
    /* Test Mozaic to implement in tile that have to be set to zero  */
    for(std::size_t istr = 0; istr < r.NStr; istr++)
      for(auto it = position.at(istr).begin(); it != position.at(istr).end(); it++)
	for(unsigned node = 0; node < NumberNodes; node++ )
	  {
	    std::size_t node_pos =  *it + node_position.at(node);
	    Latt.set_coord(node_pos);
	    r.convertCoordinates(latStr, Latt );
	    if(r.test_ghosts(Latt) == 1 &&  latStr.index != istr )
	      {
		cross_mozaic[latStr.index] = false;
		if( !any_of(cross_mozaic_indexes.begin(), cross_mozaic_indexes.end(), std::bind2nd(std::equal_to<std::size_t>(), latStr.index)))
		  cross_mozaic_indexes.push_back(latStr.index);
	      }
	  }
    
    
    for(std::size_t istr = 0; istr < r.NStr; istr++)
      std::sort( position[istr].begin() , position[istr].end() );
  }
  




  
};
