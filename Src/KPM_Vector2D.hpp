template <typename T>
class KPM_Vector <T, 2> : public KPM_VectorBasis <T,2> {
private:
  std::size_t             max[2];
  std::size_t MemIndBeg[2][2][2];
  std::size_t MemIndEnd[2][2][2];
  std::size_t        block[2][2];
  std::size_t          stride[2];
  std::size_t   stride_ghosts[2];
  LatticeStructure<2u>       & r;
  Hamiltonian<T,2u>          & h;
public:
  typedef typename extract_value_type<T>::value_type value_type;
  using KPM_VectorBasis<T,2>::simul;
  using KPM_VectorBasis<T,2>::index;
  using KPM_VectorBasis<T,2>::v;
  using KPM_VectorBasis<T,2>::memory;
  using KPM_VectorBasis<T,2>::aux_wr;
  using KPM_VectorBasis<T,2>::aux_test;
  using KPM_VectorBasis<T,2>::inc_index;
  
  KPM_Vector(int mem, Simulation<T,2> & sim) : KPM_VectorBasis<T,2>(mem, sim), r(sim.r), h(sim.h){
    unsigned d;
    Coordinates <std::size_t, 3>     z(r.Ld);
    Coordinates <int, 3> x(r.nd), dist(r.nd);

    for(std::size_t io = 0; io < r.Orb; io++)
      {
	d = 0;
	max[d] =  r.ld[1];
	stride[d]  = r.Ld[0];
	stride_ghosts[d] = 1;
	MemIndBeg[d][0][io] = z.set({std::size_t(NGHOSTS),               std::size_t(NGHOSTS), io}).index;   // From the index Starting here --- Right
	MemIndEnd[d][0][io] = z.set({std::size_t(0),                     std::size_t(NGHOSTS), io}).index;   //   To the index Starting here --- Right
	MemIndBeg[d][1][io] = z.set({std::size_t(r.Ld[0] - 2 * NGHOSTS), std::size_t(NGHOSTS), io}).index;   // From the index Starting here --- Left
	MemIndEnd[d][1][io] = z.set({std::size_t(r.Ld[0] - NGHOSTS),     std::size_t(NGHOSTS), io}).index;   //   To the index Starting here --- Left
	
	d = 1;
	max[d] =  r.Ld[0];
	stride[d] = 1;
	stride_ghosts[d] = r.Ld[0];
	MemIndBeg[d][0][io] = z.set({std::size_t(0), std::size_t(NGHOSTS),             io}).index;          // From the index Starting here --- Bottom
	MemIndEnd[d][0][io] = z.set({std::size_t(0), std::size_t(0),                   io}).index;          //   To the index Starting here --- Bottom
	MemIndBeg[d][1][io] = z.set({std::size_t(0), std::size_t(r.Ld[1] - 2*NGHOSTS), io}).index;          // From the index Starting here --- Top
	MemIndEnd[d][1][io] = z.set({std::size_t(0), std::size_t(r.Ld[1] - NGHOSTS),   io}).index;          //   To the index Starting here --- Top
      }
    
    for(d = 0 ; d < 2; d++)
      for(unsigned b  = 0 ; b < 2; b++)
	{
	  dist.set({0,0,0});
	  dist.coord[d] = int(b) * 2 - 1;
	  block[d][b] = x.set_coord( int(r.thread_id) ).add(dist).index;
	}
    initiate_vector();
  };
  
  void initiate_vector() {
    index = 0;
    Coordinates<std::size_t, 3> x(r.Ld);
    for(std::size_t io = 0; io < r.Orb; io++)
      for(std::size_t i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; i1++)
	for(std::size_t i0 = NGHOSTS; i0 < r.Ld[0] - NGHOSTS; i0++)
	  v(x.set({i0,i1,io}).index, index) = simul.rnd.init()/value_type(sqrt(r.Sizet));
  };
  
  template < unsigned MULT> 
  void Multiply() {
    /*
      Mosaic Multiplication using a TILE of STRIDE x STRIDE 
      Right Now We expect that both ld[0] and ld[1]  are multiple of STRIDE
      MULT = 0 : For the case of the Velocity/Hamiltonian
      MULT = 1 : For the case of the KPM_iteration
    */

    Coordinates<std::size_t,3> x(r.Ld);
    std::size_t i0, i1;
    inc_index();
    T * phi0 = v.col(index).data();
    T * phiM1 = v.col((memory + index - 1) % memory ).data();
    T * phiM2 = v.col((memory + index - 2) % memory ).data();
    
    // Initialize tiles that have deffects connecting elements of a previous tile
    for(auto istr = h.cross_mozaic_indexes.begin(); istr != h.cross_mozaic_indexes.end() ; istr++)
      {
	i0 = ((*istr) % (r.lStr[0]) ) * STRIDE + NGHOSTS;
	i1 = ((*istr) / r.lStr[0] ) * STRIDE + NGHOSTS;
	for(std::size_t io = 0; io < r.Orb; io++)
	  {
	    const std::size_t ip = io * x.basis[2] ;
	    const std::size_t std = x.basis[1];
	    const std::size_t j0 = ip + i0 + i1 * std;
	    const std::size_t j1 = j0 + STRIDE * std;
	    
	    for(std::size_t j = j0; j < j1; j += std )
	      for(std::size_t i = j; i < j + STRIDE ; i++)
		phi0[i] = - value_type(MULT) * phiM2[i];
	  }
      }
    
    for( i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; i1 += STRIDE  )
      for( i0 = NGHOSTS; i0 < r.Ld[0] - NGHOSTS; i0 += STRIDE )
	{
	  // Periodic component of the Hamiltonian + Anderson disorder
	  std::size_t istr = (i1 - NGHOSTS) /STRIDE * r.lStr[0] + (i0 - NGHOSTS)/ STRIDE;
	  
	  for(std::size_t io = 0; io < r.Orb; io++)
	    {
	      const std::size_t ip = io * x.basis[2] ;
	      const std::size_t dd = (h.Anderson_orb_address[io] - io)*r.Nd;
	      const std::size_t std = x.basis[1];
	      const std::size_t j0 = ip + i0 + i1 * std;
	      const std::size_t j1 = j0 + STRIDE * std;

	      // Initialize phi0
	      for(unsigned id = 0; id < h.hd.size(); id++)
		if(h.cross_mozaic.at(istr))
		  for(std::size_t j = j0; j < j1; j += std )
		    for(std::size_t i = j; i < j + STRIDE ; i++)
		      phi0[i] = - value_type(MULT) * phiM2[i];
	      

	      // Anderson disorder
	      if( h.Anderson_orb_address[io] >= 0)
		{
		  for(std::size_t j = j0; j < j1; j += std )
		    for(std::size_t i = j; i < j + STRIDE ; i++)
		      phi0[i] += value_type(MULT + 1) * phiM1[i] * h.U_Anderson.at(i - dd);
		}
	      else if (h.Anderson_orb_address[io] == - 1)
		{
		  for(std::size_t j = j0; j < j1; j += std )
		    for(std::size_t i = j; i < j + STRIDE ; i++)
		      phi0[i] += value_type(MULT + 1) * phiM1[i] * h.U_Orbital.at(io);
		}
	      // Hoppings
	      for(unsigned ib = 0; ib < h.hr.NHoppings(io); ib++)
		{
		  const std::ptrdiff_t  d1 = h.hr.distance(ib, io);
		  const T t1 =  value_type(MULT + 1) * h.hr.hopping(ib, io);
		  
		  for(std::size_t j = j0; j < j1; j += std )
		    for(std::size_t i = j; i < j + STRIDE ; i++)
		      phi0[i] += t1 * phiM1[i + d1];
		}
	    }
	  
	  // Structural disorder contribution
	  
	  for(auto id = h.hd.begin(); id != h.hd.end(); id++)
	    for(std::size_t i = 0; i <  id->position.at(istr).size(); i++)
	      {
		std::size_t ip = id->position.at(istr)[i];
		for(unsigned k = 0; k < id->hopping.size(); k++)
		  {
		    std::size_t k1 = ip + id->node_position[id->element1[k]];
		    std::size_t k2 = ip + id->node_position[id->element2[k]];
		    phi0[k1] += value_type(MULT + 1) * id->hopping[k] * phiM1[k2];
		  }
		
		for(std::size_t k = 0; k < id->U.size(); k++)
		  {
		    std::size_t k1 = ip + id->node_position[id->element[k]];
		    phi0[k1] += value_type(MULT + 1) * id->U[k] * phiM1[k1];
		  }
	      }
	}
    
    //  Broken impurities
    for(auto id = h.hd.begin(); id != h.hd.end(); id++)
      for(std::size_t i = 0; i < id->border_element1.size(); i++)
	{
	  std::size_t i1 = id->border_element1[i];
	  std::size_t i2 = id->border_element2[i];
	  phi0[i1] += value_type(MULT + 1) * id->border_hopping[i] * phiM1[i2];
	}
    
    for(auto id = h.hd.begin(); id != h.hd.end(); id++)
      for(std::size_t i = 0; i < id->border_element.size(); i++)
	{
	  std::size_t i1 = id->border_element[i];
	  phi0[i1] += value_type(MULT + 1) * id->border_U[i] * phiM1[i1];
	}
    
    Exchange_Boundaries();
  };
  
  void Velocity( T * phi0, T * phiM1, int axis) {
    /*
      Mosaic Multiplication using a TILE of STRIDE x STRIDE 
      Right Now We expect that both ld[0] and ld[1]  are multiple of STRIDE
      MULT = 0 : For the case of the Velocity/Hamiltonian
      MULT = 1 : For the case of the KPM_iteration
    */
    T zero = assign_value<T>(double(0), double(0));
    Coordinates<std::size_t,3> x(r.Ld);
    std::size_t i0, i1;
    
    // Initialize tiles that have deffects connecting elements of a previous tile
    for(auto istr = h.cross_mozaic_indexes.begin(); istr != h.cross_mozaic_indexes.end() ; istr++)
      {
	i0 = ((*istr) % (r.lStr[0]) ) * STRIDE + NGHOSTS;
	i1 = ((*istr) / r.lStr[0] ) * STRIDE + NGHOSTS;
	for(std::size_t io = 0; io < r.Orb; io++)
	  {
	    const std::size_t ip = io * x.basis[2] ;
	    const std::size_t std = x.basis[1];
	    const std::size_t j0 = ip + i0 + i1 * std;
	    const std::size_t j1 = j0 + STRIDE * std;
	    
	    for(std::size_t j = j0; j < j1; j += std )
	      for(std::size_t i = j; i < j + STRIDE ; i++)
		phi0[i] = zero;
	  }
      }
    
    for( i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; i1 += STRIDE  )
      for( i0 = NGHOSTS; i0 < r.Ld[0] - NGHOSTS; i0 += STRIDE )
	{
	  // Periodic component of the Hamiltonian + Anderson disorder
	  std::size_t istr = (i1 - NGHOSTS) /STRIDE * r.lStr[0] + (i0 - NGHOSTS)/ STRIDE;
	  
	  for(std::size_t io = 0; io < r.Orb; io++)
	    {
	      const std::size_t ip = io * x.basis[2] ;
	      const std::size_t std = x.basis[1];
	      const std::size_t j0 = ip + i0 + i1 * std;
	      const std::size_t j1 = j0 + STRIDE * std;

	      // Initialize phi0
	      for(unsigned id = 0; id < h.hd.size(); id++)
		if(h.cross_mozaic.at(istr))
		  for(std::size_t j = j0; j < j1; j += std )
		    for(std::size_t i = j; i < j + STRIDE ; i++)
		      phi0[i] = zero;
	      
	      // Hoppings
	      for(unsigned ib = 0; ib < h.hr.NHoppings(io); ib++)
		{
		  const std::ptrdiff_t  d1 = h.hr.distance(ib, io);
		  const T t1 =   h.hr.V[axis](ib, io);
		  
		  for(std::size_t j = j0; j < j1; j += std )
		    for(std::size_t i = j; i < j + STRIDE ; i++)
		      phi0[i] += t1 * phiM1[i + d1];
		}
	    }
	  
	  // Structural disorder contribution
	  
	  for(auto id = h.hd.begin(); id != h.hd.end(); id++)
	    for(std::size_t i = 0; i <  id->position.at(istr).size(); i++)
	      {
		std::size_t ip = id->position.at(istr)[i];
		for(unsigned k = 0; k < id->hopping.size(); k++)
		  {
		    std::size_t k1 = ip + id->node_position[id->element1[k]];
		    std::size_t k2 = ip + id->node_position[id->element2[k]];
		    phi0[k1] +=  (id->V[axis])[k] * phiM1[k2];
		  }
	      }
	}
    
    //  Broken impurities
    for(auto id = h.hd.begin(); id != h.hd.end(); id++)
      for(std::size_t i = 0; i < id->border_element1.size(); i++)
	{
	  std::size_t i1 = id->border_element1[i];
	  std::size_t i2 = id->border_element2[i];
	  phi0[i1] +=  (id->border_V[axis])[i] * phiM1[i2];
	}
        
    Exchange_Boundaries();
  };



  void Velocity2( T * phi0, T * phiM1, int axis1, int axis2) {
    /*
      Mosaic Multiplication using a TILE of STRIDE x STRIDE 
      Right Now We expect that both ld[0] and ld[1]  are multiple of STRIDE
      MULT = 0 : For the case of the Velocity/Hamiltonian
      MULT = 1 : For the case of the KPM_iteration
    */
    T zero = assign_value<T>(double(0), double(0));
    Coordinates<std::size_t,3> x(r.Ld);
    std::size_t i0, i1;
    
    // Initialize tiles that have deffects connecting elements of a previous tile
    for(auto istr = h.cross_mozaic_indexes.begin(); istr != h.cross_mozaic_indexes.end() ; istr++)
      {
	i0 = ((*istr) % (r.lStr[0]) ) * STRIDE + NGHOSTS;
	i1 = ((*istr) / r.lStr[0] ) * STRIDE + NGHOSTS;
	for(std::size_t io = 0; io < r.Orb; io++)
	  {
	    const std::size_t ip = io * x.basis[2] ;
	    const std::size_t std = x.basis[1];
	    const std::size_t j0 = ip + i0 + i1 * std;
	    const std::size_t j1 = j0 + STRIDE * std;
	    
	    for(std::size_t j = j0; j < j1; j += std )
	      for(std::size_t i = j; i < j + STRIDE ; i++)
		phi0[i] = zero;
	  }
      }
    
    for( i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; i1 += STRIDE  )
      for( i0 = NGHOSTS; i0 < r.Ld[0] - NGHOSTS; i0 += STRIDE )
	{
	  // Periodic component of the Hamiltonian + Anderson disorder
	  std::size_t istr = (i1 - NGHOSTS) /STRIDE * r.lStr[0] + (i0 - NGHOSTS)/ STRIDE;
	  
	  for(std::size_t io = 0; io < r.Orb; io++)
	    {
	      const std::size_t ip = io * x.basis[2] ;
	      const std::size_t std = x.basis[1];
	      const std::size_t j0 = ip + i0 + i1 * std;
	      const std::size_t j1 = j0 + STRIDE * std;

	      // Initialize phi0
	      for(unsigned id = 0; id < h.hd.size(); id++)
		if(h.cross_mozaic.at(istr))
		  for(std::size_t j = j0; j < j1; j += std )
		    for(std::size_t i = j; i < j + STRIDE ; i++)
		      phi0[i] = zero;
	      
	      // Hoppings
	      for(unsigned ib = 0; ib < h.hr.NHoppings(io); ib++)
		{
		  const std::ptrdiff_t  d1 = h.hr.distance(ib, io);
		  const T t1 =   h.hr.V2[axis1][axis2](ib, io);
		  
		  for(std::size_t j = j0; j < j1; j += std )
		    for(std::size_t i = j; i < j + STRIDE ; i++)
		      phi0[i] += t1 * phiM1[i + d1];
		}
	    }
	  
	  // Structural disorder contribution
	  
	  for(auto id = h.hd.begin(); id != h.hd.end(); id++)
	    for(std::size_t i = 0; i <  id->position.at(istr).size(); i++)
	      {
		std::size_t ip = id->position.at(istr)[i];
		for(unsigned k = 0; k < id->hopping.size(); k++)
		  {
		    std::size_t k1 = ip + id->node_position[id->element1[k]];
		    std::size_t k2 = ip + id->node_position[id->element2[k]];
		    phi0[k1] +=  (id->V2[axis1][axis2])[k]* phiM1[k2];
		  }
	      }
	}
    
    //  Broken impurities
    for(auto id = h.hd.begin(); id != h.hd.end(); id++)
      for(std::size_t i = 0; i < id->border_element1.size(); i++)
	{
	  std::size_t i1 = id->border_element1[i];
	  std::size_t i2 = id->border_element2[i];
	  phi0[i1] +=  (id->border_V2[axis1][axis2])[i] * phiM1[i2];
	}
        
    Exchange_Boundaries();
  };

  
  T VelocityInternalProduct( T * bra, T * ket, int axis)
  {
    Coordinates<std::size_t,3> x(r.Ld);
    const std::size_t STRIDE0 = 4;    
    const std::size_t STRIDE1 = 4;
    typedef typename extract_value_type<T>::value_type value_type;
    T sum;
    sum *= value_type(0.);
    
    for(std::size_t io = 0; io < r.Orb; io++)
      {
	const std::size_t ip = io * x.basis[2];
	for(std::size_t i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; i1 += STRIDE1  )
	  for(std::size_t i0 = NGHOSTS; i0 < r.Ld[0] - NGHOSTS; i0 += STRIDE0 )
	    {
	      const std::size_t std = x.basis[1];
	      const std::size_t j0 = ip + i0 + i1 * std;
	      const std::size_t j1 = j0 + STRIDE1 * std;
	      
	      for(std::size_t ib = 0; ib < h.hr.NHoppings(io); ib++)
		{
		  const  std::ptrdiff_t  d1 = h.hr.distance(ib, io);
		  const T    t1 = h.hr.V[axis](ib, io);
		  
		  for(std::size_t j = j0; j < j1; j += std )
		    for(std::size_t i = j; i < j + STRIDE0 ; i++)
		      sum += std::conj(bra[i]) * t1 * ket[i + d1];
		}
	    }
      }
    return sum;
  };

  

  template < unsigned MULT> 
  void Multiply2() {
    
    inc_index();
    T * phi0 = v.col(index).data();
    T * phiM1 = v.col((memory + index - 1) % memory ).data();
    T * phiM2 = v.col((memory + index - 2) % memory ).data();
    
    for(std::size_t io = 0; io < int(r.Orb); io++)
      for(std::size_t iy = NGHOSTS; iy < int(r.Ld[1]) - NGHOSTS; iy++)
	for(std::size_t ix = NGHOSTS; ix < int(r.Ld[0]) - NGHOSTS; ix++)
	  {
	    std::size_t i = ix + iy * r.Ld[0] + io * r.Nd;  
	    phi0[i] = - value_type(MULT) * phiM2[i];
	    for(std::size_t ib = 0; ib < h.NHoppings(io); ib++)
	      phi0[i] +=  value_type(MULT + 1) * h.t(ib, io) * phiM1[i + h.d(ib, io) ];
	  }
    
    
    Exchange_Boundaries();    
  };
  

  void Exchange_Boundaries() {
    /*
      I have four boundaries to exchange with the other threads.
      First I will copy the lines along the a[1] direction to a consecutive shared vector
    */
#pragma omp barrier    
    Coordinates<std::size_t,3u> x(r.Ld), z(r.Lt);
    T  *phi = v.col(index).data();
    
    
    for(unsigned d = 0; d < 2; d++)
      {
	
	std::size_t BSize = r.Orb * max[d] * NGHOSTS;
	T * ghosts_left = & simul.ghosts[0];
	T * ghosts_right = & simul.ghosts[BSize];
    
	for(std::size_t io = 0; io < r.Orb; io++)
	  {
	    std::size_t il = MemIndBeg[d][0][io];
	    std::size_t ir = MemIndBeg[d][1][io];
	    
	    for(std::size_t i = 0; i < max[d]; i++)
	      {
		for(unsigned ig = 0; ig < NGHOSTS; ig++)
		  {
		    ghosts_left [i + (ig + NGHOSTS*io) * max[d] ] = phi[il + ig*stride_ghosts[d]];
		    ghosts_right[i + (ig + NGHOSTS*io) * max[d] ] = phi[ir + ig*stride_ghosts[d]];
		  }
		
		il += stride[d];
		ir += stride[d];
	      }
	  }
	
	// Copy the boundaries to the shared memory
	std::copy( ghosts_left, ghosts_left + 2 * BSize, simul.Global.ghosts.begin() + 2 * BSize * r.thread_id );	  
#pragma omp barrier
	auto neigh_left = simul.Global.ghosts.begin() + 2 * block[d][0] * BSize;
	auto neigh_right  = simul.Global.ghosts.begin() + 2 * block[d][1] * BSize;
	std::copy(neigh_right,         neigh_right + BSize , ghosts_right );     // From the left to the right
	std::copy(neigh_left + BSize,  neigh_left + 2*BSize, ghosts_left  )  ;   // From the right to the left

#pragma omp barrier	
	for(std::size_t io = 0; io < r.Orb; io++)
	  {
	    std::size_t il = MemIndEnd[d][0][io];
	    std::size_t ir = MemIndEnd[d][1][io];
	    
	    for(std::size_t i = 0; i < max[d]; i++)
	      {
		for(int ig = 0; ig < NGHOSTS; ig++)
		  {
		    phi[il + ig*stride_ghosts[d]] = ghosts_left [i + (ig + NGHOSTS*io) * max[d]];
		    phi[ir + ig*stride_ghosts[d]] = ghosts_right[i + (ig + NGHOSTS*io) * max[d]];
		  }
		
		il += stride[d];
		ir += stride[d];
	      }
	  }
      }
  }
  

  /*  
  void Exchange_Boundaries() {

    const unsigned D = 2u;
    //
    // Exchange the Borders of phi.v.col(index)
    //
    Coordinates<unsigned,3> x(r.Ld); 
    T  *phi1 = v.col(index).data();
    
    
    for(unsigned d = 0; d <  D; d++)  // We will transfer the orientations perpendicular  to d
      {
	unsigned int BSize = r.Orb * max[d] * NGHOSTS;
	unsigned d1 = (d + 1) % D;
	const int st = x.basis[d1]; 
	
	for(unsigned io = 0; io < r.Orb; io++)
	  for(unsigned b = 0; b < 2u; b++)
	    {
	      const unsigned bi = (io + b * r.Orb) * max[d] * NGHOSTS;
	      unsigned I        = MemIndBeg[d][b] + io * x.basis[D] - st;
	      for(unsigned i = 0; i < max[d]; i++ )
		simul.ghosts.at(bi + i) = phi1[ I += st ];
	    }


	
	// Copy to the  shared memory
	std::copy( simul.ghosts.begin(), simul.ghosts.begin() + 2 * BSize, simul.Global.ghosts.begin() + 2 * BSize * r.thread_id );

	
#pragma omp barrier
	
	for(unsigned b = 0; b < 2; b++)
	  std::copy(simul.Global.ghosts.begin() + block[d][b] , simul.Global.ghosts.begin() + block[d][b] +  BSize, simul.ghosts.begin() + (1-b) * BSize );
#pragma omp barrier
	
	for(unsigned io = 0; io < r.Orb; io++) // Copy back from the shared memory
	  for(unsigned b = 0; b < 2; b++)
	    {
	      const unsigned bi = (io + b * r.Orb) * max[d];
	      unsigned I        = MemIndEnd[d][b] + io * x.basis[D] - st;
	      
	      for(unsigned i = 0; i < max[d]; i++ )
		phi1[ I += st ] =  simul.ghosts.at(bi + i);
	    }
      } 

  };
    */  
  void test_boundaries_system() {

    /*
      This  function tests if the boudaries exchange are well implemented 
    */
    
    Coordinates<std::size_t, 3> z(r.Lt);
    Coordinates<std::size_t, 3> x(r.Ld);
    
    for(std::size_t  io = 0; io < (std::size_t) r.Ld[2]; io++)
      for(std::size_t i1 = NGHOSTS; i1 < (std::size_t) r.Ld[1] - NGHOSTS ; i1++)
	for(std::size_t i0 = NGHOSTS; i0 < (std::size_t) r.Ld[0] - NGHOSTS ; i0++)
	  {
	    r.convertCoordinates(z, x.set({i0,i1,io}) );
	    v(x.set({i0,i1,io}).index, 0) = aux_wr(z.index);
	  }
    
    Exchange_Boundaries();

#pragma omp critical
    {

      for(std::size_t  io = 0; io < (std::size_t) r.Ld[2]; io++)
	for(std::size_t i1 = 0; i1 < (std::size_t) r.Ld[1] ; i1++)
	  {
	    for(std::size_t i0 = 0; i0 < (std::size_t) r.Ld[0]; i0++)
	      {
		r.convertCoordinates(z, x.set({i0,i1,io}) );
		T val = aux_wr(z.index); 
		if( aux_test(v(x.index , 0), val ) )
		  {
		    // std::cout << "Problems---->" << v(x.index , 0) << " " << val << std::endl;
		    std::cout << "\t wrong " << std::real(v(x.index , 0)) << " " << z.index << " " << x.index << "\t\t";
		    x.print();
		    
		  }
	      };
	  }
    }

  };

  void empty_ghosts(int mem_index) {
    /* This function takes the kpm vector that's being used, 'v' and sets to zero the part corresponding
     * to the ghosts, that is, the part of the vector that actually belongs to a different thread.
     * This is done so that when we take the dot product 'v' with another vector only terms pertraining 
     * to the current thread are considered.
     * */
  
    Coordinates<long, 3> x(r.Ld);
    
    
    // There are four sides, so set the ghosts in each side to zero individually.
    // Remember that the size of the ghost boundaries depends on NGHOSTS.
    
    for(long  io = 0; io < (long) r.Ld[2]; io++)
      for(long i0 = 0; i0 < (long) r.Ld[0]; i0++)
	for(int d = 0; d < NGHOSTS; d++)
	  v(x.set({i0,(long) d,io}).index, mem_index) *= 0;

    for(long  io = 0; io < (long) r.Ld[2]; io++)
      for(long i0 = 0; i0 < (long) r.Ld[0]; i0++)
	for(int d = 0; d < NGHOSTS; d++)
	  v(x.set({i0, (long) (r.Ld[1] - 1 - d),io}).index, mem_index) *= 0;
  
      for(long  io = 0; io < (long) r.Ld[2]; io++)
	for(long i1 = 0; i1 < (long) r.Ld[1]; i1++)
	  for(int d = 0; d < NGHOSTS; d++)
	    v(x.set({(long) d,i1,io}).index, mem_index) *= 0;

      for(long  io = 0; io < (long) r.Ld[2]; io++)
	for(long i1 = 0; i1 < (long) r.Ld[1]; i1++)
	  for(int d = 0; d < NGHOSTS; d++)
	    v(x.set({(long) (r.Ld[0] - 1 - d),i1,io}).index, mem_index) *= 0;

    };
  
};
