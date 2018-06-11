/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

#include <iomanip>
template <typename T>
class KPM_Vector <T, 2> : public KPM_VectorBasis <T,2> {
private:
  std::size_t             max[2];
  std::size_t 	*MemIndBeg[2][2];
  std::size_t 	*MemIndEnd[2][2];
  std::size_t        block[2][2];
  std::size_t          stride[2];
  std::size_t   stride_ghosts[2];
  LatticeStructure<2u>       & r;
  Hamiltonian<T,2u>          & h;
  T               ***mult_t1_mag;
  Coordinates<std::size_t,3>   x;
  T                        *phi0;
  T                       *phiM1;
  T                       *phiM2;
  const std::size_t          std;
public:
  typedef typename extract_value_type<T>::value_type value_type;
  using KPM_VectorBasis<T,2>::simul;
  using KPM_VectorBasis<T,2>::index;
  using KPM_VectorBasis<T,2>::v;
  using KPM_VectorBasis<T,2>::memory;
  using KPM_VectorBasis<T,2>::aux_wr;
  using KPM_VectorBasis<T,2>::aux_test;
  using KPM_VectorBasis<T,2>::inc_index;
  
  KPM_Vector(int mem, Simulation<T,2> & sim) : KPM_VectorBasis<T,2>(mem, sim), r(sim.r), h(sim.h), x(r.Ld), std(x.basis[1]) {
    unsigned d;
    Coordinates <std::size_t, 3>     z(r.Ld);
    Coordinates <int, 3> x(r.nd), dist(r.nd);

    mult_t1_mag = new T**[r.Orb];
    for(unsigned io = 0; io < r.Orb; io++)
      {
	mult_t1_mag[io] = new T*[h.hr.NHoppings(io)];
	for(unsigned ib = 0; ib < h.hr.NHoppings(io); ib++)
	  mult_t1_mag[io][ib] = new T[STRIDE];
      }

    for(unsigned d = 0; d < 2; d++)
      for(unsigned b = 0; b < 2; b++)
	{
	  MemIndBeg[d][b] = new std::size_t[r.Orb];
	  MemIndEnd[d][b] = new std::size_t[r.Orb];
	    
	}
    
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
  
  
  ~KPM_Vector(void){
    for(unsigned d = 0; d < 2; d++)
      for(unsigned b = 0; b < 2; b++)
	{
	  delete MemIndBeg[d][b];
	  delete MemIndEnd[d][b];
	    
	}
    for(unsigned io = 0; io < r.Orb;io++)
      {	
	for(unsigned ib = 0; ib < h.hr.NHoppings(io); ib++)
	  delete mult_t1_mag[io][ib];
	delete mult_t1_mag[io];
      }
    delete mult_t1_mag;
  }
  
  void initiate_vector() {
    std::size_t count = 0;
    for(unsigned i = 0; i < r.NStr; i++)
      count += h.hV.position.at(i).size();
    
    index = 0;
    Coordinates<std::size_t, 3> x(r.Ld);
    for(std::size_t io = 0; io < r.Orb; io++)
      for(std::size_t i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; i1++)
	for(std::size_t i0 = NGHOSTS; i0 < r.Ld[0] - NGHOSTS; i0++)
	  v(x.set({i0,i1,io}).index, index) = simul.rnd.init()/value_type(sqrt(r.Sizet - count));
    
    for(unsigned i = 0; i < r.NStr; i++)
      {
	auto & vv = h.hV.position.at(i); 
	for(unsigned j = 0; j < vv.size(); j++)
	  v(vv.at(j), index ) = 0. ;
      }
    
  };
  
  template < unsigned MULT,bool VELOCITY> 
  void build_regular_phases(int i1, unsigned axis)
  {
    unsigned l[2 + 1], count;
    Coordinates<std::ptrdiff_t, 3>  global(r.Lt);
    Coordinates<std::ptrdiff_t, 3> local1(r.Ld);
    Coordinates<std::size_t,3> x(r.Ld);
    std::fill_n(l, 2, 3);
    l[2]  = r.Orb;
    Coordinates<std::ptrdiff_t, 2 + 1> b3(l);    
    Eigen::Map<Eigen::Matrix<std::ptrdiff_t,2, 1>> vee(b3.coord); // Column vector
    
    
    for(unsigned io = 0; io < r.Orb; io++)
      {
	const std::size_t ip = io * x.basis[2] ;
	const std::size_t j0 = ip + 0 + i1 * std;
	const std::size_t j1 = j0 + STRIDE * std;

	for(unsigned ib = 0; ib < h.hr.NHoppings(io); ib++)
	  {
	    T tt  = value_type(MULT + 1) * h.hr.hopping(ib, io);
	    b3.set_coord(h.hr.dist(ib,io));        
	    vee.array() -= 1;
	    count = 0;
	    if (VELOCITY)  tt  *=  h.hr.v.at(axis)(ib,io);
	    for(std::size_t j = j0; j < j1; j += std )
	      {
		r.convertCoordinates(global, local1.set_coord(j));
		value_type phase = vee(0)*global.coord[1]*r.vect_pot(0,1);
		mult_t1_mag[io][ib][count] =  tt * h.peierls(phase);
		count++;
	      }
	  }
      }
  }
    
  template < unsigned MULT> 
  void initiate_stride(std::size_t & istr)
  {
    std::size_t i0, i1;
    const std::size_t std = x.basis[1];
    // Periodic component of the Hamiltonian + Anderson disorder
    i0 = ((istr) % (r.lStr[0]) ) * STRIDE + NGHOSTS;
    i1 = ((istr) / r.lStr[0] ) * STRIDE + NGHOSTS;
			
    for(std::size_t io = 0; io < r.Orb; io++)
      {
	const std::size_t ip = io * x.basis[2];
	const std::size_t j0 = ip + i0 + i1 * std;
	const std::size_t j1 = j0 + STRIDE * std; //j0 and j1 define the limits of the for cycle


	for(std::size_t j = j0; j < j1; j += std )
	  for(std::size_t i = j; i < j + STRIDE ; i++)
	    phi0[i] = - value_type(MULT) * phiM2[i];
      }
  }
				
  template < unsigned MULT> 
  void inline mult_local_disorder(const  std::size_t & j0, const  std::size_t & io)
  {
    const std::size_t j1 = j0 + STRIDE * std;
    const std::ptrdiff_t dd = (h.Anderson_orb_address[io] - std::ptrdiff_t(io))*r.Nd;
    // Anderson disorder
    if( h.Anderson_orb_address[io] >= 0)
      {
	for(std::size_t j = j0; j < j1; j += std )
	  for(std::size_t i = j; i < j + STRIDE ; i++)
	    phi0[i] += value_type(MULT + 1) * phiM1[i] * h.U_Anderson.at(i + dd);
      }
    else if (h.Anderson_orb_address[io] == - 1)
      {
	for(std::size_t j = j0; j < j1; j += std )
	  for(std::size_t i = j; i < j + STRIDE ; i++)
	    phi0[i] += value_type(MULT + 1) * phiM1[i] * h.U_Orbital.at(io);
      }
  }
	      
  void inline mult_regular_hoppings(const  std::size_t & j0, const  std::size_t & io)
  {
    std::size_t count;
    const std::size_t j1 = j0 + STRIDE * std;
    // Hoppings
    for(unsigned ib = 0; ib < h.hr.NHoppings(io); ib++)
      {
	const std::ptrdiff_t d1 = h.hr.distance(ib, io);
	count = 0;
	for(std::size_t j = j0; j < j1; j += std )
	  {
	    const T t1 = mult_t1_mag[io][ib][count++];
	    for(std::size_t i = j; i < j + STRIDE ; i++)
	      phi0[i] += t1 * phiM1[i + d1];								
	  }
      }
  }
			
			
  // Structural disorder contribution - iterate over the disorder models			
  template <unsigned MULT> 
  void Multiply() {
    vverbose_message("Entered Multiply");
    
    unsigned i = 0;
    /*
      Mosaic Multiplication using a TILE of STRIDE x STRIDE 
      Right Now We expect that both ld[0] and ld[1]  are multiple of STRIDE
      MULT = 0 : For the case of the Velocity/Hamiltonian
      MULT = 1 : For the case of the KPM_iteration
    */
    inc_index();
    phi0 = v.col(index).data();
    phiM1 = v.col((memory + index - 1) % memory ).data();
    phiM2 = v.col((memory + index - 2) % memory ).data();
    KPM_MOTOR<MULT, false>(phi0, phiM1, phiM2, i);
  };


  void Velocity(T * phi0,T * phiM1, unsigned axis) {
    KPM_MOTOR<0u, true>(phi0, phiM1, phiM1, axis);
  };
  
  template <unsigned MULT, bool VELOCITY>
  void KPM_MOTOR(T * phi0a, T * phiM1a, T *phiM2a, unsigned axis)
  {
    std::size_t i0, i1;    
    phi0 = phi0a;
    phiM1 = phiM1a;
    phiM2 = phiM2a;
    
    // Initialize tiles that have deffects connecting elements of a previous tile
    for(auto istr = h.cross_mozaic_indexes.begin(); istr != h.cross_mozaic_indexes.end() ; istr++)
      initiate_stride<MULT>(*istr);
    
    for( i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; i1 += STRIDE  )
      {
	build_regular_phases<MULT,VELOCITY>(i1, axis);
		
	for( i0 = NGHOSTS; i0 < r.Ld[0] - NGHOSTS; i0 += STRIDE )
	  {
		    
	    std::size_t istr = (i1 - NGHOSTS) /STRIDE * r.lStr[0] + (i0 - NGHOSTS)/ STRIDE;
	    if(h.cross_mozaic.at(istr))
	      initiate_stride<MULT>(istr);
	    // These four lines pertrain only to the magnetic field
	    for(std::size_t io = 0; io < r.Orb; io++)
	      {
		const std::size_t ip = io * x.basis[2];
		const std::size_t j0 = ip + i0 + i1 * std;
		
		// Local Energy
		if(!VELOCITY) mult_local_disorder<MULT>(j0, io);
		
		// Hoppings
		mult_regular_hoppings(j0, io);
	      }
	    for(auto id = h.hd.begin(); id != h.hd.end(); id++)
	      id->template multiply_defect<MULT, VELOCITY>(istr, phi0, phiM1, axis);
	  	    
	    // Empty the vacancies in the tile
	    auto & hV = h.hV.position.at(istr);
	    for(auto k = hV.begin(); k != hV.end(); k++)
	      phi0[*k] = 0.;

	  }
      }

    for(auto vc =  h.hV.vacancies_with_defects.begin(); vc != h.hV.vacancies_with_defects.end(); vc++)
      phi0[*vc] = 0.;

    
    /* 
       Broken Imputirities:
       The bulk domain will receive contributions from the broken defects 
       located on the neighbour domains.
       We already subtract the vacancies from these contributions 
    */
    for(auto id = h.hd.begin(); id != h.hd.end(); id++)
      id->template multiply_broken_defect<MULT,VELOCITY>(phi0, phiM1,axis);
	  
    // These four lines pertrain only to the magnetic field
    Exchange_Boundaries();
      }
    
	  // Periodic component of the Hamiltonian + Anderson disorder
	  
  

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
		    //std::cout << "\t wrong " << std::real(v(x.index , 0)) << " " << z.index << " " << x.index << "\t\t";
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
      
