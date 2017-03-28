template <typename T>
class KPM_Vector <T, 2> : public KPM_VectorBasis <T,2> {
private:
  unsigned max[2];
  unsigned MemIndBeg[2][2];
  unsigned MemIndEnd[2][2];
  unsigned block[2][2];
public:
  typedef typename extract_value_type<T>::value_type value_type;
  using KPM_VectorBasis<T,2>::r;
  using KPM_VectorBasis<T,2>::rng;
  using KPM_VectorBasis<T,2>::index;
  using KPM_VectorBasis<T,2>::v;
  using KPM_VectorBasis<T,2>::memory;
  using KPM_VectorBasis<T,2>::h;
  using KPM_VectorBasis<T,2>::Border;
  using KPM_VectorBasis<T,2>::GlobalBorder;
  using KPM_VectorBasis<T,2>::aux_wr;
  using KPM_VectorBasis<T,2>::aux_test;
  
  KPM_Vector(int mem, Hamiltonian <T,2> & h1, LatticeStructure <2>  & r1, std::vector<T> & B , std::vector<T> & GB) : KPM_VectorBasis<T,2>(mem, h1, r1, B, GB){
    unsigned d;
    Coordinates <int, 3> z(r.Ld);
    Coordinates <int, 3> x(r.nd), dist(r.nd);


    
    d = 0;
    max[d] =  r.ld[1];
    
    MemIndBeg[d][0] = z.set({1,1,0}).index;    
    MemIndEnd[d][0] = z.set({0,1,0}).index;
    MemIndBeg[d][1] = z.set({int(r.Ld[0]) - 2, 1, 0}).index;
    MemIndEnd[d][1] = z.set({int(r.Ld[0]) - 1, 1, 0}).index; 

    d = 1;
    max[d] =  r.Ld[0];
    MemIndBeg[d][0] = z.set({0,1,0}).index;    
    MemIndEnd[d][0] = z.set({0,0,0}).index;
    MemIndBeg[d][1] = z.set({0, int(r.Ld[1]) - 2, 0}).index;
    MemIndEnd[d][1] = z.set({0, int(r.Ld[1]) - 1, 0}).index;
    printf("%d %d %d %d\n", MemIndBeg[d][0], MemIndBeg[d][1], MemIndEnd[d][0], MemIndEnd[d][1]);
    
    dist.set({0,0,0});
    for(d = 0 ; d < 2; d++)
      for(unsigned b  = 0 ; b < 2; b++)
	{
	  dist.coord[d] = int(b) * 2 - 1;
	  block[d][b] = x.set_coord( int(r.thread_id) ).add(dist).index;
	  dist.coord[d] = 0;
	}
    
    for(int d = 0; d < 2; d++ )
      for(int b = 0; b < 2; b++ )
        block[d][b] = block[d][b] * 2 * r.Orb * max[d] +   b * r.Orb * max[d]; 
    initiate_vector();
  };
  
  void initiate_vector() {

    index = 0;
    Coordinates<long, 3> x(r.Ld);
    for(unsigned io = 0; io < r.Orb; io++)
      for(unsigned i1 = 1; i1 < r.Ld[1] - 1; i1++)
	for(unsigned i0 = 1; i0 < r.Ld[0] - 1; i0++)
	  v(x.set({i0,i1,io}).index, index) = rng.init()/value_type(sqrt(r.Size));
  };

  template < unsigned MULT> 
  void Multiply(const int model) {
    /*
      Mosaic Multiplication using a TILE of STRIDE x STRIDE 
      Right Now We expect that both ld[0] and ld[1]  are multiple of STRIDE
      MULT = 0 : For the case of the Velocity/Hamiltonian
      MULT = 1 : For the case of the KPM_iteration
    */
    Coordinates<unsigned,3> x(r.Ld);
    const unsigned STRIDE0 = 4;    
    const unsigned STRIDE1 = 4;
    
    T * phi0 = v.col(index).data();
    T * phiM1 = v.col((memory + index - 1) % memory ).data();
    T * phiM2 = v.col((memory + index - 2) % memory ).data();
    
    for(unsigned io = 0; io < r.Orb; io++)
      {
	const unsigned ip = io * x.basis[2] ;
	
	for(unsigned i1 = 1; i1 < r.Ld[1] - 1; i1 += STRIDE1  )
	  for(unsigned i0 = 1; i0 < r.Ld[0] - 1; i0 += STRIDE0 )
	    {
	      const unsigned std = x.basis[1];
	      const unsigned j0 = ip + i0 + i1 * std;
	      const unsigned j1 = j0 + STRIDE1 * std;
	      
	      for(unsigned j = j0; j < j1; j += std )
		for(unsigned i = j; i < j + STRIDE0 ; i++)
		  phi0[i] = - value_type(MULT) * phiM2[i];
	      
	      for(unsigned ib = 0; ib < h.NHoppings[io]; ib++)
		{
		  const int  d1 = h.d(ib, io);
		  const T    t1 =  value_type(MULT + 1) * h.t(ib, io);
		  
		  for(unsigned j = j0; j < j1; j += std )
		    {
		      T * p0 = phi0 + j;
		      T * p1 = phiM1 + j + d1;
		      
		      for(unsigned i = j; i < j + STRIDE0 ; i++)
			p0[i] += t1 * p1[i];
		    }
		}
	    }
      }
    
    Exchange_Boundaries();    
  };

  void setCoord(unsigned i0, unsigned i1, unsigned io) {
    r.coord[0] = i0;
    r.coord[1] = i1;
    r.coord[2] = io;
  }
    
  void Exchange_Boundaries() {
    const unsigned D = 2u;
    //
    // Exchange the Borders of phi.v.col(index)
    //
    Coordinates<unsigned,3> x(r.Ld); 
    T  *phi1 = v.col(index).data();

    
    for(unsigned d = 0; d <  D; d++)  // We will transfer the orientations perpendicular  to d
      {
	unsigned int BSize = r.Orb * max[d];
	unsigned d1 = (d + 1) % D;
	const int st = x.basis[d1]; 
	
	for(unsigned io = 0; io < r.Orb; io++)
	  for(unsigned b = 0; b < 2u; b++)
	    {
	      const unsigned bi = (io + b * r.Orb) * max[d];
	      unsigned I        = MemIndBeg[d][b] + io * x.basis[D] - st;

	      for(unsigned i = 0; i < max[d]; i++ )
		Border.at(bi + i) = phi1[ I += st ];
	    }
	
	// Copy to the  shared memory
#pragma omp critical
	{
	  std::cout << Border.size() << std::endl;
	  for(unsigned i = 0; i < Border.size(); i++)
	    std::cout << Border.at(i) << " ";
	  std::cout << std::endl;
	}
	std::copy( Border.begin(), Border.begin() + 2 * BSize, GlobalBorder.begin() + 2 * BSize * r.thread_id );


#pragma omp barrier
	
	for(unsigned b = 0; b < 2; b++)
	  std::copy(GlobalBorder.begin() + block[d][b] , GlobalBorder.begin() + block[d][b] +  BSize, Border.begin() + (1-b) * BSize );
#pragma omp barrier
	
	for(unsigned io = 0; io < r.Orb; io++) // Copy back from the shared memory
	  for(unsigned b = 0; b < 2; b++)
	    {
	      const unsigned bi = (io + b * r.Orb) * max[d];
	      unsigned I        = MemIndEnd[d][b] + io * x.basis[D] - st;
	      
	      for(unsigned i = 0; i < max[d]; i++ )
		phi1[ I += st ] =  Border.at(bi + i);
	    }
      } 
  };
  
  void test_boundaries_system() {

    /*
      This  function tests if the boudaries exchange are well implemented 
    */
    
    Coordinates<long, 3> z(r.Lt);
    Coordinates<long, 3> x(r.Ld);
    for(long  io = 0; io < (long) r.Ld[2]; io++)
      for(long i1 = 1; i1 < (long) r.Ld[1] - 1 ; i1++)
	for(long i0 = 1; i0 < (long) r.Ld[0] - 1 ; i0++)
	  {
	    r.buildGlobalCoordinates( z.set({i0, i1, io}));
	    v(x.set({i0,i1,io}).index, 0) = aux_wr(z.index);
	  }
    
    Exchange_Boundaries();

#pragma omp critical
    for(long  io = 0; io < (long) r.Ld[2]; io++)
      for(long i1 = 0; i1 < (long) r.Ld[1]; i1++)
	{
	  for(long i0 = 0; i0 < (long) r.Ld[0]; i0++)
	    {
	      r.buildGlobalCoordinates( z.set({i0, i1, io}) );
	      x.set({i0,i1,io});
	      T val = aux_wr(z.index); 
	      if( aux_test(v(x.index , 0), val ) )
		std::cout << "Problems---->" << v(x.index , 0) << " " << val << std::endl;
	    };
	}
  };
};
