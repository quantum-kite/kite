//
#include <iomanip>
#define v(x,L) ((L+x)%L)

template <typename T>
class Hamiltonian2D : public Hamiltonian <T> {
private:
  /*
    Regular Lattice Hamiltonian
  */
  int ith;
  int borderI[2][2];
  int borderF[2][2];
  int block[2][2]; 
  int border_index[2][2];
  int max[2][2];
  const int OStride, Stride[2];
  
  
public:
  using Hamiltonian<T>::Size;
  using Hamiltonian<T>::N;
  using Hamiltonian<T>::Orb;
  using Hamiltonian<T>::Nd;
  using Hamiltonian<T>::dim;
  using Hamiltonian<T>::Ld;
  using Hamiltonian<T>::Lt;
  using Hamiltonian<T>::ld;
  using Hamiltonian<T>::dist;
  using Hamiltonian<T>::GlobalBorders;
  using Hamiltonian<T>::MyBorders;
  using Hamiltonian<T>::nd;
  using Hamiltonian<T>::NHoppings;
  using Hamiltonian<T>::d;
  using Hamiltonian<T>::t;
  
  inline int get_lat_numb(unsigned i0, unsigned i1, unsigned io)
  {
    return (i0 + Ld[0] * i1) + io * Nd;
  };
  
  Hamiltonian2D(char *argv, class Simulation<T> & s1 ):
    Hamiltonian<T>(argv, s1 ), OStride(s1.Nd), Stride{1,s1.Ld[0]} {
    ith = omp_get_thread_num();
    int i0 = ith % nd[0], i1 = ith / nd[0];
    
    borderI[0][0] = get_lat_numb(1, 1, 0);
    borderI[0][1] = get_lat_numb(Ld[0] - 2, 1, 0);
    borderI[1][0] = get_lat_numb(0, 1, 0);
    borderI[1][1] = get_lat_numb(0, Ld[1] - 2, 0);

    borderF[0][0] = get_lat_numb(Ld[0] - 1, 1,   0);
    borderF[0][1] = get_lat_numb(0, 1, 0);
    borderF[1][0] = get_lat_numb(0, Ld[1] - 1, 0);
    borderF[1][1] = get_lat_numb(0, 0, 0);

    max[0][0] = ld[1] ;
    max[0][1] = ld[1] ;
    max[1][0] = Ld[0] ;
    max[1][1] = Ld[0] ;

    border_index[0][0] = 0;
    border_index[0][1] = Orb * ld[1];
    border_index[1][0] = 0;
    border_index[1][1] = Orb * Ld[0];

    block[0][0] = ith + v(i0 - 1, nd[0]) - i0;
    block[0][1] = ith + v(i0 + 1, nd[0]) - i0; 
    block[1][0] = ith + nd[0] * (v(i1 - 1, nd[1]) - i1);
    block[1][1] = ith + nd[0] * (v(i1 + 1, nd[1]) - i1);

    for(int d = 0; d < 2; d++ )
      for(int b = 0; b < 2; b++ )
	block[d][b] = block[d][b] * 2 * Orb * max[d][b] +   b*Orb * max[d][b]; 
  };
  
  
  

  long int get_full_lattice_number(int & i) {
    int i0 = ith % nd[0], i1 = ith / nd[0];
    int ix = i % Ld[0], iy = i/Ld[0];
    
    long int  IX = (ix - 1 + i0 * ld[0] + Lt[0]) % Lt[0];
    long int  IY = (iy - 1 + i1 * ld[1] + Lt[1]) % Lt[1];
    return IX + IY * ((long int) Lt[0]);
  };

  template <typename T1>
  void convert(T1 & x, int i) {
    x = T1(get_full_lattice_number(i) );
  }

  template <typename T1>
  void convert(std::complex<T1> & x, int i) {
    x = std::complex<T1>(get_full_lattice_number(i), 0 );
  }

  void test_boundary_exchange(KPM_Vector<T> &  phi) {
    T  *phi1 = phi.get_vector(0);
    
    for(int io = 0; io < Orb; io++)
      for(int i1 = 0; i1 < Ld[1]; i1++)
	for(int i0 = 0; i0 < Ld[0]; i0++)
	  {
	    int i =  i0 + i1 * Ld[0];
	    if(i1 == 0 || i0 == 0 || i1 == Ld[1] - 1 || i0 == Ld[0] -1)
	      phi1[i + io*OStride] = 0.;
	    else
	      convert(phi1[i + io*OStride] , i);
	  }
    
    Exchange_Boundaries(phi);
#pragma omp critical

    for(int io = 0; io < Orb; io++)
      for(int i1 = 0; i1 < Ld[1]; i1++)
	for(int i0 = 0; i0 < Ld[0]; i0++)
	  {
	    int i =   i0 + i1 * Ld[0];
	    if(fabs(sqrt(std::norm(phi1[i + io*OStride])) - get_full_lattice_number(i)) > 0.5 )
	      printf(" WRONG\n");
	  }
  };
  
};




