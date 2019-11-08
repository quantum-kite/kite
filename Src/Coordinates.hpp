#include "Generic.hpp"


template <typename T, unsigned D>
struct Coordinates;

template <typename T, unsigned D>
struct Coordinates {
  T index = 0;
  T coord[D] = {};
  unsigned (&L)[D];
  T basis[D];
  Coordinates(T x, unsigned (&b)[D]);
  Coordinates(T (&coord)[D], unsigned (&b)[D]);
  Coordinates(unsigned (&b)[D]);
  void  buildBasis();  
  void print();
  Coordinates<T,D> & set(std::initializer_list<T> a_args);

  template <typename T1>
  Coordinates<T,D> & set_index(T1 (&c)[D] ) {
    index = 0;
    for(int i = D - 1; i >= 0; i--)
      {
        coord[i] = T(c[i]);
        index += T(c[i]) * basis[i];
      }
    return *this;
  }

  Coordinates<T,D> & set_coord(T x );
  Coordinates<T,D> & add( Coordinates<T,D> & x);
  Coordinates<T,D> & subtract( Coordinates<T,D> & x);  
};





