#include "Coordinates.hpp"

template <typename T, unsigned D>
Coordinates<T,D>::Coordinates(T x, unsigned (&b)[D] ) : L(b) {
  buildBasis();
  set_coord(x);
}
  
template <typename T, unsigned D>
Coordinates<T,D>::Coordinates(T (&coord)[D], unsigned (&b)[D] ) : L(b) {
  buildBasis();
  set_index(coord);
}
  
template <typename T, unsigned D>
Coordinates<T,D>::Coordinates(unsigned (&b)[D] ) : L(b) {
  buildBasis();
  index = 0;
}

template <typename T, unsigned D>
void  Coordinates<T,D>::buildBasis() {
  T p = 1;
  
  for(unsigned i = 0; i < D; i++)
    {
      basis[i] = p;
      p *= T(L[i]);
    }
}

template <typename T, unsigned D>
void Coordinates<T,D>::print() {
    for(unsigned i = 0; i < D ; i++)
      std::cout << coord[i] << " ";
    std::cout << std::endl;
  }

template <typename T, unsigned D>
Coordinates<T,D> & Coordinates<T,D>::set(std::initializer_list<T> a_args)  {
  int k = 0;
  for (auto i: a_args) coord[k++] = T(i);
  set_index(coord);
  return *this;
}


template <typename T, unsigned D>  
Coordinates<T,D> & Coordinates<T,D>::set_coord(T x ) {
  index = x;
  for(int i = D - 1; i >= 0; i--)
    {
      coord[i] = x / basis[i];
      x = x % basis[i];
    };
  return *this;
}

template <typename T, unsigned D>  
Coordinates<T,D> & Coordinates<T,D>::add( Coordinates<T,D> & x) {
  for(int i = 0; i < int(D); i++)
    coord[i] = (L[i] + coord[i] + x.coord[i]) % L[i];
  
  set_index(coord);
  return *this;
}

template <typename T, unsigned D>
Coordinates<T,D> & Coordinates<T,D>::subtract( Coordinates<T,D> & x) {
  for(int i = 0; i < int(D); i++)
    coord[i] = coord[i] - x.coord[i];
  set_index(coord);
  return *this;
}

template struct Coordinates<std::size_t,1u>;
template struct Coordinates<std::size_t,2u>;
template struct Coordinates<std::size_t,3u>;
template struct Coordinates<std::size_t,4u>;

template struct Coordinates<long,1u>;
template struct Coordinates<long,2u>;
template struct Coordinates<long,3u>;
template struct Coordinates<long,4u>;

template struct Coordinates<int,1u>;
template struct Coordinates<int,2u>;
template struct Coordinates<int,3u>;
template struct Coordinates<int,4u>;

template struct Coordinates<unsigned,1u>;
template struct Coordinates<unsigned,2u>;
template struct Coordinates<unsigned,3u>;
template struct Coordinates<unsigned,4u>;
