/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2022, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/

#include <iostream>

// g++ -Wall -fPIC -c aux.cpp
// g++ -shared -Wl,-soname,libaux.so.1 -o libaux.so aux.o
// ln -sf libaux.so libaux.so.1

// one-liner: g++ -Wall -fPIC -c aux.cpp && g++ -shared -Wl,-soname,libaux.so.1 -o libaux.so aux.o && ln -sf libaux.so libaux.so.1

double uncorrelated_wrapper(std::size_t *vec, unsigned dim){
    return 0;
};
