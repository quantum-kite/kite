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
    std::size_t x,y;
    x = vec[0];
    y = vec[1];
    unsigned Lx,Ly;
    Lx = Ly = 128;
    double r2;
    double dx,dy;
    dx = x - Lx/2;
    dy = y - Ly/2;
    r2 = dx*dx + dy*dy;
    double V;
    V = 0;
    if(r2 < 1000){
        V = -1;
    }
    return V;
};
