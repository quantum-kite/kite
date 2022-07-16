#include <iostream>

// g++ -Wall -fPIC -c aux.cpp
// g++ -shared -Wl,-soname,libaux.so.1 -o libaux.so aux.o
// ln -sf libaux.so libaux.so.1

// one-liner: g++ -Wall -fPIC -c aux.cpp && g++ -shared -Wl,-soname,libaux.so.1 -o libaux.so aux.o && ln -sf libaux.so libaux.so.1

double uncorrelated_wrapper(double *vec, unsigned orb){
    // vec is a pointer to a vector of dimension D
    // where D is the dimensionality of the lattice
    //
    // Example: in 2D, vec[0]=x vec[1]=y 
    double x,y;
    x = vec[0];
    y = vec[1];

    // Example: in 3D, vec[0]=x vec[1]=y vec[2]=z 
    //double x,y,z;
    //x = vec[0];
    //y = vec[1];
    //z = vec[2];

    double V = 0;

    // Example for the potential: circular quantum well
    unsigned Lx,Ly;
    Lx = Ly = 512;
    double r2;
    double dx,dy;
    dx = x - Lx/2;
    dy = y - Ly/2;
    r2 = dx*dx + dy*dy;

    if(r2 < 10000 and orb == 0)
        V = -1.1;

    if(100 < x and x < 200 and 300 < y and y < 400)
        V = 0.6;

    if(orb==1 and x > 500)
        V = 1.0;

    return V;
};



