#include <iostream>

// g++ -Wall -fPIC -c aux.cpp
// g++ -shared -Wl,-soname,libaux.so.1 -o libaux.so aux.o
// ln -sf libaux.so libaux.so.1

// one-liner: g++ -Wall -fPIC -c aux.cpp && g++ -shared -Wl,-soname,libaux.so.1 -o libaux.so aux.o && ln -sf libaux.so libaux.so.1

double uncorrelated_wrapper(std::size_t *vec, unsigned dim){
    return 0;
};
