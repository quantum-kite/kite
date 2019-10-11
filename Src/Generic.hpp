/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <Eigen/Dense>
#include "H5Cpp.h"
#include <H5Group.h>
#include <complex>
#include <random>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <typeinfo>
#include <type_traits>
#include <complex>
#include <iostream>
#include <chrono>
#include <thread>
#include <cmath>
#include <math.h>
#include <initializer_list>

// Set of compilation parameters chosen in the Makefile
// MEMORY is the number of KPM vectors stored in the memory while calculating Gamma2D
// TILE is the size of the memory blocks used in the program
// COMPILE_MAIN is a flag to prevent compilation of unnecessary parts of the code when testing
#ifndef MEMORY
#define MEMORY 16
#endif

#ifndef TILE
#define TILE 64
#endif

#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef VERBOSE
#define VERBOSE 1
#endif

#ifndef ESTIMATE_TIME
#define ESTIMATE_TIME 0
#endif

#ifndef COMPILE_WAVEPACKET
#define COMPILE_WAVEPACKET 1
#endif

// other compilation parameters not set in the Makefile
// NGHOSTS is the extra length in each direction, to be used with the blocks of size TILE
#define PATTERNS  4
#define NGHOSTS   2
#define VVERBOSE 0
#define SSPRINT 0

// These are the verbose and debug messages
#define outcol "\033[1;31m"
#define outres "\033[0m"

#ifdef VERBOSE
#if VERBOSE==1
#define verbose_message(VAR)              \
  _Pragma("omp master")                   \
  {                                       \
    std::cout<<VAR<<std::flush;           \
  }                                       \
  _Pragma("omp barrier")
#else
#define verbose_message(VAR) 
#endif
#else
#define verbose_message(VAR) 
#endif


#ifdef VVERBOSE
#if VVERBOSE==1
#define vverbose_message(VAR)              \
  _Pragma("omp master")                   \
  {                                       \
    std::cout<<VAR<<std::flush;           \
  }                                       \
  _Pragma("omp barrier")
#else
#define vverbose_message(VAR) 
#endif
#else
#define vverbose_message(VAR) 
#endif

#ifdef DEBUG
#if DEBUG==1
#define debug_message(VAR) std::cout<<VAR<<std::flush; 
#else
#define debug_message(VAR) 
#endif
#else
#define debug_message(VAR) 
#endif
