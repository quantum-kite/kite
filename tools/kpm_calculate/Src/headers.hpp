#include <stdio.h>		
#include <iostream>
#include <fstream>
#include <functional>
#include <cmath>
#include <math.h>
#include <Eigen/Dense>
#include <complex>
#include <string>

#define outcol "\033[1;31m"
#define outres "\033[0m"

#ifdef VERBOSE
	#if VERBOSE==1
		#define verbose_message(VAR) std::cout<<VAR<<std::flush
	#else
		#define verbose_message(VAR) 
	#endif
#else
	#define verbose_message(VAR) 
#endif

#ifdef DEBUG
	#if DEBUG==1
		#define debug_message(VAR) std::cout << outcol << VAR << outres << std::flush
		#define verbose_message(VAR) std::cout<<VAR<<std::flush
	#else
		#define debug_message(VAR) 
	#endif
#else
	#define debug_message(VAR) 
#endif
