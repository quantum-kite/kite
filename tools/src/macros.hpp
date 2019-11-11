#define TEXTCOLOR_RED "\033[1;31m"
#define TEXTCOLOR_RESET "\033[0m"

#ifndef DEBUG
#define DEBUG 1
#endif

#ifndef VERBOSE
#define VERBOSE 1
#endif

#define PI 3.14159265
// The first scale is 1/2pi and corresponds to the conductance quantum per spin e^2/h
// The second scale is the universal conductivity of graphene e^2/4h_bar
#define scale1 1/(2*PI)
#define scale2 0.25
#define unit_scale 1/(2*PI)

#if VERBOSE==1
	#define verbose_message(VAR) std::cout<<VAR<<std::flush
#else
	#define verbose_message(VAR) 
#endif

#if DEBUG==1
	#define debug_message(VAR) std::cout << TEXTCOLOR_RED << VAR << TEXTCOLOR_RESET << std::flush
	#define verbose_message(VAR) std::cout<<VAR<<std::flush
#else
	#define debug_message(VAR) 
#endif
