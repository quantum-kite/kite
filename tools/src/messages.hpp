void print_header_message(){
#if VERBOSE
    std::cout << R"ENDSTRING(

+------------------------------------------------------------------------+
|            Chebyshev Polynomial Green's Function Approach              |
|             to Real-Space Quantum Transport Simulations                |
|                                                                        |
|                     KITE | Release 1.0                                 |
|                     Kite home: quantum-kite.com                        |
|                                                                        |
|    Created by Simao M. Joao, Joao V. Lopes (Universidade do Porto),    |
|      Tatiana G. Rappoport (Universidade Federal Rio de Janeiro),       |
|        Misa Andelkovic, Lucian Covaci (University of Antwerp)          |
|                and Aires Ferreira (University of York)                 |
|                                                                        |
|            Funded by The Royal Society | royalsociety.org              |
|                                                                        |
|  Copyright (C) 2018, 2019, M. Andelkovic, L. Covaci, A. Ferreira,      |
|                      S. M. Joao, J. V. Lopes, T. G. Rappoport          |
|                                                                        |
|  This program is free software: you can redistribute it and/or modify  |
|  it under the terms of the GNU General Public License as published by  |
|  the Free Software Foundation, either version 3 of the License, or     |
|  (at your option) any later version.                                   |
|                                                                        |
|  This program is distributed in the hope that it will be useful,       |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                  |
|  See the GNU General Public License for more details.                  |
+------------------------------------------------------------------------+
)ENDSTRING";
#endif
}

void print_flags_message(){
#if VERBOSE
    std::cout <<
        "------- FLAGS SET -------" << std::endl <<
        "DEBUG: "     << DEBUG    << std::endl <<
        "VERBOSE: "   << VERBOSE  <<  std::endl <<
        "-------------------------" << std::endl <<
	"" <<std::endl;
#endif
}

void print_info_message(){
#if VERBOSE
    std::cout << R"ENDSTRING(

------------------------------ INFORMATION ------------------------------
Linear response functions in units of e^2/h.
To stop these messages, set VERBOSE to 0 in the Makefile.
To see debug messages, set DEBUG to 1 in the Makefile.
-------------------------------------------------------------------------

)ENDSTRING";
#endif
}
