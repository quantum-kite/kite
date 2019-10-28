void print_header_message(){
#if VERBOSE
    std::cout << 
      "\n+------------------------------------------------------------------------+\n"
      "|            Chebyshev Polynomial Green's Function Approach              | \n"
      "|             to Real-Space Quantum Transport Simulations                | \n"    
      "|                                                                        | \n"                                              
      "|                     KITE | Release 1.0                                 | \n"         
      "|                     Kite home: quantum-kite.com                        | \n"
      "|                                                                        | \n"
      "|    Created by Simao M. Joao, Joao V. Lopes (Universidade do Porto),    | \n"
      "|      Tatiana G. Rappoport (Universidade Federal Rio de Janeiro),       | \n"
      "|        Misa Andelkovic, Lucian Covaci (University of Antwerp)          | \n"
      "|                and Aires Ferreira (University of York)                 | \n"
      "|                                                                        | \n"                                            
      "|            Funded by The Royal Society | royalsociety.org              | \n"
      "|                                                                        | \n"
      "|  Copyright (C) 2018, 2109, M. Andelkovic, L. Covaci, A. Ferreira,      | \n"
      "|                      S. M. Joao, J. V. Lopes, T. G. Rappoport          | \n"
      "|                                                                        | \n"
      "|  This program is free software: you can redistribute it and/or modify  | \n"
      "|  it under the terms of the GNU General Public License as published by  | \n"
      "|  the Free Software Foundation, either version 3 of the License, or     | \n"
      "|  (at your option) any later version.                                   | \n"
      "|                                                                        | \n"
      "|  This program is distributed in the hope that it will be useful,       | \n"
      "|  but WITHOUT ANY WARRANTY; without even the implied warranty of        | \n"
      "|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                  | \n"
      "|  See the GNU General Public License for more details.                  | \n"
      "+------------------------------------------------------------------------+\n";
#endif
}

void print_flags_message(){
#if VERBOSE
    std::cout <<
        "------- FLAGS SET -------\n"
        "DEBUG: "     << DEBUG    <<"\n"
        "VERBOSE: "   << VERBOSE  << "\n"
        "-------------------------\n\n";
#endif
}

void print_info_message(){
#if VERBOSE
    std::cout << 
        "\n------------------------------ INFORMATION ------------------------------\n"
        "Linear response functions in units of e^2/h.                              \n"
        "To stop these messages, set VERBOSE to 0 in the Makefile.                 \n"
        "To see debug messages, set DEBUG to 1 in the Makefile.                    \n"
        "------------------------------------------------------------------------- \n\n";
#endif
}
