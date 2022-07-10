"""       
        ##############################################################################      
        #                        KITE | Release  1.1                                 #      
        #                                                                            #      
        #                        Kite home: quantum-kite.com                         #           
        #                                                                            #      
        #  Developed by: Simao M. Joao, Joao V. Lopes, Tatiana G. Rappoport,         #       
        #  Misa Andelkovic, Lucian Covaci, Aires Ferreira, 2018-2022                 #      
        #                                                                            #      
        ##############################################################################      
"""
""" Density of states of pristine Graphene

    Units: Energy in eV
    Lattice: Honeycomb
    Configuration: L_x=L_y=32, periodic boundary conditions, double precision, automatic Hamiltonian rescaling;
    Calculation type: Average DOS;
    Last updated: 10/07/2022
"""

import kite
import numpy as np
import pybinding as pb
import os


def graphene_lattice(onsite=(0, 0)):
    # Return lattice specification for a honyecomb lattice with nearest neighbor hoppings

    theta = np.pi/3.
    t = 2.8  # eV
    a1 = np.array([1 + np.cos(theta), np.sin(theta)])
    a2 = np.array([0, 2 * np.sin(theta)])
    lat = pb.Lattice(a1=a1, a2=a2)
    
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, 0], onsite[0]),
        ('B', [1, 0], onsite[1])
    )
    lat.add_hoppings(([0, 0], 'A', 'B', - t),
                     ([-1, 0], 'A', 'B', - t),
                     ([-1, 1], 'A', 'B', - t))

    return lat

# load lattice
lattice = graphene_lattice()
# number of decomposition parts in each direction of matrix. This divides the lattice into various sections,
# each of which is calculated in parallel
nx = ny = 1
# number of unit cells in each direction.
lx = ly = 32

# - boundary conditions [mode,mode, ... ] with modes:
#   . "periodic"
#   . "open"
#   . "twist_fixed"     this option needs the extra argument ths=[phi_1,..,phi_DIM] where phi_i \in [0, 2*M_PI]  
#   . "twist_random"
# Boundary Mode
mode = "periodic"

configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[mode,mode], is_complex=True, precision=1)
# require the calculation of DOS
calculation = kite.Calculation(configuration)
calculation.dos(num_points=4000, num_moments=256, num_random=256, num_disorder=1)
# configure the *.h5 file
#os.chdir("ConfigurationFiles/")
kite.config_system(lattice, configuration, calculation, filename='graphene_lattice.h5')
