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
""" Cubic lattice with twisted boundary conditions

    Units: arbitrary (energy in units of hopping, |t| = 1)
    Lattice: Simple cubic lattice
    Configuration: L_x=L_y=L_z=32, twisted boundary conditions, double precision, automatic Hamiltonian rescaling;
    Calculation type: Average DOS;
    Last updated: 10/07/2022
"""

import kite
import numpy as np
import pybinding as pb

def cubic_lattice(onsite=(0, 0, 0)):
    # Return lattice specification for a cubic lattice with nearest neighbor hoppings
    a1 = np.array([1, 0, 0])
    a2 = np.array([0, 1, 0])
    a3 = np.array([0, 0, 1])

    # create lattice with 3 primitive vectors
    lat = pb.Lattice(a1=a1, a2=a2, a3=a3)

    # Add sublattices
    lat.add_sublattices(('A', [0, 0, 0], onsite[0]))

    # Add hoppings
    lat.add_hoppings(([1, 0, 1], 'A', 'A', - 1),
                     ([0, 1, 0], 'A', 'A', - 1),
                     ([0, 0, 1], 'A', 'A', - 1))
    return lat

# load lattice
lattice = cubic_lattice()
# number of decomposition parts [nx,ny,nz] in each direction of matrix. This divides the lattice into various sections,
# each of which is calculated in parallel
nx = ny = nz = 1
# number of unit cells in each direction.
lx = ly = lz = 32

# - boundary conditions [mode,mode, ... ] with modes:
#   . "periodic"
#   . "open"
#   . "twist_fixed" -- this option needs the extra argument ths=[phi_1,..,phi_DIM] where phi_i \in [0, 2*M_PI]  
#   . "twist_random"
# Boundary Mode
mode = "random"
# Twists in each direction
twsx = twsy = twsz = np.pi/2.0

# - specify precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
# - scaling, if None it's automatic, if present select spectrum_range=[e_min, e_max]

configuration = kite.Configuration(divisions=[nx, ny, nz], length=[lx, ly, lz], boundaries=[mode,mode,mode], is_complex=True, precision=1,angles=(twsx, twsy, twsz))
# specify calculation type
calculation = kite.Calculation(configuration)
calculation.dos(num_points=4000, num_moments=256, num_random=256, num_disorder=1)
# configure the *.h5 file
kite.config_system(lattice, configuration, calculation, filename='cubic_lattice.h5')
