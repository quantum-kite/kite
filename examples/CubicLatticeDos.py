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
""" Square lattice

    Lattice : Square 1[nm] interatomic distance and t=1[eV] hopping;
    Configuration : size of the system 512x512, with domain decomposition (nx=ny=2), periodic boundary conditions,
                    double precision, automatic scaling;
    Calculation : dos;

"""

import kite
import numpy as np
import pybinding as pb


def cubic_lattice(onsite=(0, 0, 0)):
    """Make a square lattice with nearest neighbor hopping"""
    a1 = np.array([1, 0, 0])
    a2 = np.array([0, 1, 0])
    a3 = np.array([0, 0, 1])

    # create a lattice with 3 primitive vectors
    lat = pb.Lattice(a1=a1, a2=a2, a3=a3)

    # Add sublattices
    lat.add_sublattices(('A', [0, 0, 0], onsite[0]))

    # Add hoppings
    lat.add_hoppings(([1, 0, 1], 'A', 'A', - 1),
                     ([0, 1, 0], 'A', 'A', - 1),
                     ([0, 0, 1], 'A', 'A', - 1))
    return lat


# load a square lattice
lattice = cubic_lattice()
# number of decomposition parts in each direction of matrix. This divides the lattice into various sections,
# each of which is calculated in parallel
nx = ny = nz = 1
# number of unit cells in each direction.
lx = ly = lz = 32
# Twists in Each Direction
twsx = twsy = twsz = np.pi/2.0
# Boundary Mode
mode = "random"

configuration = kite.Configuration(divisions=[nx, ny, nz], length=[lx, ly, lz], boundaries=[mode,mode,mode], is_complex=True, precision=1,angles=(twsx, twsy, twsz))
# require the calculation of DOS
calculation = kite.Calculation(configuration)
calculation.dos(num_points=4000, num_moments=256, num_random=256, num_disorder=1)
# configure the *.h5 file
kite.config_system(lattice, configuration, calculation, filename='cubic_lattice.h5')