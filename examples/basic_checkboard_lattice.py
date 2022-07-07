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
""" checkboard lattice

    Lattice : Square 1[nm] interatomic distance and t=1[eV] hopping;
    Configuration : size of the system 512x512, with domain decomposition (nx=ny=2), periodic boundary conditions,
                    double precision, automatic scaling;
    Calculation : dos;

"""

import kite
import numpy as np
import pybinding as pb

def checkboard_lattice(onsite=(0, 0)):
    """Return the lattice specification for monolayer graphene"""
    a = 1   # [nm] unit cell length
    t = 1      # [eV] nearest neighbour hopping
    delta=0
    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(
        a1=[a, 0],
        a2=[0, a]
    )

    lat.add_sublattices(
        # name and position
        ('A', [0,  0], onsite[0]),
        ('B', [a/2,  a/2], onsite[1])
    )

    lat.add_hoppings(
        # inside the main cell
        ([0,  0], 'A', 'B',t),
        # between neighboring cells
        ([-1, 0], 'A', 'B',  t),
        ([0, -1], 'A', 'B',  t),
        ([-1, -1], 'A', 'B', t)
    )

    return lat



# load a square lattice
lattice = checkboard_lattice((-0.1,0.1))
# number of decomposition parts in each direction of matrix. This divides the lattice into various sections,
# each of which is calculated in parallel
nx = ny = 2
# number of unit cells in each direction.
lx = ly = 512
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=False, precision=1)
# require the calculation of DOS
calculation = kite.Calculation(configuration)
calculation.dos(num_points=1000, num_moments=512, num_random=5, num_disorder=1)
# configure the *.h5 file
kite.config_system(lattice, configuration, calculation, filename='checkboard_lattice.h5')
