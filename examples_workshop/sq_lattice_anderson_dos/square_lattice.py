"""       
        ##############################################################################      
        #                        KITE | Release  1.0                                 #      
        #                                                                            #      
        #                        Kite home: quantum-kite.com                         #           
        #                                                                            #      
        #  Developed by: Simao M. Joao, Joao V. Lopes, Tatiana G. Rappoport,         #       
        #  Misa Andelkovic, Lucian Covaci, Aires Ferreira, 2018-2019                 #      
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


def square_lattice(onsite=(0, 0)):
    """Make a square lattice with nearest neighbor hopping"""

    a1 = np.array([1, 0])
    a2 = np.array([0, 1])

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(
        a1=a1, a2=a2
    )

    # Add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, 0], onsite[0])
    )

    # Add hoppings
    lat.add_hoppings(
        ([1, 0], 'A', 'A', - 1),
        ([0, 1], 'A', 'A', - 1)
    )

    return lat


# load a square lattice
lattice = square_lattice()
disorder = kite.Disorder(lattice)
disorder.add_disorder('A', 'Uniform', 0.0, 0.1)


# number of decomposition parts in each direction of matrix. This divides the lattice into various sections,
# each of which is calculated in parallel
nx = ny = 1
# number of unit cells in each direction.
lx = ly = 64
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=False, precision=1, spectrum_range=[-5.0,5.0])
# require the calculation of DOS
calculation = kite.Calculation(configuration)
calculation.dos(num_points=8192*4, num_moments=4*4096, num_random=1, num_disorder=1)
# configure the *.h5 file
kite.config_system(lattice, configuration, calculation, filename='square_lattice.h5', disorder=disorder)
