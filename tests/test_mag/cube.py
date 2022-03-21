"""       
        ##############################################################################      
        #                        KITE | Release  1.0                                 #      
        #                                                                            #      
        #                        Kite home: quantum-kite.com                         #           
        #                                                                            #      
        #  Developed by: Simao M. Joao, Joao V. Lopes, Tatiana G. Rappoport,         #       
        #  Misa Andelkovic, Lucian Covaci, Aires Ferreira, 2018-2021                 #      
        #                                                                            #      
        ##############################################################################      
"""
import numpy as np
import pybinding as pb
import kite


def cube(onsite=(0, 0)):

    t = 1
    a1 = np.array([1, 0, 0])
    a2 = np.array([0, 1, 0])
    a3 = np.array([0, 0, 1])
    lat = pb.Lattice( a1=a1, a2=a2, a3=a3)
    lat.add_sublattices( ('A', [0, 0, 0], onsite[0]))
    lat.add_hoppings(
        ([1, 0, 0], 'A', 'A', - t),
        ([0, 1, 0], 'A', 'A', - t),
        ([0, 0, 1], 'A', 'A', - t))

    return lat


# make a graphene lattice
lattice = cube()


# number of decomposition parts in each direction of matrix.
# This divides the lattice into various sections, each of which is calculated in parallel
nx = 2
ny = 1
nz = 2
# number of unit cells in each direction.
lx = 64
ly = 64
lz = 64
# lx = ly = lz = 64
configuration = kite.Configuration(divisions=[nx, ny, nz], length=[lx, ly, lz], boundaries=[True, True, True], is_complex=True, precision=1, spectrum_range=[-6.1,6.1])
calculation = kite.Calculation(configuration)
mod = kite.Modification(magnetic_field = 1)
calculation.dos(num_points=1000, num_moments=32768, num_random=1, num_disorder=1)
# configure the *.h5 file
# kite.config_system(lattice, configuration, calculation, filename='cube.h5')

kite.config_system(lattice, configuration, calculation, modification=mod, filename='cube.h5')
