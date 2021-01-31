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
import sys

from pybinding.repository.graphene import t_nn
t_nnn = 0.1 * t_nn
t_nnnn = 0.1 * t_nn
t = 2.8  # eV


def graphene(onsite=(0, 0), nearest_neighbors=1):
    """Make a honeycomb lattice with nearest neighbor hopping

    Parameters
    ----------
    onsite : tuple or list
        Onsite energy at different sublattices.
    """

    theta = np.pi / 3

    a1 = np.array([1 + np.cos(theta), np.sin(theta)])
    a2 = np.array([0, 2 * np.sin(theta)])
    lat = pb.Lattice(
        a1=a1, a2=a2
    )
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, 0], onsite[0]),
        ('B', [1, 0], onsite[1])
    )
    lat.add_hoppings(
        ([0, 0], 'A', 'B', - t),
        ([-1, 0], 'A', 'B', - t),
        ([-1, 1], 'A', 'B', - t)
    )

    if nearest_neighbors >= 2:
        lat.add_hoppings(
            ([0, -1], 'A', 'A', t_nn),
            ([0, -1], 'B', 'B', t_nn),
            ([1, -1], 'A', 'A', t_nn),
            ([1, -1], 'B', 'B', t_nn),
            ([1,  0], 'A', 'A', t_nn),
            ([1,  0], 'B', 'B', t_nn),
        )

    if nearest_neighbors >= 3:
        lat.add_hoppings(
            [(-2, +1), 'A', 'B', t_nnn],
            [(+0, -1), 'A', 'B', t_nnn],
            [(+0, +1), 'A', 'B', t_nnn],
        )

    if nearest_neighbors >= 4:
        lat.add_hoppings(
            ([-2, +1], 'A', 'A', t_nn),
            ([-2, +1], 'B', 'B', t_nn),
            ([1, +1], 'A', 'A', t_nn),
            ([1, +1], 'B', 'B', t_nn),
            ([1, -2], 'A', 'A', t_nn),
            ([1, -2], 'B', 'B', t_nn),
        )

    lat.min_neighbors = 2

    return lat


lattice = graphene(nearest_neighbors=3)

lx = int(sys.argv[1])
ly = int(sys.argv[2])

domain_decompose_1 = int(sys.argv[3])
domain_decompose_2 = int(sys.argv[4])

num_moments = int(sys.argv[5])

emin, emax = -3.01 * np.abs(t) - 6 * t_nn, 3.01 * np.abs(t) + 6 * t_nn
configuration = kite.Configuration(divisions=[domain_decompose_1, domain_decompose_2], length=[lx, ly],
                                   boundaries=[True, True], is_complex=False, precision=1, spectrum_range=[emin, emax])

calculation = kite.Calculation(configuration)
calculation.dos(num_points=1000, num_moments=num_moments, num_random=1, num_disorder=1)

name = 'dos_gr_nnn.h5'
kite.config_system(lattice, configuration, calculation, filename=name)
