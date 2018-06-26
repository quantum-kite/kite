"""       
        ##############################################################################      
        #                        KITE | Pre-Release version 0.1                      #      
        #                                                                            #      
        #                        Kite home: quantum-kite.com                         #           
        #                                                                            #      
        #  Developed by: Simao M. Joao, Joao V. Lopes, Tatiana G. Rappoport,         #       
        #  Misa Andelkovic, Lucian Covaci, Aires Ferreira, 2018                      #      
        #                                                                            #      
        ##############################################################################      
"""
""" Honeycomb lattice with vacancy disorder

    Lattice : Honeycomb 1[nm] interatomic distance and t=1[eV] hopping;
    Disorder : StructuralDisorder, vacancy with concentration 0.1 inside A and 0.05 inside B sublattices;
    Configuration : size of the system 512x512, without domain decomposition (nx=ny=1), periodic boundary conditions,
                    double precision, automatic scaling;
    Calculation : dos;

"""

import kite
import numpy as np
import pybinding as pb


def honeycomb_lattice(onsite=(0, 0)):
    """Make a honeycomb lattice with nearest neighbor hopping
    
    Parameters
    ----------
    onsite : tuple or list
        Onsite energy at different sublattices.
    """""

    # define lattice vectors
    theta = np.pi / 3
    a1 = np.array([1 + np.cos(theta), np.sin(theta)])
    a2 = np.array([0, 2 * np.sin(theta)])

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(
        a1=a1, a2=a2
    )

    # Add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, 0], onsite[0]),
        ('B', [1, 0], onsite[1])
    )

    # Add hoppings
    lat.add_hoppings(
        # inside the main cell, between which atoms, and the value
        ([+0, +0], 'A', 'B', - 1),
        # between neighboring cells, between which atoms, and the value
        ([-1, +0], 'A', 'B', - 1),
        ([-1, +1], 'A', 'B', - 1)
    )

    return lat


lattice = honeycomb_lattice((-0.0, 0.0))

# Add vacancy disorder as an object of a class StructuralDisorder. In this manner we can distribute vacancy disorder
# on a specific sublattice with a specific concentration.
# unless you would like the same pattern of disorder at both sublatices,
# each realisation should be specified as a separate object
struc_disorder_A = kite.StructuralDisorder(lattice, concentration=0.1)
struc_disorder_A.add_vacancy('A')
struc_disorder_B = kite.StructuralDisorder(lattice, concentration=0.1)
struc_disorder_B.add_vacancy('B')
disorder_structural = [struc_disorder_A, struc_disorder_B]
# load a honeycomb lattice and structural disorder

nx = ny = 1
lx = 512
ly = 512
# make config object which caries info about
# - the number of decomposition parts [nx, ny],
# - lengths of structure [lx, ly]
# - boundary conditions, setting True as periodic boundary conditions, and False elsewise,
# - info if the exported hopping and onsite data should be complex,
# - info of the precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
# - scaling, if None it's automatic, if present select spectrum_bound=[e_min, e_max]
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=False, precision=1)
# require the calculation of dos
num_moments = 512
calculation = kite.Calculation(configuration)
calculation.dos(num_points=1000, num_moments=num_moments, num_random=5, num_disorder=1)
# configure the *.h5 file
kite.config_system(lattice, configuration, calculation, filename='vacancies.h5',
                   disorder_structural=disorder_structural)
