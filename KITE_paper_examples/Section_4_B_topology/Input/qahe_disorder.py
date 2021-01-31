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
""" Honeycomb lattice with Kane-Mele, Rashba SOC model and Exchange """

import numpy as np
import pybinding as pb
import matplotlib.pyplot as plt
import kite
from pybinding.repository.graphene import a_cc, a
from pybinding.constants import hbar
from math import sqrt, pi

el_charge = 1.602 * 10 ** -19  #: [C] electron charge
t = -2.507  #: [eV] graphene nearest neighbor hopping
vf = 1.5 * a_cc * np.abs(t) / hbar * 10 ** -9  #: [m/s] Fermi velocity

lambda_R = 0.3*t  #: [eV] Rashba SOC in units of t
lambda_I = 0.0  #: [eV] Intrinsic SOC in units of t
exch = 0.4*t   #: [eV] Exchange in units of t
# allow sublatice dependent intrinsic SOC (spin valley)
lambda_I_A = lambda_I
lambda_I_B = lambda_I

rashba_so = lambda_R * 2.0 * 1.0j / 3.0  #: [eV] constant and geometrical factors that will define Rashba SOC

km_so = lambda_I * 1j / (3 * sqrt(3))  #: [eV] constant and geometrical factors that will define intrinsic SOC

km_so_A = lambda_I_A * 1j / (3 * sqrt(3))  #: [eV] constant and geometrical factors that will define intrinsic SOC
km_so_B = lambda_I_B * 1j / (3 * sqrt(3))  #: [eV] constant and geometrical factors that will define intrinsic SOC

a1 = np.array([+a * sqrt(3) / 2, +a / 2, 0])  #: [nm] unit cell vectors graphene
a2 = np.array([+a * sqrt(3) / 2, -a / 2, 0])

posA = np.array([-a_cc / 2, 0])
posB = np.array([+a_cc / 2, 0])

# delta1 = [acc, 0]
# delta2 = [-0.5 * acc, +sqrt(3) / 2 * acc]
# delta3 = [-0.5 * acc, -sqrt(3) / 2 * acc]

def honeycomb_lattice_with_SO(onsite=(0, 0)):
    """Make a honeycomb lattice with nearest neighbor hopping with SOC

    Parameters
    ----------
    onsite : tuple or list
        Onsite energy at different sublattices.
    """

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(
        a1=a1, a2=a2
    )

    # Add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('Aup', posA, onsite[0]+exch),
        ('Bup', posB, onsite[1]+exch),
        ('Adown', posA, onsite[0]-exch),
        ('Bdown', posB, onsite[1]-exch)
    )

    # Add hoppings-
    lat.add_hoppings(
        # ([f - i ], i , f )
        # inside the main cell, between which atoms, and the value
        ([0, 0], 'Adown', 'Bdown', t),
        ([0, 0], 'Aup', 'Bup', t),
        # between neighboring cells, between which atoms, and the value
        ([0, -1], 'Aup', 'Bup', t),
        ([0, -1], 'Adown', 'Bdown', t),

        ([-1, 0], 'Aup', 'Bup', t),
        ([-1, 0], 'Adown', 'Bdown', t)
    )

    if np.abs(lambda_R) > 0:
        lat.add_hoppings(
            # Rashba nearest neighbor, spin flip
            # inside the main cell, between which atoms, and the value
            ([0, 0], 'Aup', 'Bdown', 1j * rashba_so),  # delta1
            ([0, -1], 'Aup', 'Bdown', (+sqrt(3) / 2 - 0.5 * 1j) * rashba_so),  # delta2
            ([-1, 0], 'Aup', 'Bdown', (-sqrt(3) / 2 - 0.5 * 1j) * rashba_so),  # delta3

            ([0, 0], 'Adown', 'Bup', -1j * rashba_so),  # delta1
            ([0, -1], 'Adown', 'Bup', (+sqrt(3) / 2 + 0.5 * 1j) * rashba_so),  # delta2
            ([-1, 0], 'Adown', 'Bup', (-sqrt(3) / 2 + 0.5 * 1j) * rashba_so)  # delta3
        )

    if np.abs(lambda_I) > 0:
        # Kane-Mele SOC, same spin next-nearest
        # between neighboring cells, between which atoms, and the value
        lat.add_hoppings(
            ([0, -1], 'Aup', 'Aup', +km_so_A),
            ([0, -1], 'Adown', 'Adown', -km_so_A),

            ([-1, 0], 'Aup', 'Aup', -km_so_A),
            ([-1, 0], 'Adown', 'Adown', +km_so_A),

            ([1, -1], 'Aup', 'Aup', -km_so_A),
            ([1, -1], 'Adown', 'Adown', +km_so_A),

            ([0, -1], 'Bup', 'Bup', -km_so_B),
            ([0, -1], 'Bdown', 'Bdown', +km_so_B),

            ([-1, 0], 'Bup', 'Bup', +km_so_B),
            ([-1, 0], 'Bdown', 'Bdown', -km_so_B),

            ([1, -1], 'Bup', 'Bup', +km_so_B),
            ([1, -1], 'Bdown', 'Bdown', -km_so_B),
        )

    return lat

lattice = honeycomb_lattice_with_SO()

# Define disorder

disorder = kite.Disorder(lattice)
disorder.add_disorder(['Aup', 'Adown'], 'Uniform', +0.0, 0.1)
disorder.add_disorder(['Bup', 'Bdown'], 'Uniform', +0.0, 0.1)

# Simulation for a small system ~ 10^7 orbitals

# number of decomposition parts in each direction of matrix.
nx = 4
ny = 4 # Assumes nx * ny = 16 cores 

# number of unit cells in each direction.
lx = 2048
ly = 2048

#
# make config object which caries info about
# - the number of decomposition parts [nx, ny],
# - lengths of structure [lx, ly]
# - boundary conditions, setting True as periodic boundary conditions, and False elsewise,
# - info if the exported hopping and onsite data should be complex,
# - info of the precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
# - scaling, if None it's automatic, if present select spectrum_bound=[e_min, e_max]
#

configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=True, precision=1)
calculation = kite.Calculation(configuration)

# DC conductivity
calculation.conductivity_dc(num_points=6000, num_moments=2048, num_random=1, num_disorder=1,
                            direction='xy', temperature=1)

# Optical conductivity 
calculation.conductivity_optical(num_points=200, num_moments=2048, num_random=1, num_disorder=1,
                            direction='xx')

# configure the *.h5 file
kite.config_system(lattice, configuration, calculation, filename='qahe_disorder.h5', disorder=disorder)
