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

import numpy as np
import pybinding as pb
import kite
import sys


def graphene_initial(onsite=(0.0, 0.0)):

    theta = np.pi / 3
    a1 = np.array([2 * np.sin(theta), 0])
    a2 = np.array([np.sin(theta), 1 + np.cos(theta)])

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(
    a1=a1,
    a2=a2
    )

    # Add sublattices
    lat.add_sublattices(
    # name, position, and onsite potential
    ('A', [0, 0], onsite[0]),
    ('B', [0, 1], onsite[1])
    )

    # Add hoppings
    lat.add_hoppings(
    # inside the main cell, between which atoms, and the value
    ([0, 0], 'A', 'B', - 1),
    # between neighboring cells, between which atoms, and the value
    ([0, -1], 'A', 'B', - 1),
    ([1, -1], 'A', 'B', - 1)
    )

    return lat


lattice = graphene_initial()

divX        = 1
divY        = 1
sizeX       = 256
sizeY       = 256
NumMoments  = 512

energy_scale = 3.1
configuration = kite.Configuration(
    divisions=[divX, divY], 
    length=[sizeX, sizeY], 
    boundaries=[True, True], 
    is_complex=False, 
    precision=0, 
    spectrum_range=[-energy_scale, energy_scale]) 

calculation = kite.Calculation(configuration) 

calculation.conductivity_optical(
    num_points=256, 
    num_moments=NumMoments, 
    num_random=1, 
    num_disorder=1, 
    direction="xx", 
    temperature=0.01) 

kite.config_system(lattice, configuration, calculation, filename="graphene_condopt.h5")
