L=256
N=128
S=1
R=1
W=0.1

import numpy as np
import pybinding as pb
import kite

energy_scale = 4.0

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


disorder = kite.Disorder(lattice)
disorder.add_disorder('B', 'Uniform', 0.0, W*0.5/np.sqrt(3))
disorder.add_disorder('A', 'Uniform', 0.0, W*0.5/np.sqrt(3))

nx = ny = 2

configuration = kite.Configuration(divisions=[nx, ny], length=[L, L], boundaries=[True, True], is_complex=False, precision=0, spectrum_range=[-energy_scale, energy_scale]) 


calculation = kite.Calculation(configuration) 
calculation.singleshot_conductivity_dc(energy=[(n/100.0 - 0.5)*2 for n in range(101)], num_moments=256, num_random=1, num_disorder=1,direction='xx', eta=0.02)
kite.config_system(lattice, configuration, calculation, disorder=disorder, filename='config.h5')
