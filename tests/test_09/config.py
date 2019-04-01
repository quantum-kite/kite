import kite
import matplotlib.pyplot as plt
import numpy as np
import pybinding as pb

def graphene_initial(onsite=(0,0)):

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

nx = ny = 4
lx = ly = 256
lattice = graphene_initial()

configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=True, precision=1, spectrum_range=[-4.1,4.1])

b1, b2 = lattice.reciprocal_vectors()
Gamma = [0, 0]
M = [np.pi/np.sqrt(3), np.pi/3.0]
K = [4*np.pi/np.sqrt(27), 0.0]
moments = 32
points = [Gamma, M, K, Gamma]

k_path = pb.results.make_path(*points, step=0.1)

calculation_arpes = kite.Calculation(configuration)
calculation_arpes.arpes(k_vector=k_path, weight=[0.5,1.5], num_moments=moments, num_disorder=1)
kite.config_system(lattice, configuration, calculation_arpes, filename='config.h5')
