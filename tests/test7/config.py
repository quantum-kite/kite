import kite
import matplotlib.pyplot as plt
import numpy as np
import pybinding as pb

def square_lattice(onsite=(0, 0)):
    a1 = np.array([1, 0])
    a2 = np.array([0, 1])
    lat = pb.Lattice(a1=a1, a2=a2)
    lat.add_sublattices(('A', [0, 0], onsite[0]))
    lat.add_hoppings(
        ([1, 0], 'A', 'A', - 1),
        ([0, 1], 'A', 'A', - 1)
    )

    return lat

nx = ny = 2
lx = ly = 256
lattice = square_lattice()

configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=True, precision=1, spectrum_range=[-4.1,4.1])

b1, b2 = lattice.reciprocal_vectors()
Gamma = [0, 0]
K1 = b1[0:2] + b2[0:2]
points = [Gamma, K1]

k_path = pb.results.make_path(*points, step=1.0)

calculation_arpes = kite.Calculation(configuration)
calculation_arpes.arpes(k_vector=k_path, num_moments=16, num_disorder=1)
kite.config_system(lattice, configuration, calculation_arpes, filename='config.h5')
