"""Basic lattice specification for monolayer graphene with nearest neighbor"""


import matplotlib.pyplot as plt
import pybinding as pb
import math


a = 0.24595  # [nm] unit cell length
a_cc = 0.142  # [nm] carbon-carbon distance
t = -2.8  # [eV] nearest neighbour hopping

lat = pb.Lattice(a1=[a, 0], a2=[a / 2, a / 2 * math.sqrt(3)])  # make a lattice object
lat.add_sublattices(('A', [0, -a_cc / 2]),  # add two unequivalent lattice sites
                    ('B', [0, a_cc / 2]))
lat.add_hoppings(
    # inside the main cell
    ([0, 0], 'A', 'B', t),
    # between neighboring cells
    ([1, -1], 'A', 'B', t),
    ([0, -1], 'A', 'B', t)
)

lat.plot()
plt.show()
