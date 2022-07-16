import kite
import numpy as np
import pybinding as pb

def square_lattice(onsite=[0,0]):
    a1 = np.array([1.1, 0.1])
    a2 = np.array([0.3, 1.0])
    lat = pb.Lattice( a1=a1, a2=a2)
    lat.add_sublattices( ('A', [0, 0], onsite[0]))
    lat.add_sublattices( ('B', [0, 0.1], onsite[1]))

    # Hoppings within sublattice A
    lat.add_hoppings(([1, 0], 'A', 'A', - 1),
                     ([0, 1], 'A', 'A', - 1))

    # Hoppings within sublattice B
    lat.add_hoppings(([1, 0], 'B', 'B', - 1),
                     ([0, 1], 'B', 'B', - 1))

    # Hoppings across sublattices
    lat.add_hoppings(([0, 0], 'A', 'B', 0.1))

    return lat

lattice = square_lattice()
nx = ny = 2
lx = ly = 512

a,b = -6,6
M = 2048

configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=["periodic", "periodic"], is_complex=False, precision=1, spectrum_range=[a,b], custom_local=True, custom_local_print=True)
calculation = kite.Calculation(configuration)
calculation.dos(num_points=1000, num_moments=M, num_random=1, num_disorder=1)
kite.config_system(lattice, configuration, calculation, filename='config.h5')
