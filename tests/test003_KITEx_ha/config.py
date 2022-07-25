import kite
import numpy as np
import pybinding as pb

def haldane(onsite=(0, 0), t=1):
    a = 0.24595  # [nm] unit cell length
    a_cc = 0.142  # [nm] carbon-carbon distance
    t2 = t/10

    a1 = a * np.array([a, 0])
    a2 = a * np.array([1 / 2, 1 / 2 * np.sqrt(3)])

    lat = pb.Lattice(a1=a1, a2=a2)
    lat.add_sublattices(
        ('A', [0, -a_cc/2], onsite[0]),
        ('B', [0,  a_cc/2], onsite[1])
    )

    lat.add_hoppings(
        ([0,  0], 'A', 'B', -t),
        ([1, -1], 'A', 'B', -t),
        ([0, -1], 'A', 'B', -t),
        ([1, 0], 'A', 'A', -t2 * 1j),
        ([0, -1], 'A', 'A', -t2 * 1j),
        ([-1, 1], 'A', 'A', -t2 * 1j),
        ([1, 0], 'B', 'B', -t2 * -1j),
        ([0, -1], 'B', 'B', -t2 * -1j),
        ([-1, 1], 'B', 'B', -t2 * -1j)
    )
    return lat

lattice = haldane()
nx = ny = 2
lx = ly = 64

mode = "periodic"
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[mode, mode], is_complex=True, precision=0, spectrum_range=[-10, 10])
calculation = kite.Calculation(configuration)
calculation.dos(num_points=1000, num_moments=64, num_random=1, num_disorder=1)
kite.config_system(lattice, configuration, calculation, filename="config.h5")
