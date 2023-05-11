import numpy as np
import pybinding as pb
import kite


def graphene(onsite=(0, 0)):
    theta = np.pi / 3
    t = 1  # eV
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

    return lat


lattice = graphene()
nx = ny = 2
lx = ly = 512
M = 128
NPoints = 1000

configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=["periodic", "periodic"], is_complex=False, precision=1, spectrum_range=[-5, 5])
calculation = kite.Calculation(configuration)
calculation.dos(num_points=NPoints, num_moments=M, num_random=1, num_disorder=1)
kite.config_system(lattice, configuration, calculation, filename='config_sym.h5')
