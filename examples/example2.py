import matplotlib.pyplot as plt
import kite
import numpy as np
import pybinding as pb

energy_scale = 3.06


def graphene_initial(onsite=(0, 0)):
    theta = np.pi / 3
    a1 = np.array([1 + np.cos(theta), np.sin(theta)])
    a2 = np.array([0, 2 * np.sin(theta)])
    lat = pb.Lattice(
        a1=a1,
        a2=a2
    )
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, 0], onsite[0]),
        ('B', [1, 0], onsite[1])
    )
    lat.add_hoppings(
        ([0, 0], 'A', 'B', - 1),
        ([-1, 0], 'A', 'B', - 1),
        ([-1, 1], 'A', 'B', - 1)
    )

    return lat


lattice = graphene_initial()

nx = ny = 1

lx = 128
ly = 128

configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=False, precision=1, spectrum_range=[-energy_scale, energy_scale])

calculation = kite.Calculation(configuration)
calculation.conductivity_optical(num_points=1000, num_disorder=1, num_random=1, num_moments=128, direction='xx')

modification = kite.Modification(magnetic_field=False)

kite.config_system(lattice, configuration, calculation, modification, 'example2.h5')
