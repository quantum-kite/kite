import matplotlib.pyplot as plt
import export_lattice as ex
import numpy as np
import pybinding as pb

from export_lattice import Configuration, Calculation, Modification, export_lattice

energy_scale = 4.02
def graphene_initial(onsite=(0, 0)):
    a1 = np.array([1,0])
    a2 = np.array([0, 1])

    lat = pb.Lattice(
        a1=a1,
        a2=a2
    )

    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, 0], onsite[0])
    )

    lat.add_hoppings(
        ([1, 0], 'A', 'A', - 1),
        ([0, 1], 'A', 'A', - 1)
    )

    return lat

lattice = graphene_initial()

nx = ny = 2
lx = 512
ly = 512

configuration = Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                 is_complex=False, precision=1, energy_scale=energy_scale)

calculation = Calculation(configuration)
calculation.dos(num_points=1000, num_moments=256, num_random=1, num_disorder=1)

modification = ex.Modification(magnetic_field=False)
export_lattice(lattice, configuration, calculation, modification, 'example3.h5')


