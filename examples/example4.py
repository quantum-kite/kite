import kite_config as kite
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
        ([-1, 0], 'A', 'B', - 1),
        ([-1, 1], 'A', 'B', - 1)
    )

    node0 = [[+0, +0], 'A']
    node1 = [[+0, +0], 'B']
    struc_disorder_one = kite.StructuralDisorder(lat, concentration=1.00)
    struc_disorder_one.add_structural_disorder(
        (*node0, *node1, -1)
    )
    return lat, [struc_disorder_one]


lattice, disorded_structural = graphene_initial()

nx = ny = 1
lx = 256
ly = 256

configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=False, precision=1, spectrum_range=[-energy_scale, energy_scale])
calculation = kite.Calculation(configuration)
calculation.dos(num_points=1000, num_random=1, num_disorder=1, num_moments=256)

modification = kite.Modification(magnetic_field=False)
kite.export_lattice(lattice, configuration, calculation, modification, 'example4.h5',
                    disorded_structural=disorded_structural)
