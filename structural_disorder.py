import matplotlib.pyplot as plt
import numpy as np
import pybinding as pb
from export_lattice import Configuration, Calculation, Modification, Disorder, StructuralDisorder, \
    export_lattice, make_pybinding_model

# define lattice of monolayer graphene with 1[nm] interatomic distance and t=1/3[eV] hopping,
# EnergyScale is the scaling factor of the hopping parameters, important for the rescaling of the spectral quantity.
#  INFO: other examples are defined in define_lattice.py script


def graphene_initial(onsite=(0, 0)):
    """Return the basic lattice specification for monolayer graphene with nearest neighbor"""

    theta = np.pi / 3
    t=2.8
    a1 = np.array([1 + np.cos(theta), np.sin(theta)])
    a2 = np.array([0, 2 * np.sin(theta)])

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(
        a1=a1,
        a2=a2
    )

    # Add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, 0], onsite[0]),
        ('B', [1, 0], onsite[1])
    )

    # Add hoppings
    lat.add_hoppings(
        # inside the main cell, between which atoms, and the value
        ([0, 0], 'A', 'B', - t),
        # between neighboring cells, between which atoms, and the value
        ([-1, 0], 'A', 'B', - t),
        ([-1, 1], 'A', 'B', - t)
    )

    # Add bond disorder as an object of a class StructuralDisorder. In this manner we can add onsite and bond defects
    # with a specific concentration, which will be added to the simulated system. The procedure for adding is same
    # as adding the hopping, with the difference that the bond disorded is not bounded to one site in the [0, 0]
    # unit cell.
    node0 = [[+0, +0], 'A']
    node1 = [[+0, +0], 'B']
    node2 = [[+1, +0], 'A']
    node3 = [[+0, +1], 'B']
    node4 = [[+0, +1], 'A']
    node5 = [[-1, +1], 'B']

    struc_disorder_one = StructuralDisorder(lat, concentration=0.05)
    struc_disorder_one.add_structural_disorder(
        # add bond disorder in the form [from unit cell], 'sublattice_from', [to_unit_cell], 'sublattice_to', value:
        (*node0, *node1, 1),
        (*node1, *node2, 1),
        (*node2, *node3, 1),
        (*node3, *node4, 1),
        (*node4, *node5, 1),
        (*node5, *node0, 1),
        # in this way we can add onsite disorder in the form [unit cell], 'sublattice', value
        ([+0, +0], 'B', 0.3)
    )
    # It is possible to add multiple different disorder type which should be forwarded to the export_lattice function
    # as a list.
    struc_disorder_two = StructuralDisorder(lat, concentration=0.2)
    struc_disorder_two.add_structural_disorder(
        (*node0, *node1, 0.4),
        (*node4, *node5, 0.4),
        (*node5, *node0, 0.4),
        ([+0, +0], 'B', 0.4)
    )
    struc_disorder_two.add_vacancy('B')

    struc_disorder_three = StructuralDisorder(lat, concentration=0.01)
    struc_disorder_three.add_vacancy('A')

    # if there is disorder it should be returned separately from the lattice
    return lat, [struc_disorder_one, struc_disorder_two, struc_disorder_three]


lattice, disorded_structural = graphene_initial()
# number of decomposition parts in each direction of matrix.
# This divides the lattice into various sections, each of which is calculated in parallel
nx = ny = 1

# number of unit cells in each direction.
lx = 512
ly = 512

configuration = Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                              is_complex=False, precision=1,spectrum_range=[-10, 10])

calculation = Calculation(configuration)
calculation.dos(num_moments=1024, num_random=1, num_disorder=1, num_points=1000)

modification = Modification(magnetic_field=False)

export_lattice(lattice, configuration, calculation, modification, 'structural_disorder.h5',
               disorded_structural=disorded_structural)

# plotting the lattice
lattice.plot()
# plt.show()
