import kite
import numpy as np
import pybinding as pb

# define lattice of monolayer graphene with 1[nm] interatomic distance and t=1/3[eV] hopping,
# EnergyScale is the scaling factor of the hopping parameters, important for the rescaling of the spectral quantity.
#  INFO: other examples are defined in define_lattice.py script
energy_scale = 3.06


def graphene_initial(onsite=(0, 0)):
    """Return the basic lattice specification for monolayer graphene with nearest neighbor"""

    theta = np.pi / 3
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
        ([0, 0], 'A', 'B', - 1),
        # between neighboring cells, between which atoms, and the value
        ([-1, 0], 'A', 'B', - 1),
        ([-1, 1], 'A', 'B', - 1)
    )

    return lat


lattice = graphene_initial()

# number of decomposition parts in each direction of matrix.
# This divides the lattice into various sections, each of which is calculated in parallel
nx = ny = 1

# number of unit cells in each direction.
lx = 256
ly = 256

configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=True, precision=1, spectrum_range=[-energy_scale, energy_scale])

calculation = kite.Calculation(configuration)
calculation.dos(num_moments=1024, num_random=1, num_disorder=1, num_points=1000)

modification = kite.Modification(magnetic_field=False)

kite.config_system(lattice, configuration, calculation, modification, 'example8.h5')

# plotting the lattice
lattice.plot()
# plt.show()
