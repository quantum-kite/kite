import kite
import numpy as np
import pybinding as pb


def square_lattice(onsite):
    """Make a square lattice with nearest neighbor hopping"""

    a1 = np.array([1, 0])
    a2 = np.array([0, 1])

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(
        a1=a1, a2=a2
    )

    # Add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, 0], onsite[0])
    )

    # Add hoppings
    lat.add_hoppings(
        ([1, 0], 'A', 'A', - 1),
        ([0, 1], 'A', 'A', - 1)
    )

    return lat


# load a square lattice
lattice = square_lattice([0])
# number of decomposition parts in each direction of matrix. This divides the lattice into various sections,
# each of which is calculated in parallel
nx = ny = 2
# number of unit cells in each direction.
lx = 512
ly = 256
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=False, precision=1, spectrum_range=[-4.1,4.1])

calculation = kite.Calculation(configuration)
calculation.dos(num_points=1000, num_moments=256, num_disorder=1, num_random=1)
kite.config_system(lattice, configuration, calculation, filename='config1.h5')
