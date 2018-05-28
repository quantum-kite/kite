import matplotlib.pyplot as plt
import export_lattice as ex
import numpy as np
import pybinding as pb
import sys
from export_lattice import Configuration, Calculation, Modification, export_lattice

# define lattice of monolayer graphene with 1[nm] interatomic distance and t=1/3[eV] hopping,
# EnergyScale is the scaling factor of the hopping parameters, important for the rescaling of the spectral quantity.
#  INFO: other examples are defined in define_lattice.py script
energy_scale = 3.06

NTHREADS = int(sys.argv[1])
mod = sys.argv[2]
NMOMENTS = int(sys.argv[3])
L = int(sys.argv[4])

def graphene_initial(onsite=(0.0, 0.0)):
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

# number of decomposition parts in each direction of matrix. This divides the lattice into various sections,
# each of which is calculated in parallel
nx = ny = NTHREADS

# number of unit cells in each direction.


lx = L
ly = L

configuration = Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                      is_complex=False, precision=0, energy_scale=energy_scale)

# direction ID 'xx': 0, 'yy': 1, 'zz': 2, 'xy': 3, 'xz': 4, 'yx': 3, 'yz': 5, 'zx': 4, 'zy': 5
modification = Modification(magnetic_field=False)


calculation = Calculation(configuration)
calculation.conductivity_optical(num_points=1000, num_random=1, num_disorder=1, num_moments=NMOMENTS, direction='xx')
export_lattice(lattice, configuration, calculation, modification, mod+'optical_cond_xx.h5')

calculation = Calculation(configuration)
calculation.conductivity_optical(num_points=1000, num_random=1, num_disorder=1, num_moments=NMOMENTS, direction='xy')
export_lattice(lattice, configuration, calculation, modification, mod+'optical_cond_xy.h5')

calculation = Calculation(configuration)
calculation.conductivity_optical_nonlinear(num_points=1000, num_random=1, num_disorder=1, num_moments=NMOMENTS, 
        direction='xxx', temperature=1.0, special=1)
export_lattice(lattice, configuration, calculation, modification, mod+'nonlinear_optical_cond_xxx.h5')

calculation = Calculation(configuration)
calculation.dos(num_points=1000, num_random=1, num_disorder=1, num_moments=NMOMENTS)
export_lattice(lattice, configuration, calculation, modification, mod+'dos.h5')

calculation = Calculation(configuration)
calculation.singleshot_conductivity_dc(energy=[(n/100.0 - 0.5)*energy_scale*2 for n in range(101)], 
        num_moments=NMOMENTS, num_random=1, num_disorder=1, direction='xx', gamma=0.01)
export_lattice(lattice, configuration, calculation, modification, mod+'singleshot_xx.h5')

