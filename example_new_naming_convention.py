import matplotlib.pyplot as plt
import export_lattice as ex
import numpy as np
import pybinding as pb

from export_lattice import Configuration, Calculation, Modification, export_lattice

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

# number of decomposition parts in each direction of matrix. This divides the lattice into various sections,
# each of which is calculated in parallel
nx = ny = 1

# number of unit cells in each direction.
lx = 256
ly = 256

# make config object which caries info about
# - the number of decomposition parts [nx, ny],
# - lengths of structure [lx, ly]
# - boundary conditions, setting True as periodic boundary conditions, and False elsewise,
# - info if the exported hopping and onsite data should be complex,
# - info of the precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
configuration = Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                              is_complex=False, precision=0, energy_scale=energy_scale)

# make calculation object which caries info about
# - the desired function:

#   dos - denstity of states,
#   conductivity_optical - optical conductivity in linear response,
#   conductivity_dc - zero frequency conductivity in linear response,
#   conductivity_optical_nonlinear - optical conductivity nonlinear response,
#   singleshot_conductivity_dc - zero frequency singleshot conductivity.

# direction ID 'xx': 0, 'yy': 1, 'zz': 2, 'xy': 3, 'xz': 4, 'yx': 3, 'yz': 5, 'zx': 4, 'zy': 5
calculation = Calculation(configuration)
calculation.dos(num_points=1000, num_random=1, num_moments=1024)
calculation.conductivity_optical(num_points=1000, num_random=1, num_moments=256, direction='xx')
calculation.conductivity_dc(num_points=1000, num_random=2, num_moments=256, direction='xx', temperature=300)
calculation.conductivity_dc(num_points=1000, num_random=1, num_moments=256, direction='xy', temperature=300)
calculation.singleshot_conductivity_dc(energy=[0, 0.1], num_random=1, num_moments=256, direction='xx', gamma=0.1)
calculation.singleshot_conductivity_dc(energy=[0, 0.1], num_random=1, num_moments=256, direction='yy', gamma=0.1)
calculation.conductivity_optical_nonlinear(num_points=1000, num_random=1, num_moments=256, direction='xx',
                                           temperature=0, special=1)

# make modification object which caries info about (TODO: Other modifications can be added here)
# - magnetic field can be set to True. Default case is False. In exported file it's converted to 1 and 0.
modification = Modification(magnetic_field=False)

# export the lattice from the lattice object, config and calculation object and the name of the file
# the disorder is optional. If there is disorder in the lattice for now it should be given separately
export_lattice(lattice, configuration, calculation, modification, 'example_new_naming_convention.h5')

# plotting the lattice
lattice.plot()
# plt.show()
