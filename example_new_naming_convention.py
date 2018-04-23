import matplotlib.pyplot as plt
import export_lattice as ex
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
        ([0, 0], 'A', 'B', - 1 / energy_scale),
        # between neighboring cells, between which atoms, and the value
        ([-1, 0], 'A', 'B', - 1 / energy_scale),
        ([-1, 1], 'A', 'B', - 1 / energy_scale)
    )

    return lat


lattice = graphene_initial()

# number of decomposition parts in each direction of matrix. This divides the lattice into various sections, each of which is calculated in parallel
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
configuration = ex.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                 is_complex=False, precision=1, energy_scale=energy_scale)


# make calculation object which caries info about
# - the name of the function, case insensitive choice

#   dos - denstity of states == 1,
#   conductivity_optical - conductivity in xx direction == 2,
#   conductivity_dc - conductivity in xy direction == 11,
#   conductivity_optical_nonlinear - optical conductivity == 20
#   singleshot_conductivity_dc - singleshot conductivity == 29

# - number of moments for the calculation,
# - number of different random vector realisations,
# - number of disorder realisations,
# - number of points for evaluating full spectrum calculation,
# - temperature at which we want to calculate full spectrum conductivity,
# - energy and gamma for single energy calculations.

# adding a direction number:
# 'xx': 0,
# 'xy': 1,
# 'xz': 2,
# 'yx': 3,
# 'yy': 4,
# 'yz': 5,
# 'zx': 6
# 'zy': 7
# 'zz': 8

# number that will distinguish the quantity is fun_number[type] + avail_dir/avail_dir_spec[direction]
# for example conductivity_optical_nonlinear in yx direction is 20 + 3 = 23

calculation = ex.Calculation(fname=['DOS', 'conductivity_dc', 'conductivity_optical_nonlinear',
                                    'singleshot_conductivity_dc'],
                             num_moments=[1024, 1, 1, 1], num_random=[1, 1, 1, 1], temperature=[100, 0],
                             num_disorder=[1, 1, 1, 1], direction=['xx', 'xy', 'zz'], energy=[0,0], gamma=[0.1], special=1)

# make modification object which caries info about (TODO: Other modifications can be added here)
# - magnetic field can be set to True. Default case is False. In exported file it's converted to 1 and 0.
modification = ex.Modification(magnetic_field=False)

# export the lattice from the lattice object, config and calculation object and the name of the file
# the disorder is optional. If there is disorder in the lattice for now it should be given separately
ex.export_lattice(lattice, configuration, calculation, modification, 'example_new_naming_convention.h5')

# plotting the lattice
lattice.plot()
# plt.show()
