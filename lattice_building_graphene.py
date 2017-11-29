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

    # Add disorder
    # Each sublattice can have different disorder. If there are multiple orbitals at one sublattice, one needs to add
    # disorder vector of the same size as the number of orbitals. Type of disorder available are Gaussian,
    # Deterministic and Uniform. Each of the needs the have mean value, and standard deviation, where standard deviation
    # of deterministic disorder should be 0.

    disorder = ex.Disorder(lat)
    disorder.add_disorder('A', 'Deterministic', 0.0, 0)
    disorder.add_disorder('B', 'Deterministic', 0.0, 0)

    # Add bond disorder as an object of a class StructuralDisorder. In this manner we can add onsite and bond defects
    # with a specific concentration, which will be added to the simulated system. The procedure for adding is same
    # as adding the hopping, with the difference that the bond disorded is not bounded to one site in the [0, 0]
    # unit cell.
    
    # if there is disorder it should be returned separately from the lattice
    return lat, disorder, []


lattice, disorder, disorded_structural = graphene_initial()
# number of decomposition parts in each direction of matrix.

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
# - the name of the function
#   DOS - denstity of states == 1,
#   CondXX - conductivity in xx direction == 2,
#   CondXY - conductivity in xy direction == 3,
#   OptCond - optical conductivity == 4
#   SpinCond - spin conductivity == 5
#   SingleCondXX - single energy XX conductivity == 6
#   SingleCondXY - single energy XY conductivity == 7
# - number of moments for the calculation,
# - number of different random vector realisations,
# - number of disorder realisations.
# - energy and gamma for single energy calculations.
calculation = ex.Calculation(fname=['DOS', 'CondXX', 'SingleCondXX', 'SingleCondXY'],
                             num_moments=[256, 32, 256, 256], num_random=[1, 1, 1, 1],
                             num_disorder=[1, 1, 1, 1], energy=[0, 0], gamma=[1e-2, 1e-2])

# make modification object which caries info about (TODO: Other modifications can be added here)
# - magnetic field can be set to True. Default case is False. In exported file it's converted to 1 and 0.
modification = ex.Modification(magnetic_field=True)

# export the lattice from the lattice object, config and calculation object and the name of the file
# the disorder is optional. If there is disorder in the lattice for now it should be given separately
ex.export_lattice(lattice, configuration, calculation, modification, 'test_f.h5',
                  disorder=disorder, disorded_structural=disorded_structural)

# plotting the lattice
lattice.plot()
plt.show()
