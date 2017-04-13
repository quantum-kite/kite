import matplotlib.pyplot as plt
import define_lattice as dl
import export_lattice as ex


# define lattice of monolayer graphene by setting the number of nearest neighbors and onsite energies at atoms.
# - onsite energies can be a number, array, or matrix, where for the last case off-diagonal terms represent orbital
# hoppings at the same lattice site.
lattice = dl.monolayer_graphene(nearest_neighbors=1, onsite=([[1, 1], [1, 2]], [[1, 1], [1, 2]]))

# plotting the lattice
lattice.plot()
plt.show()

# number of decomposition parts in each direction of matrix.
nx = ny = 1
# number of unit cells in each direction.
lx = ly = 1
# make config object which caries info about
# - the number of decomposition parts [nx, ny],
# - lengths of structure [lx, ly]
# - boundary conditions, setting True as periodic boundary conditions, and False elsewise,
# - info if the exported hopping and onsite data should be complex,
# - info of the precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
config = ex.Config(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True], is_complex=False, precision=1)

# make calculation object which caries info about
# - the name of the function
#   DOS - denstity of states == 1,
#   CondXX - conductivity in xx direction == 2,
#   CondXY - conductivity in xy direction == 3,
#   OptCond - optical conductivity == 4
#   SpinCond - spin conductivity == 5
# - number of moments for the calculation  == 6,
# - number of different random vector realisations  == 7,
# - number of disorder realisations == 8.
calculation = ex.Calculation(fname='DOS', num_moments=1024, num_random=50, num_disorder=1)

# make modification object which caries info about
# - disorder realisation, for now you can chose between 'rectangular' and 'gaussian' distribution and choose its mean
# value and width
# - magnetic field can be set to True. Default case is False. In exported file it's converted to 1 and 0.
modification = ex.Modifications(disorder={'name': 'rectangular', 'width': 2, 'mean_value': 1})

# export the lattice from the lattice object, config and calculation object and the name of the file
ex.export_lattice(lattice, config, calculation, modification, 'test.h5')
