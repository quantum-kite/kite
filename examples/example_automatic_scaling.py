import kite_config as kite
import matplotlib.pyplot as plt
import numpy as np
import pybinding as pb

from pybinding.repository import graphene

lattice = graphene.monolayer()

# add Disorder
disorder = kite.Disorder(lattice)
disorder.add_disorder('B', 'Deterministic', -0.3)
disorder.add_disorder('A', 'Uniform', +0.5, 0.5)

# add vacancy StructuralDisorder
disorder_struc = kite.StructuralDisorder(lattice, concentration=0.05)
disorder_struc.add_vacancy('A')

# node0 = [[+0, +0], 'A']
# node1 = [[+0, +0], 'B']
# disorder_struc = kite.StructuralDisorder(lattice, concentration=0.1)
# disorder_struc.add_structural_disorder(
#     # add onsite StructuralDisorder
#     (*node0, 2),
# )

# # if bond disorder is selected error is raised!
# node0 = [[+0, +0], 'A']
# node1 = [[+0, +0], 'B']
# node2 = [[+1, +0], 'A']
# node3 = [[+0, +1], 'B']
# node4 = [[+0, +1], 'A']
# node5 = [[-1, +1], 'B']
# disorder_struc = kite.StructuralDisorder(lattice, concentration=0.1)
# disorder_struc.add_structural_disorder(
#     # add bond disorder in the form [from unit cell], 'sublattice_from', [to_unit_cell], 'sublattice_to', value:
#     (*node0, *node3, 1),
#     (*node1, *node2, 1),
#     (*node2, *node3, 1),
#     (*node3, *node4, 1),
#     (*node4, *node5, 1),
#     (*node5, *node0, 1),
# )

nx = 2
ny = 2

lx = 1024
ly = 1024

# spectrum_range=[e_min, e_max] manually select Hamiltonian bounds.
# if spectrum_range is not selected automatic scaling of a Pybinding model with equivalent disorder strength
# is done in the background.

# example: automated scaling
# configuration = Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True], is_complex=False,
#                               precision=1)

# example: manual scaling
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=False, precision=1, spectrum_range=[-10, 10])

calculation = kite.Calculation(configuration)
calculation.dos(num_points=10000, num_moments=1024, num_random=1, num_disorder=1)

# make modification object which caries info about (TODO: Other modifications can be added here)
# - magnetic field can be set to True. Default case is False. In exported file it's converted to 1 and 0.
modification = kite.Modification(magnetic_field=False)

# export the lattice from the lattice object, config and calculation object and the name of the file
# the disorder is optional. If there is disorder in the lattice for now it should be given separately
kite.export_lattice(lattice, configuration, calculation, modification, 'auto_scaling.h5', disorder=disorder,
                    disorded_structural=disorder_struc)

# it's possible to make roughly an equivalent model using pybinding (more disorder realisations you select, more similar
# it will be)
# calculate DOS to confirm
model = kite.make_pybinding_model(lattice, disorder=disorder, disorder_structural=disorder_struc)
# if one wants to plot onsite map
# model.onsite_map.plot()
# plt.show()
kpm = pb.kpm(model)
dos = kpm.calc_dos(energy=np.linspace(-10, 10, 10000), broadening=1e-2, num_random=1)
dos.plot()
print(kpm.report())
plt.show()
