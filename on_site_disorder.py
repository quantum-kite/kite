import matplotlib.pyplot as plt
import numpy as np
import pybinding as pb
from pybinding.repository import graphene

from export_lattice import Configuration, Calculation, Modification, Disorder, StructuralDisorder, \
    export_lattice, make_pybinding_model


lattice = graphene.monolayer()

# add Disorder
disorder = Disorder(lattice)
disorder.add_disorder('B', 'Deterministic', -1.0)
disorder.add_disorder('A', 'Uniform', +1.5, 1.0)



nx = 1
ny = 1

lx = 512
ly = 512

# spectrum_range=[e_min, e_max] manually select Hamiltonian bounds.
# if spectrum_range is not selected automatic scaling of a Pybinding model with equivalent disorder strength
# is done in the background.
configuration = Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                              is_complex=False, precision=1)

calculation = Calculation(configuration)
calculation.dos(num_points=5000, num_moments=512, num_random=1, num_disorder=1)

# make modification object which caries info about (TODO: Other modifications can be added here)
# - magnetic field can be set to True. Default case is False. In exported file it's converted to 1 and 0.
modification = Modification(magnetic_field=False)

# export the lattice from the lattice object, config and calculation object and the name of the file
# the disorder is optional. If there is disorder in the lattice for now it should be given separately
export_lattice(lattice, configuration, calculation, modification, 'on_site_disorder.h5', disorder=disorder)
