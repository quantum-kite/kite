import kite

from pybinding.repository import graphene

lattice = graphene.monolayer()

nx = 1
ny = 1

lx = 128
ly = 128

# spectrum_range=[e_min, e_max] manually select Hamiltonian bounds.
# if spectrum_range is not selected automatic scaling of a Pybinding model with equivalent disorder strength
# is done in the background.
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=False, precision=1)

calculation = kite.Calculation(configuration)
calculation.dos(num_points=10000, num_moments=1024, num_random=1, num_disorder=1)

# magnetic field can be set either as:
#  - a value of magnetic_field, the closest commensurate value is returned
modification = kite.Modification(magnetic_field=400)
#  - an integer multiple of flux quantum
# modification = Modification(flux=1)

# export the lattice from the lattice object, config and calculation object and the name of the file
# the disorder is optional. If there is disorder in the lattice for now it should be given separately
kite.config_system(lattice, configuration, calculation, modification, 'magnetic_field.h5')
# adding the same magnetic field to the pybinding model won't work due to different size of the system!!!
