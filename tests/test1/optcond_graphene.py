import kite

from pybinding.repository import graphene
lattice = graphene.monolayer()
nx = ny = 2
lx = ly = 512
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=False, precision=1)
calculation = kite.Calculation(configuration)
# calculation.dos(num_points=4096, num_moments=512, num_random=10, num_disorder=1)
calculation.conductivity_optical(num_points=1024, num_disorder=1, num_random=1, num_moments=256, direction='xx')

kite.config_system(lattice, configuration, calculation, filename='graphene_optical.h5')
