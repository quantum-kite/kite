""" Magnetic field

    Lattice : Monolayer graphene;
    Disorder : None;
    Configuration : size of the system 256x256, without domain decomposition (nx=ny=1), periodic boundary conditions,
                    double precision, manual scaling;
    Calculation : dos;
    Modification : magnetic field is On;

    Note : automatic scaling is not taking into account the magnetic field, an assumption is that the spectrum bounds
    will not change!
"""

import kite

from pybinding.repository import graphene

# load a monolayer graphene lattice
lattice = graphene.monolayer()
# number of decomposition parts in each direction of matrix. This divides the lattice into various sections,
# each of which is calculated in parallel
nx = ny = 1
# number of unit cells in each direction.
lx = ly = 256
# make config object which caries info about
# - the number of decomposition parts [nx, ny],
# - lengths of structure [lx, ly]
# - boundary conditions, setting True as periodic boundary conditions, and False elsewise,
# - info if the exported hopping and onsite data should be complex,
# - info of the precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
# - scaling, if None it's automatic, if present select spectrum_bound=[e_min, e_max]
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=False, precision=1)
# require the calculation of DOS
calculation = kite.Calculation(configuration)
calculation.dos(num_points=1000, num_moments=1024, num_random=1, num_disorder=1)
# magnetic field can be set either as:
#  - a value of magnetic_field, the closest commensurate value is returned
modification = kite.Modification(magnetic_field=400)
#  - an integer multiple of flux quantum
# modification = Modification(flux=1)

# configure the *.h5 file
kite.config_system(lattice, configuration, calculation, modification, 'magnetic_field.h5')
# adding the same magnetic field to the pybinding model won't work due to different size of the system!!!
