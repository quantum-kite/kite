""" Honeycomb lattice conductivity optical

    Lattice : Honeycomb 1[nm] interatomic distance and t=1[eV] hopping;
    Disorder : None;
    Configuration : size of the system 128x128, without domain decomposition (nx=ny=1), periodic boundary conditions,
                    double precision, automatic scaling;
    Calculation : conductivity_optical along the 'xx' direction;
    Modification : magnetic field is off;

"""

import kite
import numpy as np
import pybinding as pb


def honeycomb_lattice(onsite=(0, 0)):
    """Make a honeycomb lattice with nearest neighbor hopping"""

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
        ([-1, 1], 'A', 'B', - 1),
    )
    return lat


# load a honeycomb lattice
lattice = honeycomb_lattice()
# number of decomposition parts in each direction of matrix.
# This divides the lattice into various sections, each of which is calculated in parallel
nx = ny = 1
# number of unit cells in each direction.
lx = ly = 128
# make config object which caries info about
# - the number of decomposition parts [nx, ny],
# - lengths of structure [lx, ly]
# - boundary conditions, setting True as periodic boundary conditions, and False elsewise,
# - info if the exported hopping and onsite data should be complex,
# - info of the precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
# - scaling, if None it's automatic, if present select spectrum_bound=[e_min, e_max]
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=False, precision=1)
# require the calculation of optical conductivity
calculation = kite.Calculation(configuration)
calculation.conductivity_optical(num_points=1000, num_disorder=1, num_random=1, num_moments=128, direction='xx')
# make modification object which caries info about
# - magnetic field can be set to True. Default case is False
modification = kite.Modification(magnetic_field=False)
kite.config_system(lattice, configuration, calculation, modification, 'honeycomb_lat_cond_opt.h5')
