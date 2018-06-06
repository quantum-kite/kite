""" Square lattice

    Lattice : Square 1[nm] interatomic distance and t=1[eV] hopping;
    Disorder : Disorder class Gaussian and Uniform at different sublattices,
               StructuralDisorder class vacancy and bond disorder;
    Configuration : size of the system 512x512, with domain decomposition (nx=ny=2), periodic boundary conditions,
                    double precision, automatic scaling;
    Calculation : dos;
    Modification : magnetic field is off;

"""

import kite
import numpy as np
import pybinding as pb


def square_lattice(onsite=(0, 0)):
    """Make a square lattice with nearest neighbor hopping"""

    a1 = np.array([1, 0])
    a2 = np.array([0, 1])

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(
        a1=a1, a2=a2
    )

    # Add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, 0], onsite[0])
    )

    # Add hoppings
    lat.add_hoppings(
        ([1, 0], 'A', 'A', - 1),
        ([0, 1], 'A', 'A', - 1)
    )

    return lat


# load a square lattice
lattice = square_lattice()
# number of decomposition parts in each direction of matrix. This divides the lattice into various sections,
# each of which is calculated in parallel
nx = ny = 2
# number of unit cells in each direction.
lx = ly = 512
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=False, precision=1)
# require the calculation of DOS
calculation = kite.Calculation(configuration)
calculation.dos(num_points=1000, num_moments=256, num_random=1, num_disorder=1)
# make modification object which caries info about
# - magnetic field can be set to True. Default case is False
modification = kite.Modification(magnetic_field=False)
kite.config_system(lattice, configuration, calculation, modification, 'square_lat.h5')
