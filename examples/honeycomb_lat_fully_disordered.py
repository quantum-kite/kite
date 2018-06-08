""" Honeycomb lattice where one bond is added as structural disorder

    Lattice : Honeycomb 1[nm] interatomic distance and t=1[eV] hopping;
    Disorder : StructuralDisorder adding a bond inside the full lattice as a disorder with concentration 1;
    Configuration : size of the system 256x256, without domain decomposition (nx=ny=1), periodic boundary conditions,
                    double precision, manual scaling;
    Calculation : dos;
    Modification : magnetic field is off;

"""

import kite
import numpy as np
import pybinding as pb


def honeycomb_lattice(onsite=(0, 0)):
    """Make a honeycomb lattice with nearest neighbor hopping"""

    # define lattice vectors
    theta = np.pi / 3
    a1 = np.array([1 + np.cos(theta), np.sin(theta)])
    a2 = np.array([0, 2 * np.sin(theta)])

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(
        a1=a1, a2=a2
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
        ([+0, +0], 'A', 'B', - 1),
        # between neighboring cells, between which atoms, and the value
        ([-1, +0], 'A', 'B', - 1),
        ([-1, +1], 'A', 'B', - 1)
    )

    # Add bond disorder as an object of a class StructuralDisorder. In this manner we can add onsite and bond defects
    # with a specific concentration, which will be added to the simulated system. The procedure for adding is same
    # as adding the hopping, with the only difference that the bond disorded is not bounded to one site in the [0, 0]
    # unit cell.
    node0 = [[+0, +0], 'A']
    node1 = [[+0, +0], 'B']

    struc_disorder = kite.StructuralDisorder(lat, concentration=1.00)
    struc_disorder.add_structural_disorder(
        (*node0, *node1, -1)
    )
    return lat, struc_disorder


# load a honeycomb lattice and structural disorder
lattice, disorder_structural = honeycomb_lattice()
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
e_min, e_max = -3.06, 3.06
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=False, precision=1, spectrum_range=[e_min, e_max])
# require the calculation of DOS
calculation = kite.Calculation(configuration)
calculation.dos(num_points=1000, num_random=1, num_disorder=1, num_moments=256)
# configure the *.h5 file
kite.config_system(lattice, configuration, calculation, filename='honeycomb_lat_fully_disordered.h5',
                   disorder_structural=disorder_structural)
