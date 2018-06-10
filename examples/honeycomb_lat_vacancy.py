""" Honeycomb lattice with vacancy disorder

    Lattice : Honeycomb 1[nm] interatomic distance and t=1[eV] hopping;
    Disorder : StructuralDisorder, vacancy with concentration 0.004 inside both A and B sublattices;
    Configuration : size of the system 256x256, without domain decomposition (nx=ny=1), periodic boundary conditions,
                    double precision, automatic scaling;
    Calculation : dos and singleshot conductivity at different energies;
    Modification : magnetic field is off;

"""

import kite
import numpy as np
import pybinding as pb


def honeycomb_lattice(onsite=(0, 0)):
    """Make a honeycomb lattice with nearest neighbor hopping
    
    Parameters
    ----------
    onsite : tuple or list
        Onsite energy at different sublattices.
    """""

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

    # Add vacancy disorder as an object of a class StructuralDisorder. In this manner we can distribute vacancy disorder
    # on a specific sublattice with a specific concentration.
    struc_disorder = kite.StructuralDisorder(lat, concentration=0.004)
    struc_disorder.add_vacancy('A')
    struc_disorder.add_vacancy('B')

    return lat, struc_disorder


# load a honeycomb lattice and structural disorder
lattice, disorder_structural = honeycomb_lattice()

nx = ny = 1
lx = 256
ly = 256
# make config object which caries info about
# - the number of decomposition parts [nx, ny],
# - lengths of structure [lx, ly]
# - boundary conditions, setting True as periodic boundary conditions, and False elsewise,
# - info if the exported hopping and onsite data should be complex,
# - info of the precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
# - scaling, if None it's automatic, if present select spectrum_bound=[e_min, e_max]
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=False, precision=1)
# require the calculation of dos and single shot conductivity
num_points = 10
num_moments = 512
eta = 0.1
calculation = kite.Calculation(configuration)
calculation.dos(num_points=1000, num_moments=num_moments, num_random=1, num_disorder=1)
calculation.singleshot_conductivity_dc(direction='xx', num_moments=num_moments, num_random=1, num_disorder=1,
                                       energy=[0.3 / num_points * i for i in range(num_points)], eta=eta)
# configure the *.h5 file
kite.config_system(lattice, configuration, calculation, filename='honeycomb_lat_vacancy.h5',
                   disorder_structural=disorder_structural)
