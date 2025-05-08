""" Honeycomb lattice with vacancy disorder

    ##########################################################################
    #                         Copyright 2020/22, KITE                        #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in units of hopping, |t| = 1
    Lattice: Honeycomb
    Configuration: Periodic boundary conditions, double precision,
                    automatic scaling, size of the system 512x512, without domain decomposition (nx=ny=1)
    Disorder: StructuralDisorder, vacancy with concentration 0.1 inside A and 0.05 inside B sublattices
    Calculation type: Average DOS
    Last updated: 08/05/2025
"""

__all__ = ["main"]

import kite
import numpy as np
import pybinding as pb


def honeycomb_lattice(onsite=(0, 0), t=1):
    """Return lattice specification for a honeycomb lattice with nearest neighbor hoppings"""

    # parameters
    a = 1  # [nm] unit cell length
    a_cc = 1 / np.sqrt(3)  # [nm] carbon-carbon distance
    t = 2.8  # eV

    # define lattice vectors
    a1 = a * np.array([1, 0])
    a2 = a * np.array([1 / 2, 1 / 2 * np.sqrt(3)])

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(a1=a1, a2=a2)

    # add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, -a_cc/2], onsite[0]),
        ('B', [0,  a_cc/2], onsite[1])
    )

    # Add hoppings
    lat.add_hoppings(
        # inside the main cell, between which atoms, and the value
        ([0, 0], 'A', 'B', -t),
        # between neighboring cells, between which atoms, and the value
        ([1, -1], 'A', 'B', -t),
        ([0, -1], 'A', 'B', -t)
    )
    return lat


def main(onsite=(0, 0), t=1):
    """Prepare the input file for KITEx"""
    # load lattice
    lattice = honeycomb_lattice(onsite, t)

    # add vacancy StructuralDisorder
    # In this manner we can distribute vacancy disorder
    # on a specific sublattice with a specific concentration.
    # unless you would like the same pattern of disorder at both sublatices,
    # each realisation should be specified as a separate object
    struc_disorder_A = kite.StructuralDisorder(lattice, concentration=0.1)
    struc_disorder_A.add_vacancy('A')
    struc_disorder_B = kite.StructuralDisorder(lattice, concentration=0.05)
    struc_disorder_B.add_vacancy('B')
    disorder_structural = [struc_disorder_A, struc_disorder_B]

    # number of decomposition parts [nx,ny] in each direction of matrix.
    # This divides the lattice into various sections, each of which is calculated in parallel
    nx = ny = 1
    # number of unit cells in each direction.
    lx = ly = 512

    # make config object which caries info about
    # - the number of decomposition parts [nx, ny],
    # - lengths of structure [lx, ly]
    # - boundary conditions [mode,mode, ... ] with modes:
    #   . "periodic"
    #   . "open"
    #   . "twisted" -- this option needs the extra argument angles=[phi_1,..,phi_DIM] where phi_i \in [0, 2*M_PI]
    #   . "random"

    # Boundary Mode
    mode = "periodic"

    # - specify precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
    # - scaling, if None it's automatic, if present select spectrum_range=[e_min, e_max]
    configuration = kite.Configuration(
        divisions=[nx, ny],
        length=[lx, ly],
        boundaries=[mode, mode],
        is_complex=False,
        precision=1
    )

    # specify calculation type
    calculation = kite.Calculation(configuration)
    calculation.dos(
        num_points=1000,
        num_moments=512,
        num_random=1,
        num_disorder=1
    )

    # configure the *.h5 file
    output_file = "vacancies-output.h5"
    kite.config_system(lattice, configuration, calculation, filename=output_file, disorder_structural=disorder_structural)

    # for generating the desired output from the generated HDF5-file, run
    # ../build/KITEx vacancies-output.h5
    # ../build/KITE-tools vacancies-output.h5

    # returning the name of the created HDF5-file
    return output_file


if __name__ == "__main__":
    main()
