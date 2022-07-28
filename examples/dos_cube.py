""" Density of states of a cubic lattice

    ##########################################################################
    #                         Copyright 2022, KITE                           #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in eV
    Lattice: cubic
    Configuration: Periodic boundary conditions, double precision,
                    manual rescaling, size of the system 512x512, with domain decomposition (nx=ny=2)
    Calculation type: Average DOS
    Last updated: 28/07/2022
"""

__all__ = ["main"]

import kite
import numpy as np
import pybinding as pb


def cube(onsite=(0, 0)):
    """Return lattice specification for a cube lattice with nearest neighbor hoppings"""

    # parameters
    t = 1


    # define lattice vectors
    a1 = np.array([1, 0, 0])
    a2 = np.array([0, 1, 0])
    a3 = np.array([0, 0, 1])

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(a1=a1, a2=a2)

    # add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, 0, 0], onsite[0])
    )

    # Add hoppings
    lat.add_hoppings(
        # between neighboring cells, between which atoms, and the value
        ([1, 0, 0], 'A', 'A', -t),
        ([0, 1, 0], 'A', 'A', -t),
        ([0, 0, 1], 'A', 'A', -t)
    )
    return lat


def main(onsite=(0, 0)):
    """Prepare the input file for KITEx"""
    # load lattice
    lattice = cube(onsite)

    # number of decomposition parts [nx,ny] in each direction of matrix.
    # This divides the lattice into various sections, each of which is calculated in parallel
    nx = ny = nz = 1
    # number of unit cells in each direction.
    lx = ly = lz = 256

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
        divisions=[nx, ny, nz],
        length=[lx, ly, lz],
        boundaries=[mode, mode, mode],
        is_complex=False,
        precision=1,
        spectrum_range=[-6.1, 6.1]
    )

    # specify calculation type
    calculation = kite.Calculation(configuration)
    calculation.dos(
        num_points=1000,
        num_moments=256,
        num_random=1,
        num_disorder=1
    )

    # configure the *.h5 file
    output_file = "cube-output.h5"
    kite.config_system(lattice, configuration, calculation, filename=output_file)

    # for generating the desired output from the generated HDF5-file, run
    # ../build/KITEx cube-output.h5
    # ../tools/build/KITE-tools cube-output.h5

    # returning the name of the created HDF5-file
    return output_file


if __name__ == "__main__":
    main()
