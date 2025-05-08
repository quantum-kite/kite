""" Density of states of a square lattice (twisted boundary conditions)

    ##########################################################################
    #                         Copyright 2020/22, KITE                        #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in units of hopping, |t| = 1
    Lattice: Square lattice
    Configuration: Twisted boundary conditions, double precision, automatic rescaling
    Calculation type: Average DOS
    Last updated: 08/05/2025
"""

__all__ = ["main"]

import kite
import numpy as np
import pybinding as pb


def square_lattice(onsite=(0, 0), t=1):
    """Return lattice specification for a square lattice with nearest neighbor hoppings"""

    # define lattice vectors
    a1 = np.array([1, 0])
    a2 = np.array([0, 1])

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(a1=a1, a2=a2)

    # add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, 0], onsite[0])
    )

    # Add hoppings
    lat.add_hoppings(
        # between neighboring cells, between which atoms, and the value
        ([1, 0], 'A', 'A', -t),
        ([0, 1], 'A', 'A', -t)
    )
    return lat


def main(onsite=(0, 0), t=1):
    """Prepare the input file for KITEx"""
    # load lattice
    lattice = square_lattice(onsite, t)

    # number of decomposition parts [nx,ny] in each direction of matrix.
    # This divides the lattice into various sections, each of which is calculated in parallel
    nx = ny = 2
    # number of unit cells in each direction.
    lx = ly = 64

    # make config object which caries info about
    # - the number of decomposition parts [nx, ny],
    # - lengths of structure [lx, ly]
    # - boundary conditions [mode,mode, ... ] with modes:
    #   . "periodic"
    #   . "open"
    #   . "twisted" -- this option needs the extra argument angles=[phi_1,..,phi_DIM] where phi_i \in [0, 2*M_PI]
    #   . "random"

    # Boundary Mode
    mode = "twisted"

    # Twists in each direction
    twsx = twsy = np.pi/2.0

    # - specify precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
    # - scaling, if None it's automatic, if present select spectrum_range=[e_min, e_max]
    configuration = kite.Configuration(
        divisions=[nx, ny],
        length=[lx, ly],
        boundaries=[mode, mode],
        is_complex=True,
        precision=1,
        angles=[twsx, twsy]
    )

    # specify calculation type
    calculation = kite.Calculation(configuration)
    calculation.dos(
        num_points=4000,
        num_moments=512,
        num_random=32,
        num_disorder=1
    )

    # configure the *.h5 file
    output_file = "square_lattice_twisted_bc-output.h5"
    kite.config_system(lattice, configuration, calculation, filename=output_file)

    # for generating the desired output from the generated HDF5-file, run
    # ../build/KITEx square_lattice_twisted_bc-output.h5
    # ../build/KITE-tools square_lattice_twisted_bc-output.h5

    # returning the name of the created HDF5-file
    return output_file


if __name__ == "__main__":
    main()
