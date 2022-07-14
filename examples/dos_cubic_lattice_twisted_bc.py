""" Density of states of a cubic lattice model (twisted boundary conditions)

    ##########################################################################
    #                         Copyright 2022, KITE                           #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in units of hopping, |t| = 1
    Lattice: Simple cubic lattice
    Configuration: twisted boundary conditions, double precision, automatic rescaling
    Calculation type: Average DOS
    Last updated: 14/07/2022
"""

__all__ = ["main"]

import kite
import numpy as np
import pybinding as pb


def cubic_lattice(onsite=(0, 0, 0), t=1):
    """Return lattice specification for a cubic lattice with nearest neighbor hoppings"""

    # define lattice vectors
    a1 = np.array([1, 0, 0])
    a2 = np.array([0, 1, 0])
    a3 = np.array([0, 0, 1])

    # create lattice with 3 primitive vectors
    lat = pb.Lattice(a1=a1, a2=a2, a3=a3)

    # Add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, 0, 0], onsite[0])
    )

    # Add hoppings
    lat.add_hoppings(
        # between neighboring cells, between which atoms, and the value
        ([1, 0, 1], 'A', 'A', -t),
        ([0, 1, 0], 'A', 'A', -t),
        ([0, 0, 1], 'A', 'A', -t)
    )
    return lat


def main(onsite=(0, 0, 0), t=1):
    """Prepare the input file for KITEx"""
    # load lattice
    lattice = cubic_lattice(onsite, t)

    # number of decomposition parts [nx,ny,nz] in each direction of matrix.
    # This divides the lattice into various sections, each of which is calculated in parallel
    nx = ny = nz = 1
    # number of unit cells in each direction.
    lx = ly = lz = 32

    # make config object which caries info about
    # - the number of decomposition parts [nx, ny, nz],
    # - lengths of structure [lx, ly, lz]
    # - boundary conditions [mode,mode, ... ] with modes:
    #   . "periodic"
    #   . "open"
    #   . "twist_fixed" -- this option needs the extra argument ths=[phi_1,..,phi_DIM] where phi_i \in [0, 2*M_PI]
    #   . "twist_random"
    # Boundary Mode
    mode = "random"

    # - specify precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
    # - scaling, if None it's automatic, if present select spectrum_range=[e_min, e_max]
    configuration = kite.Configuration(divisions=[nx, ny, nz],
                                       length=[lx, ly, lz],
                                       boundaries=[mode, mode, mode],
                                       is_complex=True)

    # specify calculation type
    calculation = kite.Calculation(configuration)
    calculation.dos(num_points=4000,
                    num_moments=256,
                    num_random=256,
                    num_disorder=1)

    # configure the *.h5 file
    output_file = "cubic_lattice-output.h5"
    kite.config_system(lattice, configuration, calculation, filename=output_file)

    # for generating the desired output from the generated HDF5-file, run
    # ../build/KITEx cubic_lattice-output.h5
    # ../tools/build/KITE-tools cubic_lattice-output.h5

    # returning the name of the created HDF5-file
    return output_file


if __name__ == "__main__":
    main()
