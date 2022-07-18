""" Optical conductivity of disordered graphene lattice

    ##########################################################################
    #                         Copyright 2022, KITE                           #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in eV
    Configuration: Periodic boundary conditions, double precision,
                    given scaling, size of the system 512x512, with domain decomposition (nx=ny=2)
    Calculation type: Average DOS
    Last updated: 18/07/2022
"""

__all__ = ["main"]

import kite
import numpy as np
import pybinding as pb


def make_lattice(onsite=(0, 0), t=1):
    """Return lattice specification for a lattice with next nearest neighbor hoppings"""

    # define lattice vectors
    a1 = np.array([1.1, 0.1])
    a2 = np.array([0.3, 1.0])

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(a1=a1, a2=a2)

    # add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, 0], onsite[0]),
        ('B', [0, 0.1], onsite[1])
    )

    # Add hoppings
    lat.add_hoppings(
        # between neighboring cells, between which atoms, and the value
        ([1, 0], 'A', 'A', -t),
        ([0, 1], 'A', 'A', -t),
        ([1, 0], 'B', 'B', -t),
        ([0, 1], 'B', 'B', -t),
        ([0, 0], 'A', 'B', t / 10)
    )
    return lat


def main(onsite=(0, 0), t=1):
    """Prepare the input file for KITEx"""
    # load lattice
    lattice = make_lattice(onsite, t)

    # number of decomposition parts [nx,ny] in each direction of matrix.
    # This divides the lattice into various sections, each of which is calculated in parallel
    nx = ny = 2
    # number of unit cells in each direction.
    lx = ly = 512

    # make config object which caries info about
    # - the number of decomposition parts [nx, ny],
    # - lengths of structure [lx, ly]
    # - boundary conditions [mode,mode, ... ] with modes:
    #   . "periodic"
    #   . "open"
    #   . "twist_fixed" -- this option needs the extra argument ths=[phi_1,..,phi_DIM] where phi_i \in [0, 2*M_PI]
    #   . "twist_random"
    # Boundary Mode
    mode = "periodic"

    # - specify precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
    # - scaling, if None it's automatic, if present select spectrum_range=[e_min, e_max]
    e_min, e_max = -6, 6
    configuration = kite.Configuration(divisions=[nx, ny],
                                       length=[lx, ly],
                                       boundaries=[mode, mode],
                                       is_complex=False,
                                       precision=1,
                                       spectrum_range=[e_min, e_max],
                                       custom_local=True,
                                       custom_local_print=True)

    calculation = kite.Calculation(configuration)
    moments = 2048
    calculation.dos(num_points=1000,
                    num_moments=moments,
                    num_random=1,
                    num_disorder=1)

    # configure the *.h5 file
    output_file = "config.h5"
    kite.config_system(lattice, configuration, calculation, filename=output_file)

    # returning the name of the created HDF5-file
    return output_file


if __name__ == "__main__":
    main()
