""" Density of states and DC conductivity of a square lattice

    ##########################################################################
    #                         Copyright 2020/2022, KITE                      #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in units of hopping, |t| = 1, Length in units of lattice parameter, |a| = 1
    Lattice: Square lattice
    Configuration: Periodic boundary conditions, double precision, manual rescaling,
                   size of the system 128x128, with domain decomposition (nx=ny=2)
    Calculation type: Average DOS, Conductivity DC
    Last updated: 28/07/2022
"""

__all__ = ["main"]

import kite
import numpy as np
import pybinding as pb


def square_lattice(onsite=[0], t=1):
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


def main(onsite=[0], t=1):
    """Prepare the input file for KITEx"""
    # load lattice
    lattice = square_lattice(onsite, t)

    # number of decomposition parts [nx,ny,nz] in each direction of matrix.
    # This divides the lattice into various sections, each of which is calculated in parallel
    nx = ny = 2
    # number of unit cells in each direction.
    lx = ly = 128

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
        precision=1,
        spectrum_range=[-4.1, 4.1]
    )

    # specify calculation type
    calculation = kite.Calculation(configuration)
    calculation.dos(
        num_points=4096,
        num_moments=128,
        num_random=10,
        num_disorder=1
    )
    calculation.conductivity_dc(
        num_points=1000,
        num_moments=256,
        num_random=10,
        direction='xy',
        temperature=0.001
    )

    mod = kite.Modification(magnetic_field=40)

    # configure the *.h5 file
    output_file = "square_lattice_condDC-output.h5"
    kite.config_system(lattice, configuration, calculation, modification=mod, filename=output_file)

    # for generating the desired output from the generated HDF5-file, run
    # ../build/KITEx square_lattice_condDC-output.h5
    # ../build/KITE-tools square_lattice_condDC-output.h5

    # returning the name of the created HDF5-file
    return output_file


if __name__ == "__main__":
    main()
