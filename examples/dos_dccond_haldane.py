""" Density of States and DC conductivity of the Haldane model

    ##########################################################################
    #                         Copyright 2022, KITE                           #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in units of hopping, |t| = 1
    Lattice: Honeycomb
    Configuration: Periodic boundary conditions, double precision,
                    automatic scaling, size of the system 256x256, with domain decomposition (nx=ny=2)
    Disorder: Disorder class Uniform at different sublattices
    Calculation type: Average DOS and DC conductivity (xy)
    Last updated: 28/07/2022
"""

__all__ = ["main"]

import kite
import numpy as np
import pybinding as pb


def haldane(onsite=(0, 0), t=1):
    """Return lattice specification for Haldane model"""

    # parameters
    a = 0.24595  # [nm] unit cell length
    a_cc = 0.142  # [nm] carbon-carbon distance
    t2 = t/10

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
        ([0, -1], 'A', 'B', -t),
        ([1, 0], 'A', 'A', -t2 * 1j),
        ([0, -1], 'A', 'A', -t2 * 1j),
        ([-1, 1], 'A', 'A', -t2 * 1j),
        ([1, 0], 'B', 'B', -t2 * -1j),
        ([0, -1], 'B', 'B', -t2 * -1j),
        ([-1, 1], 'B', 'B', -t2 * -1j)
    )
    return lat


def main(onsite=(0, 0), t=1):
    """Prepare the input file for KITEx"""
    # load lattice
    lattice = haldane(onsite, t)

    # add Disorder
    disorder = kite.Disorder(lattice)
    disorder.add_disorder('A', 'Uniform', +0.0, 0.4)
    disorder.add_disorder('B', 'Uniform', +0.0, 0.4)

    # number of decomposition parts [nx,ny] in each direction of matrix.
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
        is_complex=True,
        precision=0,
        spectrum_range=[-10, 10]
    )

    # specify calculation type
    calculation = kite.Calculation(configuration)
    calculation.dos(
        num_points=1000,
        num_moments=256,
        num_random=1,
        num_disorder=1
    )

    # require the calculation conductivity_dc
    calculation.conductivity_dc(
        num_points=1000,
        num_moments=256,
        num_random=1,
        num_disorder=1,
        direction='xy',
        temperature=0.05
         )

    # configure the *.h5 file
    output_file = "haldane-output.h5"
    kite.config_system(lattice, configuration, calculation, filename=output_file,
                       disorder=disorder)

    # for generating the desired output from the generated HDF5-file, run
    # ../build/KITEx haldane-output.h5
    # ../build/KITE-tools haldane-output.h5

    # note: to generate  the conductivity data file for a desired window of Fermi energies, please use
    # ../build/KITE-tools h5_file.h --CondDC -F Emin Emax NumPoints
    # Run ../build/KITE-tools --help for more options

    # returning the name of the created HDF5-file
    return output_file


if __name__ == "__main__":
    main()
