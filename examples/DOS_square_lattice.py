""" Density of states of a square lattice

    ##############################################################################
    #                        Copyright 2022, KITE                                #
    #                        Home page: quantum-kite.com                         #
    ##############################################################################

    Units: Energy in units of hopping, |t| = 1
    Lattice: Square lattice
    Configuration: Periodic boundary conditions, double precision, automatic rescaling
    Calculation type: Average DOS
    Last updated: 13/07/2022
"""

import kite
import numpy as np
import pybinding as pb

def square_lattice(onsite=(0, 0)):
    # Return lattice specification for a square lattice with nearest neighbor hoppings

    # parameters
    t = 1

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

if __name__ == "__main__":
    # load lattice
    lattice = square_lattice()

    # number of decomposition parts [nx,ny,nz] in each direction of matrix.
    # This divides the lattice into various sections, each of which is calculated in parallel
    nx = ny = 1
    # number of unit cells in each direction.
    lx = ly = 32

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
    configuration = kite.Configuration(divisions=[nx, ny],
                                       length=[lx, ly],
                                       boundaries=[mode, mode],
                                       is_complex=False,
                                       precision=1)

    # specify calculation type
    calculation = kite.Calculation(configuration)
    calculation.dos(num_points=4000,
                    num_moments=256,
                    num_random=256,
                    num_disorder=1)

    # configure the *.h5 file
    kite.config_system(lattice, configuration, calculation, filename='square_lattice-data.h5')

    # for generating the desired output from the generated HDF5-file, run
    # ../build/KITEx square_lattice-data.h5
    # ../tools/build/KITE-tools square_lattice-data.h5
