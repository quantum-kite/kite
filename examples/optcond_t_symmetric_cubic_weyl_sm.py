""" Optical Conductivity of a T-Symmetric Cubic Weyl Semimetal with Anderson Disorder

    ##########################################################################
    #                         Copyright 2020/22, KITE                        #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in units of nearest neighbour hopping, |t| = 1
    Lattice: simple cubic lattice (two orbitals per site)
    Configuration: random twisted boundary conditions, double precision, manual rescaling
    Calculation type: Linear Optical Conductivity (sigma_xx) [Heavily Paralellized Example]
    Last updated: 08/05/2025
"""

__all__ = ["main"]

import kite
import numpy as np
import pybinding as pb


def weyl_semimetal(t=1):
    """Return lattice specification for a Weyl Semimetal lattice"""

    # define lattice vectors
    a1 = np.array([1, 0, 0])
    a2 = np.array([0, 1, 0])
    a3 = np.array([0, 0, 1])

    # create a lattice with 3 primitive vectors
    lat = pb.Lattice(a1=a1, a2=a2, a3=a3)

    # add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, 0, 0], 0.0),
        ('B', [0, 0, 0], 0.0)
    )

    # Add hoppings
    lat.add_hoppings(
        # between neighboring cells, between which atoms, and the value
        ([1, 0, 0], 'A', 'B', 0.5 * t),
        ([1, 0, 0], 'B', 'A', 0.5 * t),
        ([0, 1, 0], 'A', 'B', 0.5 * t * 1j),
        ([0, 1, 0], 'B', 'A', -0.5 * t * 1j),
        ([0, 0, 1], 'A', 'A', 0.5 * t),
        ([0, 0, 1], 'B', 'B', -0.5 * t)
    )
    return lat


def main(anderson_w=0):
    """Prepare the input file for KITEx"""
    # load lattice
    lattice = weyl_semimetal()

    # add scalar on-site Disorder (Box distribution in [W/2,W/2])
    # Number of Samples controlled by num_disorder
    disorder = kite.Disorder(lattice)
    disorder.add_disorder('A', 'Uniform', 0.0, anderson_w/np.sqrt(12))
    disorder.add_disorder('B', 'Uniform', 0.0, anderson_w/np.sqrt(12))

    # number of decomposition parts [nx,ny,nz] in each direction of matrix.
    # This divides the lattice into various sections, each of which is calculated in parallel
    nx = 2
    ny = nz = 4
    # number of unit cells in each direction.
    lx = ly = lz = 128

    # make config object which caries info about
    # - the number of decomposition parts [nx, ny,nz],
    # - lengths of structure [lx, ly,lz]
    # - boundary conditions [mode,mode, ... ] with modes:
    #   . "periodic"
    #   . "open"
    #   . "twisted" -- this option needs the extra argument angles=[phi_1,..,phi_DIM] where phi_i \in [0, 2*M_PI]
    #   . "random" -- The number of random twists is fixed by num_disorder

    # Boundary Mode
    mode = "random"

    # - specify precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
    # - scaling, if None it's automatic, if present select spectrum_range=[e_min, e_max]
    configuration = kite.Configuration(
        divisions=[nx, ny, nz],
        length=[lx, ly, lz],
        boundaries=[mode, mode, mode],
        is_complex=True,
        precision=1,
        spectrum_range=[-1.8 - anderson_w, 1.8 + anderson_w]
    )

    # specify calculation type
    calculation = kite.Calculation(configuration)
    calculation.conductivity_optical(
        num_points=512,
        num_moments=1024,
        num_random=1,
        num_disorder=1,
        direction="xx",
        temperature=0.01
    )

    # configure the *.h5 file
    output_file = "weyl_optical-output.h5"
    kite.config_system(lattice, configuration, calculation, filename=output_file, disorder=disorder)

    # for generating the desired output from the generated HDF5-file, run
    # ../build/KITEx weyl_optical-output.h5
    # ../build/KITE-tools weyl_optical-output.h5 --CondOpt -O 0.02 1 401 -E 1000 -F 0.25 -S 0.005

    # returning the name of the created HDF5-file
    return output_file


if __name__ == "__main__":
    main()
