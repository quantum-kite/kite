""" Optical conductivity of disordered graphene lattice

    ##########################################################################
    #                         Copyright 2022, KITE                           #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in eV
    Lattice: Honeycomb
    Configuration: Periodic boundary conditions, double precision,
                    automatic scaling, size of the system 512x512, with domain decomposition (nx=ny=1)
    Disorder: Disorder class Gaussian at different sublattices
    Calculation type: Average DOS and optical conductivity (xx)
    Last updated: 18/07/2022
"""

__all__ = ["main"]

import kite

from pybinding.repository import graphene


def main(onsite=(0, 0)):
    """Prepare the input file for KITEx"""
    # load a monolayer graphene lattice
    lattice = graphene.monolayer(onsite=onsite)

    # add Disorder
    disorder = kite.Disorder(lattice)
    disorder.add_disorder('B', 'Gaussian', 0.0, 0.6)
    disorder.add_disorder('A', 'Gaussian', 0.0, 0.6)

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
    #   . "twisted" -- this option needs the extra argument ths=[phi_1,..,phi_DIM] where phi_i \in [0, 2*M_PI]
    #   . "random"

    # Boundary Mode
    mode = "periodic"

    # - specify precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
    # - scaling, if None it's automatic, if present select spectrum_range=[e_min, e_max]
    configuration = kite.Configuration(divisions=[nx, ny],
                                       length=[lx, ly],
                                       boundaries=[mode, mode],
                                       is_complex=False,
                                       precision=1,
                                       spectrum_range=[-10, 10])

    # specify calculation type
    calculation = kite.Calculation(configuration)
    calculation.dos(num_points=4096,
                    num_moments=512,
                    num_random=10,
                    num_disorder=1)
    # require the calculation conductivity_optical (xx)
    calculation.conductivity_optical(num_points=1024,
                                     num_disorder=1,
                                     num_random=5,
                                     num_moments=256,
                                     direction='xx')

    # configure the *.h5 file
    output_file = "optcond_gaussian_disorder-output.h5"
    kite.config_system(lattice, configuration, calculation, filename=output_file,
                       disorder=disorder)

    # for generating the desired output from the generated HDF5-file, run
    # ../build/KITEx optcond_gaussian_disorder-output.h5
    # ../tools/build/KITE-tools optcond_gaussian_disorder-output.h5

    # returning the name of the created HDF5-file
    return output_file


if __name__ == "__main__":
    main()
