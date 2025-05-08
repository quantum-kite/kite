""" Density of states of graphene with On-site and vacancy disorder

    ##########################################################################
    #                         Copyright 2020/22, KITE                        #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in eV
    Lattice: Honeycomb
    Configuration: Periodic boundary conditions, double precision,
                    manual scaling, size of the system 256x256, without domain decomposition (nx=ny=1),
    Disorder: Disorder class Deterministic and Uniform at different sublattices,
               StructuralDisorder class vacancy and bond disorder
    Calculation type: Average DOS
    Note: automatic scaling is not supported when bond disorder is present
    Last updated: 08/05/2025
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
    disorder.add_disorder('B', 'Deterministic', -0.7)
    disorder.add_disorder('A', 'Uniform', 0.0, 0.5)

    # add vacancy StructuralDisorder
    disorder_structural = kite.StructuralDisorder(lattice, concentration=0.05)
    disorder_structural.add_vacancy('A')

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
        precision=1,
        spectrum_range=[-10, 10]
    )

    # specify calculation type
    calculation = kite.Calculation(configuration)
    calculation.dos(
        num_points=10000,
        num_moments=1024,
        num_random=1,
        num_disorder=1
    )

    # configure the *.h5 file
    output_file = "mixed_disorder-output.h5"
    kite.config_system(lattice, configuration, calculation, filename=output_file,
                       disorder=disorder, disorder_structural=disorder_structural)

    # for generating the desired output from the generated HDF5-file, run
    # ../build/KITEx mixed_disorder-output.h5
    # ../build/KITE-tools mixed_disorder-output.h5

    # returning the name of the created HDF5-file
    return output_file


if __name__ == "__main__":
    main()
