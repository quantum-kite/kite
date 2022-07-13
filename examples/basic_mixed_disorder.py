""" Density of states of graphene with On-site and vacancy disorder

    ##############################################################################
    #                        Copyright 2022, KITE                                #
    #                        Home page: quantum-kite.com                         #
    ##############################################################################

    Units: Energy in eV
    Lattice: Honeycomb
    Configuration: Periodic boundary conditions, double precision,
                    manual scaling, size of the system 256x256, without domain decomposition (nx=ny=1),
    Disorder: Disorder class Deterministic and Uniform at different sublattices,
               StructuralDisorder class vacancy and bond disorder
    Calculation type: Average DOS
    Note: automatic scaling is not supported when bond disorder is present
    Last updated: 13/07/2022
"""

__all__ = ["main"]

import kite

from pybinding.repository import graphene


def main(onsite=(0, 0)):
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
                                       precision=1,
                                       spectrum_range=[-10, 10])

    # specify calculation type
    calculation = kite.Calculation(configuration)
    calculation.dos(num_points=10000,
                    num_moments=2048,
                    num_random=1,
                    num_disorder=1)

    # configure the *.h5 file
    kite.config_system(lattice, configuration, calculation, filename='mixed_disorder-data.h5',
                       disorder=disorder, disorder_structural=disorder_structural)

    # for generating the desired output from the generated HDF5-file, run
    # ../build/KITEx mixed_disorder-data.h5
    # ../tools/build/KITE-tools mixed_disorder-data.h5


if __name__ == "__main__":
    main()
