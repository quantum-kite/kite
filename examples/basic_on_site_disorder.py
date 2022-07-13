""" Density of states of graphene with On-site disorder

    ##############################################################################
    #                        Copyright 2022, KITE                                #
    #                        Home page: quantum-kite.com                         #
    ##############################################################################

    Units: Energy in eV
    Lattice: Honeycomb
    Configuration: Periodic boundary conditions, double precision,
                    manual scaling, size of the system 256x256, without domain decomposition (nx=ny=1),
    Disorder: Disorder class Deterministic and Uniform at different sublattices
    Calculation type: Average DOS
    Last updated: 13/07/2022
"""

import kite

from pybinding.repository import graphene


if __name__ == "__main__":
    # load a monolayer graphene lattice
    lattice = graphene.monolayer()

    # add Disorder
    disorder = kite.Disorder(lattice)
    disorder.add_disorder('B', 'Deterministic', -1.0)
    disorder.add_disorder('A', 'Uniform', +1.5, 1.0)

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
                                       precision=1)

    # specify calculation type
    calculation = kite.Calculation(configuration)
    calculation.dos(num_points=4000,
                    num_moments=512,
                    num_random=5,
                    num_disorder=1)

    # configure the *.h5 file
    kite.config_system(lattice, configuration, calculation, filename='on_site_disorder-data.h5',
                       disorder=disorder)

    # for generating the desired output from the generated HDF5-file, run
    # ../build/KITEx on_site_disorder-data.h5
    # ../tools/build/KITE-tools on_site_disorder-data.h5
