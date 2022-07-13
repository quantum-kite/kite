""" Twisted bilayer graphene, loading a predefined lattice

    ##############################################################################
    #                        Copyright 2022, KITE                                #
    #                        Home page: quantum-kite.com                         #
    ##############################################################################

    Units: Energy in eV
    Lattice: Twisted bilayer graphene
    Disorder: None
    Configuration: Periodic boundary conditions, double precision,
                    automatic scaling, size of the system  flexible, with domain decomposition (nx=ny=1)
    Calculation type: Average DOS
    Last updated: 13/07/2022
"""

import pybinding as pb
import numpy as np
import kite


def twisted_bilayer_lattice(angle_index=0):
    # Return lattice specification for a twisted bilayer lattice with nearest neighbor hoppings

    # define the angle
    angle = np.array([2.005, 7.341, 13.174, 21.787])[angle_index]

    # define the name of the pb.Lattice object
    name = 'lattice_twisted_bilayer/lattice_tblg_{:.3f}'.format(angle)

    # load a predefined lattice object, the lattice can be saved with pb.save(lattice, name)
    lat = pb.load(name)
    return lat


if __name__ == "__main__":
    # select the twist angle, choose between
    #   0 ->  2.005
    #   1 ->  7.341
    #   2 -> 13.174
    #   3 -> 21.787
    angle_index = 0

    # load lattice
    lattice = twisted_bilayer_lattice(angle_index)

    # number of decomposition parts [nx,ny] in each direction of matrix.
    # This divides the lattice into various sections, each of which is calculated in parallel
    nx = ny = 1
    # number of unit cells in each direction.
    lx = ly = 64, 64

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
    calculation.dos(num_moments=1000,
                    num_random=1,
                    num_disorder=1,
                    num_points=1000)

    # configure the *.h5 file
    kite.config_system(lattice, configuration, calculation, filename='tblg_{:d}-data.h5'.format(angle_index))

    # for generating the desired output from the generated HDF5-file, run (replace n with the value os angle_index)
    # ../build/KITEx tbl_n.h5
    # ../tools/build/KITE-tools tbl_n-data.h5

    # Note: for a quick check, make a Pybinding model and check the DOS
    #   model = pb.Model(lattice, pb.rectangle(100, 100), pb.translational_symmetry(a1=50, a2=50))
    # To specify Disorder, use the function that converts the KITE disorder object to a Pybinding model
    #   model = kite.make_pybinding_model(lattice)
    #   dos = pb.kpm(model).calc_dos(np.linspace(-4, 4, 2000), broadening=1e-2, num_random=1)
    #   dos.plot()
