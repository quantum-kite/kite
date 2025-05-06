""" Density of states for twisted bilayer graphene

    ##########################################################################
    #                         Copyright 2022, KITE                           #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in eV
    Lattice: Twisted bilayer graphene
    Disorder: None
    Configuration: Periodic boundary conditions, double precision,
                   automatic scaling, size of the system  flexible, with domain decomposition (nx=ny=2)
    Calculation type: Average DOS
    Last updated: 28/07/2022
"""

__all__ = ["main"]

import pybinding as pb
import numpy as np
import kite


def twisted_bilayer_lattice(angle_index=0):
    """Return lattice specification for a twisted bilayer lattice with nearest neighbor hoppings"""

    # define the angle
    angle = np.array([2.005, 7.341, 13.174, 21.787])[angle_index]

    # define the name of the pb.Lattice object
    name = '.lattice_twisted_bilayer/lattice_tblg_{:.3f}'.format(angle)

    # load a predefined lattice object, the lattice can be saved with pb.save(lattice, name)
    lat = pb.load(name)
    return lat


def main(angle_index=3):
    """Prepare the input file for KITEx"""
    # select the twist angle, choose between
    #   angle_index     twist angle (degrees)
    #     0         ->       2.005
    #     1         ->       7.341
    #     2         ->      13.174
    #     3         ->      21.787

    # load lattice
    lattice = twisted_bilayer_lattice(angle_index)

    # number of decomposition parts [nx,ny] in each direction of matrix.
    # This divides the lattice into various sections, each of which is calculated in parallel
    nx = ny = 2
    # number of unit cells in each direction.
    lx = ly = 128,128

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
        num_moments=2000,
        num_random=1,
        num_disorder=1,
        num_points=1000
    )

    # configure the *.h5 file
    output_file = "tblg_{0}-output.h5".format(str(angle_index))
    kite.config_system(lattice, configuration, calculation, filename=output_file)

    # for generating the desired output from the generated HDF5-file, run (replace 0 with the value of angle_index)
    # ../build/KITEx tblg_0-output.h5
    # ../build/KITE-tools tblg_0-output.h5

    # Note: for a quick check, make a Pybinding model and check the DOS
    #   model = pb.Model(lattice, pb.rectangle(100, 100), pb.translational_symmetry(a1=50, a2=50))
    # To specify Disorder, use the function that converts the KITE disorder object to a Pybinding model
    #   model = kite.make_pybinding_model(lattice)
    #   dos = pb.kpm(model).calc_dos(np.linspace(-4, 4, 2000), broadening=1e-2, num_random=1)
    #   dos.plot()

    # returning the name of the created HDF5-file
    return output_file


if __name__ == "__main__":
    main()
