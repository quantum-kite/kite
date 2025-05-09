""" Density of states for twisted bilayer graphene

    ##########################################################################
    #                         Copyright 2020/22, KITE                        #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in eV
    Lattice: Twisted bilayer graphene
    Disorder: None
    Configuration: Periodic boundary conditions, double precision,
                    given rescaling, size of the system flexible, with domain decomposition (nx=ny=1)
    Calculation type: Average DOS
    Last updated: 08/05/2025
"""

__all__ = ["main"]

import pybinding as pb
import numpy as np
import os
import re
import copy
import kite


def _load_ovito_lattice(name):
    """Load a lattice from a Ovito .xyz file. At the moment, works for files that have following format
       Type_id X_coord Y_coord Z_coord Data
       The point needs to have all 3 coordinates, and the approach is not general.

    Parameters
    ----------

    name : str
        Name of the xyz file

    Return
    ----------

        x, y, z coordinates, atom types and l1, l2, l3 unit cell vectors
    """
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

    load_f = open(os.path.join(__location__, name), 'r')

    # read the atom number
    num_atoms = load_f.readline()

    # read the lattice vectors
    vectors = load_f.readline()
    vec = re.findall(r"[-+]?\d*\.*\d+", vectors)
    vec = list(map(float, vec))

    space_size = int(np.sqrt(len(vec)))

    _l1, _l2, _l3 = vec[0:space_size], vec[space_size:2 * space_size], vec[2 * space_size:3 * space_size]

    _atom_type = []
    _x_coord = []
    _y_coord = []
    _z_coord = []
    for line in load_f:
        atom = []
        for u in line.split():
            u = u.strip('\'')
            atom.append(u)
        _atom_type.append(atom[0])
        _x_coord.append(float(atom[1]) * 0.1)
        _y_coord.append(float(atom[2]) * 0.1)
        _z_coord.append(float(atom[3]) * 0.1)

    _x_coord = np.asarray(_x_coord)
    _y_coord = np.asarray(_y_coord)
    _z_coord = np.asarray(_z_coord)
    _atom_type = np.asarray(_atom_type)

    return _x_coord, _y_coord, _z_coord, _atom_type, _l1, _l2, _l3


def twisted_bilayer_lattice(name='relaxed_tblg_1.050.xyz'):
    # WARNING: the calculation ran with parameters given here will run for at least 72h
    # load an xyz file, relaxed or unrelaxed
    lat = pb.load('lattice_' + name[:-4])
    return lat


def twisted_bilayer_period(name='relaxed_tblg_1.050.xyz'):
    # WARNING: the calculation ran with parameters given here will run for at least 72h
    # load an xyz file, relaxed or unrelaxed
    # load coordinates, atom types and lattice vectors
    _, _, _, _, l1, l2, _ = _load_ovito_lattice(name)
    l1 = np.array(l1[0:2]) / 10.
    l2 = np.array(l2[0:2]) / 10.
    return l1, l2


def main():
    """Prepare the input file for KITEx"""
    # load lattice
    name='relaxed_tblg_1.050.xyz'
    lattice = twisted_bilayer_lattice(name)
    l1, l2 = twisted_bilayer_period(name)
    print('Done loading ', flush=True)

    # make the pybidning model
    print('Making the model', flush=True)
    # check for the bonds
    model = pb.Model(lattice,
                     pb.force_double_precision(),
                     pb.translational_symmetry(a1=np.linalg.norm(l1), a2=np.linalg.norm(l2))
                     )
    model.eval()
    print(model.report(), flush=True)

    # calculation the scaling factors as proposed by pybinding
    kpm = pb.kpm(model, kernel=pb.jackson_kernel(), matrix_format="CSR", optimal_size=False, interleaved=False)
    print(kpm.scaling_factors)

    # number of decomposition parts [nx,ny] in each direction of matrix.
    # This divides the lattice into various sections, each of which is calculated in parallel
    nx = ny = 5, 4
    # number of unit cells in each direction.
    lx = ly = 640, 512
    # energy scale
    e_min, e_max = -11.7, 9.3
    num_moments = 12000

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
                                       spectrum_range=[e_min, e_max])

    # specify calculation type
    calculation = kite.Calculation(configuration)
    calculation.dos(num_points=5000,
                    num_moments=num_moments,
                    num_random=1,
                    num_disorder=1)

    # configure the *.h5 file
    output_file = "dos_" + name[:-4] + "_{}_moments.h5".format(num_moments)
    kite.config_system(lattice, configuration, calculation, filename=output_file)

    # for generating the desired output from the generated HDF5-file, run
    # ../build/KITEx [output-file]
    # ../tools/build/KITE-tools [output-file]

    # returning the name of the created HDF5-file
    return output_file


if __name__ == "__main__":
    main()
