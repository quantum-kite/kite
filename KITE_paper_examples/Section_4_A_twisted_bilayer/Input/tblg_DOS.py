"""       
        ##############################################################################      
        #                        KITE | Release  1.0                                 #      
        #                                                                            #      
        #                        Kite home: quantum-kite.com                         #           
        #                                                                            #      
        #  Developed by: Simao M. Joao, Joao V. Lopes, Tatiana G. Rappoport,         #       
        #  Misa Andelkovic, Lucian Covaci, Aires Ferreira, 2018-2019                 #      
        #                                                                            #      
        ##############################################################################      
"""
import os
import pybinding as pb
import numpy as np
import re
import copy
import kite

from math import sqrt



def load_ovito_lattice(name):
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

    space_size = int(sqrt(len(vec)))

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


# WARNING: the calculation ran with parameters given here will run for at least 72h 
# load an xyz file, relaxed or unrelaxed
name = 'relaxed_tblg_1.050.xyz'


# load coordinates, atom types and lattice vectors
x_coord, y_coord, z_coord, atom_type, l1, l2, l3 = load_ovito_lattice(name)
l1 = np.array(l1[0:2]) / 10.
l2 = np.array(l2[0:2]) / 10.
positions = np.column_stack((x_coord, y_coord, z_coord))
positions_to = copy.deepcopy(positions)

# load a PB lattice
complete_lattice = pb.load('lattice_' + name[:-4])
num_x = 1
num_y = 1
print('Done loading ', flush=True)
print('Making the model', flush=True)

# check for the bonds
model = pb.Model(complete_lattice,
                 pb.force_double_precision(),
                 pb.translational_symmetry(a1=num_x * np.linalg.norm(l1), a2=num_y * np.linalg.norm(l2))
                 )
model.eval()
print(model.report(), flush=True)
kpm = pb.kpm(model, kernel=pb.jackson_kernel(), matrix_format="CSR", optimal_size=False, interleaved=False)
print(kpm.scaling_factors)

nx = 5
ny = 4
e_min = - 11.7
e_max = + 9.3

num_a1 = 640
num_a2 = 512

num_moments = 12000

# make config object which caries info about
# - the number of decomposition parts [nx, ny],
# - lengths of structure [lx, ly]
# - boundary conditions, setting True as periodic boundary conditions, and False elsewise,
# - info if the exported hopping and onsite data should be complex,
# - info of the precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
configuration = kite.Configuration(divisions=[nx, ny], length=[num_a1, num_a2], boundaries=[True, True],
                                   is_complex=False, precision=0, spectrum_range=[e_min, e_max])

calculation = kite.Calculation(configuration)
calculation.dos(num_points=5000, num_moments=num_moments, num_random=1, num_disorder=1)

# export the lattice from the lattice object, config and calculation object and the name of the file
# the disorder is optional. If there is disorder in the lattice for now it should be given separately
kite.config_system(complete_lattice, configuration, calculation,
                   filename='dos_' + name[:-4] + '_{}_moments.h5'.format(num_moments))
