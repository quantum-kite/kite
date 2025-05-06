""" Local Density of states of graphene

    ##########################################################################
    #                         Copyright 2022, KITE                           #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in eV
    Lattice: Honeycomb
    Configuration: Periodic boundary conditions, double precision,
                    manual rescaling, size of the system 128x128, with domain decomposition (nx=ny=1)
    Calculation type: Average DOS
    Last updated: 28/07/2022
"""

__all__ = ["main"]

import kite
import numpy as np
import pybinding as pb


def graphene_lattice(onsite=(0, 0)):
    """Return lattice specification for a honeycomb lattice with nearest neighbor hoppings"""

    # parameters
    a = 0.24595  # unit cell length
    a_cc = 0.142  # carbon-carbon distance
    t = 2.8  # eV

    # define lattice vectors
    a1 = a * np.array([1, 0])
    a2 = a * np.array([1 / 2, 1 / 2 * np.sqrt(3)])

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(a1=a1, a2=a2)

    # add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, -a_cc/2], onsite[0]),
        ('B', [0,  a_cc/2], onsite[1])
    )

    # Add hoppings
    lat.add_hoppings(
        # inside the main cell, between which atoms, and the value
        ([0, 0], 'A', 'B', -t),
        # between neighboring cells, between which atoms, and the value
        ([1, -1], 'A', 'B', -t),
        ([0, -1], 'A', 'B', -t)
    )
    return lat


def analyze_results(filename, lattice):
    """Analyze the results from the LDOS-calculation"""
    import matplotlib.pyplot as plt

    # extraxt the positions of the sublattices
    sublattices = lattice.sublattices
    sub_list = [(sublattices[sub_name].unique_id, sub_name) for sub_name in [*sublattices]]
    sub_list.sort(key=lambda x: x[0])
    sublattice_position = np.array([sublattices[sub_name[1]].position for sub_name in sub_list])

    # read the input file
    file_content = np.loadtxt(filename).T
    pos = file_content[:2, :] if file_content.shape[0] == 4 else file_content[:3, :]  # 2D or 3D lattice
    orbs = file_content[-2, :]
    values = file_content[-1, :]
    n_pos = len(values)
    n_orb = len(set(orbs))

    # fetch the positions of the sublattices
    orb_tmp = np.zeros((n_pos, n_orb))
    for i, i_o in enumerate(set(orbs)):
        orb_tmp[np.array(orbs) == i_o,  i] = 1

    # find the position of the atoms, unit cell + pos. of subl. in unit cell
    xs = np.dot(pos.T, lattice.vectors) + np.dot(orb_tmp, sublattice_position)

    # create colors for the different LDOS values
    values_max = np.max(np.abs(values))
    
    fig = plt.figure(figsize=(8, 5))
    ax_sctr = fig.add_axes([0.1, 0.05, 0.7, 0.9])
    sctr = ax_sctr.scatter(xs[:, 0], xs[:, 1], c=values, s=70, vmin=-values_max, vmax=values_max)
    plt.colorbar(sctr, cax=fig.add_axes([0.85, 0.1, 0.03, 0.8]))
    ax_sctr.set_aspect('equal')
    ax_sctr.set_title(filename)
    fig.savefig(filename[:-4]+".pdf")
    plt.close(fig)


def main(onsite=(0, 0)):
    """Prepare the input file for KITEx"""
    # load lattice
    lattice = graphene_lattice(onsite)

    # add StructuralDisorder
    node0 = [[+0, +0], 'A']
    node1 = [[+0, +0], 'B']
    node2 = [[-1, +0], 'B']
    node3 = [[-1, +1], 'B']

    timp = -0.0
    limp = -1.0

    struc_disorder = kite.StructuralDisorder(lattice, position=[[32, 32]])
    struc_disorder.add_structural_disorder(
        (*node0, *node1, timp),
        (*node0, *node2, timp),
        (*node0, *node3, timp),
        (*node0, limp)
    )

    # number of decomposition parts [nx,ny] in each direction of matrix.
    # This divides the lattice into various sections, each of which is calculated in parallel
    nx = ny = 1
    # number of unit cells in each direction.
    lx = ly = 128

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
        spectrum_range=[-2.8, 2.8]
    )

    # specify calculation type
    calculation = kite.Calculation(configuration)

    pos_matrix = []
    sub_matrix = []
    d1 = 32
    d2 = 32
    n = 64

    for di in range(-5, 5):
        for dj in range(-5, 5):
            pos_matrix.append([d1 + di, d2 + dj])
            pos_matrix.append([d1 + di, d2 + dj])
            sub_matrix.append('A')
            sub_matrix.append('B')

    calculation.ldos(
        energy=np.linspace(-0.4, 0.4, 41),
        num_moments=64,
        num_disorder=1,
        position=pos_matrix,
        sublattice=sub_matrix
    )

    # configure the *.h5 file
    output_file = "graphene_lattice_ldos-output.h5"
    kite.config_system(lattice, configuration, calculation, disorder_structural=struc_disorder,
                       filename=output_file)

    # for generating the desired output from the generated HDF5-file, run
    # ../build/KITEx graphene_lattice_ldos-output.h5
    # ../build/KITE-tools graphene_lattice_ldos-output.h5

    # returning the name of the created HDF5-file
    return output_file


if __name__ == "__main__":
    main()
