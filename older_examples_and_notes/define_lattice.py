import numpy as np
import pybinding as pb
import export_lattice as ex

from math import sqrt
# a = 0.24595  # [nm] unit cell length
# a_cc = 0.142  # [nm] carbon-carbon distance
# t = -2.8  # [eV] nearest neighbour hopping
from pybinding.repository.graphene import a_cc, a, t, t_nn

# INFO: This script contains some examples of predefined lattices


def square_lattice():
    d = 0.2  # [nm] unit cell length
    t = 1  # [eV] hopping energy

    # create a simple 2D lattice with vectors a1 and a2
    lat = pb.Lattice(a1=[d, 0], a2=[0, d])
    lat.add_sublattices(
        ('A', [0, 0])  # add an atom called 'A' at position [0, 0]
    )

    lat.add_hoppings(
        # (relative_index, from_sublattice, to_sublattice, energy)
        ([0, 1], 'A', 'A', t),
        ([1, 0], 'A', 'A', t)
    )

    return lat


def graphene_basic(onsite=(0, 0)):
    """Return the basic lattice specification for monolayer graphene with nearest neighbor"""

    # without importing anything
    a = 0.24595  # [nm] unit cell length
    a_cc = 0.142  # [nm] carbon-carbon distance
    t = -2.8  # [eV] nearest neighbour hopping

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(
        a1=[a, 0],
        a2=[a / 2, a / 2 * sqrt(3)]
    )

    # Add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, -a_cc / 2], onsite[0]),
        ('B', [0, a_cc / 2], onsite[1])
    )

    # Add hoppings
    lat.add_hoppings(
        # inside the main cell, between which atoms, and the value
        ([0, 0], 'A', 'B', t),
        # between neighboring cells, between which atoms, and the value
        ([1, -1], 'A', 'B', t),
        ([0, -1], 'A', 'B', t)
    )

    return lat


def monolayer_graphene(onsite=(0, 0), **kwargs):
    """Return the lattice specification for monolayer graphene"""

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(
        a1=[a, 0],
        a2=[a / 2, a / 2 * sqrt(3)]
    )
    # num of orbitals at each site is equal to the number of rows of the onsite object
    # the number of orbitals can be different, and the only limitation is that the hopping element between the two sites
    #  needs to be of the size num_orbitals_a x num_orbitals_b.
    num_orbitals_a = np.asarray(onsite[0]).shape[0]
    num_orbitals_b = np.asarray(onsite[1]).shape[0]

    # register name for hoppings
    lat.register_hopping_energies({
        't': kwargs.get('t', t * np.ones((num_orbitals_a, num_orbitals_b))),
        't_nn': kwargs.get('t_nn', t_nn * np.ones((num_orbitals_a, num_orbitals_b))),
        't_nnn': kwargs.get('t_nnn', 0.05 * np.ones((num_orbitals_a, num_orbitals_b))),
    })

    coord_A = [0, -a_cc / 2]
    coord_B = [0, a_cc / 2]

    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', coord_A, onsite[0]),
        ('B', coord_B, onsite[1])
    )

    lat.add_hoppings(
        # inside the main cell
        ([0, 0], 'A', 'B', 't'),
        # between neighboring cells
        ([1, -1], 'A', 'B', 't'),
        ([0, -1], 'A', 'B', 't')
    )

    lat.min_neighbors = 2

    disorder = ex.Disorder(lat)
    disorder.add_disorder('A', ['Gaussian', 'Uniform'], [1., 2.], [2., 5])
    disorder.add_disorder('B', ['Uniform', 'Gaussian'], [1., 3.], [2., 4.])

    return lat, disorder


def multi_orbital_example():
    lat = pb.Lattice([1, 0], [0, 1])
    lat.add_sublattices(
        ("A", [0.0, 0.0], 0.5),  # single-orbital: scalar
        ("B", [0.0, 0.2], [[1.5, 2],  # two-orbital: 2x2 Hermitian matrix
                           [2, 1.5]]),
        ("C", [0.3, 0.1], [[0, 0], [0, 0]]),  # two-orbital: zero onsite term
        ("D", [0.1, 0.0], [[4, 0, 0],  # three-orbital: only diagonal
                           [0, 5, 0],
                           [0, 0, 6]]),
        ("E", [0.2, 0.2], [44, 55, 66])  # three-orbital: only diagonal, terse notation
    )
    lat.add_hoppings(
        ([0, 1], "A", "A", 1.2),  # scalar
        ([0, 1], "B", "B", [[1, 2],  # 2x2
                            [3, 4]]),
        ([0, 0], "A", "B", [[1.2, 10.0]]),  # scalar

        ([0, 0], "B", "C", [[2, 0],  # 2x2
                            [1, 0]]),
        ([0, 0], "A", "D", [[1, 2, 3]]),  # 1x3
        ([0, 1], "D", "A", [[7],  # 3x1
                            [8],
                            [9]]),
        ([0, 0], "B", "D", [[1, 0, 0],  # 2x3
                            [2, 0, 3]])
    )

    lat.min_neighbors = 2

    disorder = ex.Disorder(lat)
    disorder.add_disorder('A', 'Gaussian', 1., 2.)
    disorder.add_disorder('B', ['Uniform', 'Gaussian'], [1., 3.], [2., 4.])

    return lat, disorder
