import numpy as np
import pybinding as pb

from math import sqrt
# a = 0.24595  # [nm] unit cell length
# a_cc = 0.142  # [nm] carbon-carbon distance
# t = -2.8  # [eV] nearest neighbour hopping
from pybinding.repository.graphene import a_cc, a, t, t_nn

pb.pltutils.use_style()


def monolayer_graphene(nearest_neighbors=1, onsite=(0, 0), **kwargs):
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
    num_orbitals_b = np.asarray(onsite[0]).shape[0]

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

    if nearest_neighbors >= 2:
        lat.add_hoppings(
            ([0, -1], 'A', 'A', 't_nn'),
            ([0, -1], 'B', 'B', 't_nn'),
            ([1, -1], 'A', 'A', 't_nn'),
            ([1, -1], 'B', 'B', 't_nn'),
            ([1, 0], 'A', 'A', 't_nn'),
            ([1, 0], 'B', 'B', 't_nn'),
        )

    if nearest_neighbors >= 3:
        lat.add_hoppings(
            [(1, -2), 'A', 'B', 't_nnn'],
            [(1, 0), 'A', 'B', 't_nnn'],
            [(-1, 0), 'A', 'B', 't_nnn'],
        )

    if nearest_neighbors >= 4:
        raise RuntimeError('No more')

    lat.min_neighbors = 2

    available_distr = {'rectangular', 'gaussian'}
    onsite_label = kwargs.get('onsite_labels', '[]')
    if len(onsite_label):
        if onsite_label not in available_distr:
            print('Available distributions are \n', list(available_distr))
            raise SystemExit('Selected distribution not available for the calculation! ')

    lat.onsite_labels = onsite_label
    lat.magnetic_field = kwargs.get('magnetic_field', 'False')
    return lat