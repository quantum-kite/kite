import numpy as np
import h5py as hp

from scipy.sparse import coo_matrix


# Class that introduces Structural Disorder into the initially built lattice.
# The exported dataset StructuralDisorder has the following groups:
# - Concentration: concentration of disorder,
# - NumBondDisorder: number of bond changes,
# - NumOnsiteDisorder: number of onsite energy,
# The disorder is represented through nodes, where each node represents and maps a single orbital.
# - NumNodes: total number of orbitals included in the disorded,
# - NodePosition: orbital associated with the node,
# - NodeFrom, and NodeTo: bond from and to with 2 columns and NumBondDisorder rows, where the second column is complex
# conjugated hopping,
# - Hopping: values of the new hopping, with 2 columns and NumBondDisorder rows where the hoppings in different
# columns are conjugated values,
# - NodeOnsite: array of nodes that have onsite disorder,
# - U0: value of the onsite disorder.

class StructuralDisorder:
    def __init__(self, lattice, concentration):

        self._concentration = concentration
        self._num_bond_disorder_per_type = 0
        self._num_onsite_disorder_per_type = 0

        self._orbital_from = []
        self._orbital_to = []
        self._orbital_onsite = []

        self._idx_node = 0
        self._disorder_hopping = []
        self._disorder_onsite = []

        self._nodes_from = []
        self._nodes_to = []
        self._nodes_onsite = []

        self._orbital_vacancy = []

        self._num_nodes = 0
        self._nodes_map = dict()
        self._node_orbital = []
        num_orbitals = np.zeros(lattice.nsub, dtype=np.uint64)
        for name, sub in lattice.sublattices.items():
            # num of orbitals at each sublattice is equal to size of onsite energy
            num_energies = np.asarray(sub.energy).shape[0]
            num_orbitals[sub.alias_id] = num_energies
        self._num_orbitals_total = np.sum(np.asarray(num_orbitals))
        self._num_orbitals = np.asarray(num_orbitals)
        self._num_orbitals_before = np.cumsum(np.asarray(num_orbitals)) - num_orbitals
        self._lattice = lattice

    def add_vacancy(self, *disorder):

        num_vacancy_disorder = 0
        for dis in disorder:
            # check if it's just concentration or sublatt
            if len(dis) == 1:
                num_vacancy_disorder += 1
                self.add_local_vacancy_disorder(*dis)
            else:
                raise SystemExit('Vacancy disorder should be added in a form:'
                                 '\n sublattice name,'
                                 '\n or in a form of disorder onsite energy:'
                                 '\n ([rel. unit cell], sublattice_name, '
                                 'onsite energy)')

    def add_structural_disorder(self, *disorder):
        self._nodes_map = dict()

        num_bond_disorder_per_type = 0
        num_onsite_disorder_per_type = 0
        for dis in disorder:
            if len(dis) == 5:
                num_bond_disorder_per_type += 1
                self.add_local_bond_disorder(*dis)
            else:
                if len(dis) == 3:
                    num_onsite_disorder_per_type += 1
                    self.add_local_onsite_disorder(*dis)
                else:
                    raise SystemExit('Disorder should be added in a form of bond disorder:'
                                     '\n([rel. unit cell from], sublattice_from, [rel. unit cell to], sublattice_to, '
                                     'value),'
                                     '\n or in a form of disorder onsite energy:'
                                     '\n ([rel. unit cell], sublattice_name, '
                                     'onsite energy)')

        self._num_bond_disorder_per_type = num_bond_disorder_per_type
        self._num_onsite_disorder_per_type = num_onsite_disorder_per_type
        sorted_node_orb = sorted(self._nodes_map, key=lambda x: self._nodes_map[x])
        sorted_nodes = [self._nodes_map[x] for x in sorted_node_orb]

        sorted_dict = dict(zip(sorted_node_orb, sorted_nodes))
        self._nodes_map = sorted_dict
        self._node_orbital = sorted_node_orb

    def map_the_orbital(self, orb, nodes_map):
        idx_node = len(nodes_map)
        if not (orb in nodes_map):
            nodes_map[orb] = idx_node
            idx_node += 1
        self._nodes_map = nodes_map

        return idx_node

    def add_local_vacancy_disorder(self, sub):
        orbital_vacancy = []

        names, sublattices = zip(*self._lattice.sublattices.items())

        if sub not in names:
            raise SystemExit('Desired initial sublattice doesnt exist in the chosen lattice! ')

        indx = names.index(sub)
        lattice_sub = sublattices[indx]

        sub_id = lattice_sub.alias_id

        it = np.nditer(lattice_sub.energy, flags=['multi_index'])
        while not it.finished:
            orbit = self._num_orbitals_before[sub_id] + it.multi_index[0]
            if orbit not in orbital_vacancy:
                orbital_vacancy.append(orbit)
            it.iternext()

        self._orbital_vacancy.append(orbital_vacancy)

    def add_local_bond_disorder(self, relative_index_from, from_sub, relative_index_to, to_sub, hoppings):

        orbital_from = []
        orbital_to = []
        orbital_hop = []

        vectors = np.asarray(self._lattice.vectors)
        space_size = vectors.shape[0]

        distance_relative = np.asarray(relative_index_from) - np.asarray(relative_index_to)

        if np.linalg.norm(distance_relative) > space_size:
            raise SystemExit('Currently only the next nearest distances are supported, make the bond of the bond '
                             'disorder shorter! ')

        names, sublattices = zip(*self._lattice.sublattices.items())

        if from_sub not in names:
            raise SystemExit('Desired initial sublattice doesnt exist in the chosen lattice! ')
        if to_sub not in names:
            raise SystemExit('Desired final sublattice doesnt exist in the chosen lattice! ')

        indx_from = names.index(from_sub)
        lattice_sub_from = sublattices[indx_from]

        indx_to = names.index(to_sub)
        lattice_sub_to = sublattices[indx_to]

        from_sub_id = lattice_sub_from.alias_id
        to_sub_id = lattice_sub_to.alias_id

        nodes_map = self._nodes_map

        nodes_from = []
        nodes_to = []

        h = np.nditer(hoppings, flags=['multi_index'])
        while not h.finished:
            relative_move_from = np.dot(np.asarray(relative_index_from) + 1,
                                        3 ** np.linspace(0, space_size - 1, space_size, dtype=np.int32))
            relative_move_to = np.dot(np.asarray(relative_index_to) + 1,
                                      3 ** np.linspace(0, space_size - 1, space_size, dtype=np.int32))

            if isinstance(hoppings, np.ndarray):
                orb_from = int(relative_move_from +
                               (self._num_orbitals_before[from_sub_id] + h.multi_index[0]) * 3 ** space_size)
                orb_to = int(relative_move_to +
                             (self._num_orbitals_before[to_sub_id] + h.multi_index[1]) * 3 ** space_size)

                self.map_the_orbital(orb_from, nodes_map)
                self.map_the_orbital(orb_to, nodes_map)

                orbital_from.append(orb_from)
                orbital_to.append(orb_to)

                nodes_from.append(nodes_map[orb_from])
                nodes_to.append(nodes_map[orb_to])

                # conjugate
                orbital_from.append(orb_to)
                orbital_to.append(orb_from)

                nodes_from.append(nodes_map[orb_to])
                nodes_to.append(nodes_map[orb_from])

                orbital_hop.append(h[0])
                orbital_hop.append(np.conj(h[0]))

            else:
                orb_from = int(relative_move_from + self._num_orbitals_before[from_sub_id] * 3 ** space_size)
                orb_to = int(relative_move_to + self._num_orbitals_before[to_sub_id] * 3 ** space_size)

                self.map_the_orbital(orb_from, nodes_map)
                self.map_the_orbital(orb_to, nodes_map)

                orbital_from.append(orb_from)
                orbital_to.append(orb_to)

                nodes_from.append(nodes_map[orb_from])
                nodes_to.append(nodes_map[orb_to])

                # conjugate
                orbital_from.append(orb_to)
                orbital_to.append(orb_from)

                nodes_from.append(nodes_map[orb_to])
                nodes_to.append(nodes_map[orb_from])

                orbital_hop.append(h[0])
                orbital_hop.append((h[0].conjugate()))

            h.iternext()

        self._orbital_from.append(orbital_from)
        self._orbital_to.append(orbital_to)

        self._disorder_hopping.append(orbital_hop)

        self._nodes_from.append(nodes_from)
        self._nodes_to.append(nodes_to)

        if len(nodes_map) > self._num_nodes:
            self._num_nodes = len(nodes_map)

    def add_local_onsite_disorder(self, relative_index, sub, value):
        orbital_onsite = []
        orbital_onsite_en = []

        nodes_map = self._nodes_map

        vectors = np.asarray(self._lattice.vectors)
        space_size = vectors.shape[0]

        names, sublattices = zip(*self._lattice.sublattices.items())

        if sub not in names:
            raise SystemExit('Desired initial sublattice doesnt exist in the chosen lattice! ')

        indx_sub = names.index(sub)
        lattice_sub = sublattices[indx_sub]

        sub_id = lattice_sub.alias_id

        nodes_onsite = []

        h = np.nditer(value, flags=['multi_index'])
        while not h.finished:
            relative_move = np.dot(np.asarray(relative_index) + 1,
                                   3 ** np.linspace(0, space_size - 1, space_size, dtype=np.int32))

            if isinstance(value, np.ndarray):
                orb = int(relative_move + (self._num_orbitals_before[sub_id] + h.multi_index[0]) * 3 ** space_size)

                self.map_the_orbital(orb, nodes_map)

                orbital_onsite_en.append(h[0])

                nodes_onsite.append(nodes_map[orb])
                orbital_onsite.append(orb)
            else:
                orb = int(relative_move + self._num_orbitals_before[sub_id] * 3 ** space_size)

                self.map_the_orbital(orb, nodes_map)

                orbital_onsite_en.append(h[0])
                nodes_onsite.append(nodes_map[orb])
                orbital_onsite.append(orb)
            h.iternext()

        self._orbital_onsite.append(orbital_onsite)

        self._disorder_onsite.append(orbital_onsite_en)

        self._nodes_onsite.append(nodes_onsite)

        if len(nodes_map) > self._num_nodes:
            self._num_nodes = len(nodes_map)


# Class that introduces Disorder into the initially built lattice.
# The informations about the disorder are the type, mean value, and standard deviation. The function that you could use
# in the bulding of the lattice is add_disorder. The class method takes care of the shape of the disorder chosen (it
# needs to be same as the number of orbitals at a given atom), and takes care of the conversion to the c++ orbital-only
# format.
class Disorder:
    def __init__(self, lattice):
        # type of the disorder, can be 'Gaussian', 'Uniform' and 'Deterministic'.
        self._type = []
        # type_id of the disorder, can be 'Gaussian': 1, 'Uniform': 2 and 'Deterministic': 3.
        self._type_id = []
        # mean value of the disorder.
        self._mean = []
        # standard deviation of the disorder.
        self._stdv = []
        # orbital that has the chosen disorder.
        self._orbital = []

        num_orbitals = np.zeros(lattice.nsub, dtype=np.uint64)
        for name, sub in lattice.sublattices.items():
            # num of orbitals at each sublattice is equal to size of onsite energy
            num_energies = np.asarray(sub.energy).shape[0]
            num_orbitals[sub.alias_id] = num_energies
        self._num_orbitals_total = np.sum(np.asarray(num_orbitals))
        self._num_orbitals = np.asarray(num_orbitals)
        self._num_orbitals_before = np.cumsum(np.asarray(num_orbitals)) - num_orbitals
        self._lattice = lattice

    # class method that introduces the disorder to the lattice
    def add_disorder(self, sublattice, dis_type, mean_value, standard_deviation):
        if isinstance(dis_type, list):
            if isinstance(sublattice, list):
                for indx, name in enumerate(sublattice):
                    self.add_local_disorder(name, dis_type[indx], mean_value[indx], standard_deviation[indx])
            else:
                self.add_local_disorder(sublattice, dis_type, mean_value, standard_deviation)

        else:
            self.add_local_disorder(sublattice, [dis_type], [mean_value], [standard_deviation])

    def add_local_disorder(self, sublattice_name, dis_type, mean_value, standard_deviation):

        vectors = np.asarray(self._lattice.vectors)
        space_size = vectors.shape[0]

        names, sublattices = zip(*self._lattice.sublattices.items())

        if sublattice_name not in names:
            raise SystemExit('Desired sublattice doesnt exist in the chosen lattice! ')
        indx = names.index(sublattice_name)
        lattice_sub = sublattices[indx]
        size_orb = self._num_orbitals[lattice_sub.alias_id]

        hopping = {'relative_index': np.zeros(space_size, dtype=np.int32), 'from_id': lattice_sub.alias_id,
                   'to_id': lattice_sub.alias_id, 'mean_value': lattice_sub.energy}

        # number of orbitals before i-th sublattice, where is is the array index
        orbitals_before = self._num_orbitals_before

        orbital_from = self._orbital
        orbital_to = []
        orbital_dis_mean = self._mean
        orbital_dis_stdv = self._stdv
        orbital_dis_type = self._type
        orbital_dis_type_id = self._type_id

        dis_number = {'Gaussian': 1, 'Uniform': 2, 'Deterministic': 3, 'gaussian': 1, 'uniform': 2, 'deterministic': 3}
        for index, it in enumerate(mean_value):
            relative_move = np.dot(hopping['relative_index'] + 1,
                                   3 ** np.linspace(0, space_size - 1, space_size, dtype=np.int32))
            if len(mean_value) > 1:
                orbital_from.append(orbitals_before[hopping['from_id']] + index)
                orbital_to.append(relative_move + (orbitals_before[hopping['to_id']] + index) * 3 ** space_size)
            else:
                orbital_from.append(hopping['from_id'])
                orbital_to.append(relative_move + hopping['to_id'] * 3 ** space_size)
            orbital_dis_mean.append(it)
            orbital_dis_stdv.append(standard_deviation[index])
            orbital_dis_type.append(dis_type[index])
            if dis_type[index] in dis_number:
                orbital_dis_type_id.append(dis_number[dis_type[index]])
                if dis_type[index] == 'Deterministic' or dis_type[index] == 'deterministic':
                    if standard_deviation[index] != 0:
                        raise SystemExit(
                            'Standard deviation of deterministic disorder must be 0.')
            else:
                raise SystemExit(
                    'Disorder not present! Try between Gaussian, Deterministic, and Uniform case insensitive ')

        if not (all(np.asarray(i).shape == size_orb for i in [dis_type, mean_value, standard_deviation])):
            print('Shape of disorder', len(dis_type), len(mean_value), len(standard_deviation),
                  'is different than the number of orbitals at sublattice ', sublattice_name, 'which is', size_orb,
                  '\n')
            raise SystemExit('All parameters should have the same length! ')

        self._type = orbital_dis_type
        self._type_id = orbital_dis_type_id
        self._mean = orbital_dis_mean
        self._stdv = orbital_dis_stdv
        self._orbital = orbital_from


class Modification:
    def __init__(self, **kwargs):
        self._magnetic_field = kwargs.get('magnetic_field', False)

    @property
    def magnetic_field(self):  # magnetic_field:
        """Returns true if magnetic field is on, else off."""
        return self._magnetic_field


class Calculation:

    @property
    def get_dos(self):
        """Returns the requested DOS functions."""
        return self._dos

    @property
    def get_conductivity_dc(self):
        """Returns the requested DC conductivity functions."""
        return self._conductivity_dc

    @property
    def get_conductivity_optical(self):
        """Returns the requested optical conductivity functions."""
        return self._conductivity_optical

    @property
    def get_conductivity_optical_nonlinear(self):
        """Returns the requested nonlinear optical conductivity functions."""
        return self._conductivity_optical_nonlinear

    @property
    def get_singleshot_conductivity_dc(self):
        """Returns the requested singleshot DC conductivity functions."""
        return self._singleshot_conductivity_dc

    def __init__(self, configuration=None):

        if configuration is not None and not isinstance(configuration, Configuration):
            raise TypeError("You're forwaring a wrong type!")

        self._scaling_factor = configuration.energy_scale
        self._dos = []
        self._conductivity_dc = []
        self._conductivity_optical = []
        self._conductivity_optical_nonlinear = []
        self._singleshot_conductivity_dc = []

        self._avail_dir_full = {'xx': 0, 'yy': 1, 'zz': 2, 'xy': 3, 'xz': 4, 'yx': 3, 'yz': 5, 'zx': 4, 'zy': 5}
        self._avail_dir_sngl = {'xx': 0, 'yy': 1, 'zz': 2}

    def dos(self, num_points, num_moments, num_random):
        """Calculate the density of states as a function of energy

        Parameters
        ----------
        num_points : int
            Number of energy point inside the spectrum at which the DOS will be calculated.
        num_moments : int
            Number of polynomials in the Chebyshev expansion.
        num_random : int
            Number of random vectors to use for the stochastic evaluation of trace.
        """

        self._dos.append({'num_points': num_points, 'num_moments': num_moments, 'num_random': num_random})

    def conductivity_dc(self, direction, num_points, num_moments, num_random=1, temperature=0):
        """Calculate the density of states as a function of energy

        Parameters
        ----------
        direction : string
            direction in xyz coordinates along which the conductivity is calculated.
            Supports 'xx', 'yy', 'zz', 'xy', 'xz', 'yx', 'yz', 'zx', 'zy'.
        num_points : int
            Number of energy point inside the spectrum at which the DOS will be calculated.
        num_moments : int
            Number of polynomials in the Chebyshev expansion.
        num_random : int
            Number of random vectors to use for the stochastic evaluation of trace.
        temperature : float
            Value of the temperature at which we calculate the response.
        """
        if direction not in self._avail_dir_full:
            print('The desired direction is not available. Choose from a following set: \n',
                  self._avail_dir_full.keys())
            raise SystemExit('Invalid direction!')
        else:
            self._conductivity_dc.append(
                {'direction': self._avail_dir_full[direction], 'num_points': num_points, 'num_moments': num_moments,
                 'num_random': num_random, 'temperature': temperature / self._scaling_factor})

    def conductivity_optical(self, direction, num_points, num_moments, num_random=1, temperature=0):
        """Calculate the density of states as a function of energy

        Parameters
        ----------
        direction : string
            direction in xyz coordinates along which the conductivity is calculated.
            Supports 'xx', 'yy', 'zz', 'xy', 'xz', 'yx', 'yz', 'zx', 'zy'.
        num_points : int
            Number of energy point inside the spectrum at which the DOS will be calculated.
        num_moments : int
            Number of polynomials in the Chebyshev expansion.
        num_random : int
            Number of random vectors to use for the stochastic evaluation of trace.
        temperature : float
            Value of the temperature at which we calculate the response.
        """
        if direction not in self._avail_dir_full:
            print('The desired direction is not available. Choose from a following set: \n',
                  self._avail_dir_full.keys())
            raise SystemExit('Invalid direction!')
        else:
            self._conductivity_optical.append(
                {'direction': self._avail_dir_full[direction], 'num_points': num_points, 'num_moments': num_moments,
                 'num_random': num_random, 'temperature': temperature / self._scaling_factor})

    def conductivity_optical_nonlinear(self, direction, num_points, num_moments, num_random=1, temperature=0, **kwargs):
        """Calculate the density of states as a function of energy

        Parameters
        ----------
        direction : string
            direction in xyz coordinates along which the conductivity is calculated.
            Supports 'xx', 'yy', 'zz', 'xy', 'xz', 'yx', 'yz', 'zx', 'zy'.
        num_points : int
            Number of energy point inside the spectrum at which the DOS will be calculated.
        num_moments : int
            Number of polynomials in the Chebyshev expansion.
        num_random : int
            Number of random vectors to use for the stochastic evaluation of trace.
        temperature : float
            Value of the temperature at which we calculate the response.

            Optional parameters, forward special, a parameter that can simplify the calculation for some materials.
        """

        if direction not in self._avail_dir_full:
            print('The desired direction is not available. Choose from a following set: \n',
                  self._avail_dir_full.keys())
            raise SystemExit('Invalid direction!')
        else:
            special = kwargs.get('special', 0)

            self._conductivity_optical_nonlinear.append(
                {'direction': self._avail_dir_full[direction], 'num_points': num_points,
                 'num_moments': num_moments, 'num_random': num_random,
                 'temperature': temperature / self._scaling_factor, 'special': special})

    def singleshot_conductivity_dc(self, energy, direction, gamma, num_moments, num_random=1):
        """Calculate the density of states as a function of energy

        Parameters
        ----------
        energy : ndarray or float
            Array or a single value of energies at which singleshot_conductivity_dc will be calculated.
        direction : string
            direction in xyz coordinates along which the conductivity is calculated.
            Supports 'xx', 'yy', 'zz'.
        gamma : Float
            Parameter that affects the broadening of the kernel function.
        num_moments : int
            Number of polynomials in the Chebyshev expansion.
        num_random : int
            Number of random vectors to use for the stochastic evaluation of trace.
        """

        if direction not in self._avail_dir_sngl:
            print('The desired direction is not available. Choose from a following set: \n',
                  self._avail_dir_sngl.keys())
            raise SystemExit('Invalid direction!')
        else:
            self._singleshot_conductivity_dc.append(
                {'energy': np.array(energy) / self._scaling_factor, 'direction': self._avail_dir_sngl[direction],
                 'gamma': gamma / self._scaling_factor, 'num_moments': num_moments,
                 'num_random': num_random})


class Configuration:

    def __init__(self, divisions=(1, 1), length=(1, 1), boundaries=(False, False), is_complex=False, precision=1,
                 energy_scale=1):
        """Define basic parameters used in the calculation

       Parameters
       ----------
       divisions : int, tuple(int, int), tuple(int, int, int)
           Number of decomposition parts of the system.
       length : int, tuple(int, int), tuple(int, int, int)
           Number of unit cells in each direction.
       boundaries : int, tuple(int, int), tuple(int, int, int)
           Periodic boundary conditions each direction.
       is_complex : bool
           Boolean that reflects whether the type of Hamiltonian is complex or not.
       precision : int
            Integer which defines the precision of the number used in the calculation. Float - 0, double - 1,
            long double - 2.
        energy_scale : float
            Energy scale which defines the scaling factor of all the energy related parameters. The scaling is done
            automatically in the background after this definition.
       """

        self._energy_scale = energy_scale
        self._is_complex = int(is_complex)
        self._precision = precision
        self._divisions = divisions
        self._boundaries = np.asarray(boundaries).astype(int)

        self._length = length
        self._htype = np.float32
        self.set_type()

    def set_type(self, ):
        if self._is_complex == 0:
            if self._precision == 0:
                self._htype = np.float32
            elif self._precision == 1:
                self._htype = np.float64
            elif self._precision == 2:
                self._htype = np.float128
            else:
                raise SystemExit('Precision should be 0, 1 or 2')
        else:
            if self._precision == 0:
                self._htype = np.complex64
            elif self._precision == 1:
                self._htype = np.complex128
            elif self._precision == 2:
                self._htype = np.complex256

    @property
    def energy_scale(self):
        """Returns the energy scale of the hopping parameters."""
        return self._energy_scale

    @property
    def comp(self):  # -> is_complex:
        """Returns 0 if hamiltonian is real and 1 elsewise."""
        return self._is_complex

    @property
    def prec(self):  # -> precision:
        """Returns 0, 1, 2 if precision if float, double, and long double respectively."""
        return self._precision

    @property
    def div(self):  # -> divisions:
        """Returns the number of decomposed elements of matrix in x, y and/or z direction. Their product gives the total
        number of threads spawn."""
        return self._divisions

    @property
    def bound(self):  # -> boundaries:
        """Returns the boundary conditions in each direction, 0 - no boundary condtions, 1 - peridoc bc. """
        return self._boundaries

    @property
    def leng(self):  # -> length:
        """Return the number of unit cell repetitions in each direction. """
        return self._length

    @property
    def type(self):  # -> type:
        """Return the type of the Hamiltonian complex or real, and float, double or long double. """
        return self._htype


def export_lattice(lattice, config, calculation, modification, filename, **kwargs):
    """Export the lattice and related parameters to the *.h5 file

    Parameters
    ----------
    lattice : pb.Lattice
        Pybinding lattice object that carries the info about the unit cell vectors, unit cell cites, hopping terms and
        onsite energies.
    config : Configuration
        Configuration object, basic parameters defining size, precision, energy scale and number of decomposition parts
        in the calculation.
    calculation : Calculation
        Calculation object that defines the requested functions for the calculation.
    modification : Modification
        Modification object, has the magnetic field operator which can be True or False.
    filename : string
        Output filename.
    **kwargs: Optional arguments like Disorder or Disorder_structural.

    """

    # get the lattice vectors and set the size of space (1D, 2D or 3D) as the total number of vectors.
    disorder = kwargs.get('disorder', None)
    disorded_structural = kwargs.get('disorded_structural', None)

    vectors = np.asarray(lattice.vectors)
    space_size = vectors.shape[0]
    vectors = vectors[:, 0:space_size]

    # hamiltonian is complex 1 or real 0
    complx = int(config.comp)

    # iterate through all the sublattices and add onsite energies to hoppings list
    # count num of orbitals and read the atom positions.
    hoppings = []
    # get number of orbitals at each atom.
    num_orbitals = np.zeros(lattice.nsub, dtype=np.int64)
    # get all atom positions to the position array.
    position_atoms = np.zeros([lattice.nsub, space_size], dtype=np.float64)
    for name, sub in lattice.sublattices.items():
        # num of orbitals at each sublattice is equal to size of onsite energy
        num_energies = np.asarray(sub.energy).shape[0]
        num_orbitals[sub.alias_id] = num_energies
        # position_atoms is a list of vectors of size space_size
        position_atoms[sub.alias_id, :] = sub.position[0:space_size]
        # define hopping dict from relative hopping index from and to id (relative number of sublattice in relative
        # index lattice) and onsite
        hopping = {'relative_index': np.zeros(space_size, dtype=np.int32), 'from_id': sub.alias_id,
                   'to_id': sub.alias_id, 'hopping_energy': sub.energy}
        hoppings.append(hopping)

    # num_orbitals = np.asarray(num_orbitals)
    position_atoms = np.array(position_atoms)
    # repeats the positions of atoms based on the number of orbitals
    position = np.repeat(position_atoms, num_orbitals, axis=0)

    # iterate through all the hoppings and add hopping energies to hoppings list
    for name, hop in lattice.hoppings.items():
        hopping_energy = hop.energy
        for term in hop.terms:
            hopping = {'relative_index': term.relative_index[0:space_size], 'from_id': term.from_id,
                       'to_id': term.to_id, 'hopping_energy': hopping_energy}
            hoppings.append(hopping)
            # if the unit cell is [0, 0]
            if np.linalg.norm(term.relative_index[0:space_size]) == 0:
                hopping = {'relative_index': term.relative_index[0:space_size], 'from_id': term.to_id,
                           'to_id': term.from_id, 'hopping_energy': np.conj(hopping_energy)}
                hoppings.append(hopping)
            # if the unit cell [i, j] is different than [0, 0] and also -[i, j] hoppings with opposite direction
            if np.linalg.norm(term.relative_index[0:space_size]):
                hopping = {'relative_index': -term.relative_index[0:space_size], 'from_id': term.to_id,
                           'to_id': term.from_id, 'hopping_energy': np.conj(hopping_energy)}
                hoppings.append(hopping)

    orbital_from = []
    orbital_to = []
    orbital_hop = []
    # number of orbitals before i-th sublattice, where is is the array index
    orbitals_before = np.cumsum(num_orbitals) - num_orbitals
    # iterate through all hoppings, and define unique orbital hoppings
    # orbital_to in unit cell [i, j] is defined  as [i, j] x [1, 3] + relative_orbital_num*3**2 2D
    # orbital_to in unit cell [i, j, k] is defined  as [i, j, k] x [1, 3, 9] + relative_orbital_num*3**3 3D
    # relative index of orbital_from is unique as only hoppings from the orbitals in the initial unit cell are exported
    for h in hoppings:
        hopping_energy = h['hopping_energy']
        it = np.nditer(hopping_energy, flags=['multi_index'])
        while not it.finished:
            relative_move = np.dot(h['relative_index'] + 1,
                                   3 ** np.linspace(0, space_size - 1, space_size, dtype=np.int32))
            # if hopping_energy.size > 1:
            orbital_from.append(orbitals_before[h['from_id']] + it.multi_index[0])
            orbital_to.append(relative_move + (orbitals_before[h['to_id']] + it.multi_index[1]) * 3 ** space_size)
            # else:
            #     orbital_from.append(h['from_id'])
            #     orbital_to.append(relative_move + h['to_id'] * 3 ** space_size)

            orbital_hop.append(it[0] if complx else np.real(it[0]))

            it.iternext()

    # extract t - hoppings where each row corresponds to hopping from row number orbital and d - for each hopping it's
    # unique identifier
    t_list = []
    d_list = []
    # make a sparse matrix from orbital_hop, and (orbital_from, orbital_to) as it's easier to take nonzero hoppings from
    # sparse matrix
    matrix = coo_matrix((orbital_hop, (orbital_from, orbital_to)),
                        shape=(np.max(orbital_from) + 1, np.max(orbital_to) + 1))
    # num_hoppings is a vector where each value corresponds to num of hoppings from orbital equal to it's index
    num_hoppings = np.zeros(matrix.shape[0])
    # iterate through all rows of matrix, number of row = number of orbital from

    for i in range(matrix.shape[0]):
        # all hoppings from orbital i
        row_mat = matrix.getrow(i)
        # number of hoppings from orbital i
        num_hoppings[i] = row_mat.size

        t_list.append(row_mat.data)
        d_list.append(row_mat.indices)

    # fix the size of hopping and distance matrices, where the number of columns is max number of hoppings
    max_hop = int(np.max(num_hoppings))
    d = np.zeros((matrix.shape[0], max_hop))
    t = np.zeros((matrix.shape[0], max_hop), dtype=matrix.data.dtype)
    for i_row, d_row in enumerate(d_list):
        t_row = t_list[i_row]
        d[i_row, :len(d_row)] = d_row
        t[i_row, :len(t_row)] = t_row

    f = hp.File(filename, 'w')

    f.create_dataset('IS_COMPLEX', data=complx, dtype='u4')
    # precision of hamiltonian float, double, long double
    f.create_dataset('PRECISION', data=config.prec, dtype='u4')
    # number of repetitions in each of the directions
    f.create_dataset('L', data=config.leng, dtype='u4')
    # periodic boundary conditions, 0 - no, 1 - yes.
    bound = config.bound
    print('Periodic boundary conditions are set along x and y direction respec. to', bool(bound[0]), 'and ',
          bool(bound[1]), 'At the moment, no other boundary conditions are possible.')
    f.create_dataset('Boundaries', data=bound, dtype='u4')
    # number of divisions of the in each direction of hamiltonian. nx x ny = num_threads
    print('Chosen number of decomposition parts is:', config.div[0], 'x', config.div[1],
          'INFO: this product will correspond to the total number of threads. '
          'You should choose at most the number of processor cores you have.')
    f.create_dataset('Divisions', data=config.div, dtype='u4')
    # space dimension of the lattice 1D, 2D, 3D
    f.create_dataset('DIM', data=space_size, dtype='u4')
    # lattice vectors. Size is same as DIM
    f.create_dataset('LattVectors', data=vectors, dtype=np.float64)
    # position for each atom
    f.create_dataset('OrbPositions', data=position, dtype=np.float64)
    # total number of orbitals
    f.create_dataset('NOrbitals', data=np.sum(num_orbitals), dtype='u4')
    # scaling factor for the hopping parameters
    f.create_dataset('EnergyScale', data=config.energy_scale, dtype=np.float64)
    # Hamiltonian group
    grp = f.create_group('Hamiltonian')
    # Hamiltonian group
    grp.create_dataset('NHoppings', data=num_hoppings, dtype='u4')
    # distance
    grp.create_dataset('d', data=d, dtype='i4')

    if complx:
        # hoppings
        grp.create_dataset('Hoppings', data=t.astype(config.type) / config.energy_scale)
    else:
        # hoppings
        grp.create_dataset('Hoppings', data=t.real.astype(config.type) / config.energy_scale)

    # magnetic field
    if modification.magnetic_field:
        grp.create_dataset('MagneticField', data=int(modification.magnetic_field), dtype='u4')

    grp_dis = grp.create_group('Disorder')

    if disorder:
        grp_dis.create_dataset('OnsiteDisorderModelType', data=disorder._type_id, dtype=np.int32)
        grp_dis.create_dataset('OrbitalNum', data=disorder._orbital, dtype=np.int32)
        grp_dis.create_dataset('OnsiteDisorderMeanValue', data=disorder._mean / config.energy_scale, dtype=np.float64)
        grp_dis.create_dataset('OnsiteDisorderMeanStdv', data=disorder._stdv / config.energy_scale, dtype=np.float64)
    else:
        grp_dis.create_dataset('OnsiteDisorderModelType', (1, 0))
        grp_dis.create_dataset('OrbitalNum', (1, 0))
        grp_dis.create_dataset('OnsiteDisorderMeanValue', (1, 0))
        grp_dis.create_dataset('OnsiteDisorderMeanStdv', (1, 0))

    grp_dis_vac = grp.create_group('Vacancy')
    idx_vacancy = 0
    grp_dis = grp.create_group('StructuralDisorder')

    if disorded_structural:

        if isinstance(disorded_structural, list):
            num_dis = len(disorded_structural)
        else:
            num_dis = 1
            disorded_structural = [disorded_structural]

        for idx in range(num_dis):

            disorded_struct = disorded_structural[idx]

            num_orb_vac = len(disorded_struct._orbital_vacancy)
            if num_orb_vac > 0:
                grp_dis_type = grp_dis_vac.create_group('Type{val}'.format(val=idx_vacancy))
                grp_dis_type.create_dataset('Orbitals', data=np.asarray(disorded_struct._orbital_vacancy),
                                            dtype=np.int32)
                grp_dis_type.create_dataset('Concentration', data=disorded_struct._concentration,
                                            dtype=np.float64)
                grp_dis_type.create_dataset('NumOrbitals', data=num_orb_vac, dtype=np.int32)
                idx_vacancy += 1

            if disorded_struct._num_bond_disorder_per_type or disorded_struct._num_onsite_disorder_per_type:
                # Type idx
                grp_dis_type = grp_dis.create_group('Type{val}'.format(val=idx))
                # Concentration of this type
                grp_dis_type.create_dataset('Concentration', data=np.asarray(disorded_struct._concentration),
                                            dtype=np.float64)
                # Number of bond disorder
                grp_dis_type.create_dataset('NumBondDisorder',
                                            data=2 * np.asarray(disorded_struct._num_bond_disorder_per_type),
                                            dtype=np.int32)
                # Number of onsite disorder
                grp_dis_type.create_dataset('NumOnsiteDisorder',
                                            data=np.asarray(disorded_struct._num_onsite_disorder_per_type),
                                            dtype=np.int32)

                # Node of the bond disorder from
                grp_dis_type.create_dataset('NodeFrom', data=np.asarray(disorded_struct._nodes_from).flatten(),
                                            dtype=np.int32)
                # Node of the bond disorder to
                grp_dis_type.create_dataset('NodeTo', data=np.asarray(disorded_struct._nodes_to).flatten(),
                                            dtype=np.int32)
                # Node of the onsite disorder
                grp_dis_type.create_dataset('NodeOnsite', data=np.asarray(disorded_struct._nodes_onsite),
                                            dtype=np.int32)

                # Num nodes
                grp_dis_type.create_dataset('NumNodes', data=disorded_struct._num_nodes, dtype=np.int32)
                # Orbital mapped for this node
                grp_dis_type.create_dataset('NodePosition', data=np.asarray(disorded_struct._node_orbital),
                                            dtype=np.uint32)

                # Onsite disorder energy
                grp_dis_type.create_dataset('U0',
                                            data=np.asarray(disorded_struct._disorder_onsite).real.astype(
                                                config.type)) / config.energy_scale
                # Bond disorder hopping
                disorder_hopping = disorded_struct._disorder_hopping
                if complx:
                    # hoppings
                    grp_dis_type.create_dataset('Hopping',
                                                data=np.asarray(disorder_hopping).astype(
                                                    config.type).flatten() / config.energy_scale)
                else:
                    # hoppings
                    grp_dis_type.create_dataset('Hopping',
                                                data=np.asarray(disorder_hopping).real.astype(
                                                    config.type).flatten() / config.energy_scale)

    # Calculation function defined with num_moments, num_random vectors, and num_disorder etc. realisations
    grpc = f.create_group('Calculation')
    if calculation.get_dos:
        grpc_p = grpc.create_group('dos')

        for single_dos in calculation.get_dos:
            grpc_p.create_dataset('NumMoments', data=single_dos['num_moments'], dtype=np.int32)
            grpc_p.create_dataset('NumRandoms', data=single_dos['num_random'], dtype=np.int32)
            grpc_p.create_dataset('NumPoints', data=single_dos['num_points'], dtype=np.int32)

    if calculation.get_conductivity_dc:
        grpc_p = grpc.create_group('conductivity_dc')

        moments, random, point, temp, direction = [], [], [], [], []
        for single_cond_dc in calculation.get_conductivity_dc:
            moments.append(single_cond_dc['num_moments'])
            random.append(single_cond_dc['num_random'])
            point.append(single_cond_dc['num_points'])
            temp.append(single_cond_dc['temperature'])
            direction.append(single_cond_dc['direction'])

        grpc_p.create_dataset('NumMoments', data=np.asarray(moments), dtype=np.int32)
        grpc_p.create_dataset('NumRandoms', data=np.asarray(random), dtype=np.int32)
        grpc_p.create_dataset('NumPoints', data=np.asarray(point), dtype=np.int32)
        grpc_p.create_dataset('Temperature', data=np.asarray(temp), dtype=np.float64)
        grpc_p.create_dataset('Direction', data=np.asarray(direction), dtype=np.int32)

    if calculation.get_conductivity_optical:
        grpc_p = grpc.create_group('conductivity_optical')

        moments, random, point, temp, direction = [], [], [], [], []
        for single_cond_opt in calculation.get_conductivity_optical:
            moments.append(single_cond_opt['num_moments'])
            random.append(single_cond_opt['num_random'])
            point.append(single_cond_opt['num_points'])
            temp.append(single_cond_opt['temperature'])
            direction.append(single_cond_opt['direction'])

        grpc_p.create_dataset('NumMoments', data=np.asarray(moments), dtype=np.int32)
        grpc_p.create_dataset('NumRandoms', data=np.asarray(random), dtype=np.int32)
        grpc_p.create_dataset('NumPoints', data=np.asarray(point), dtype=np.int32)
        grpc_p.create_dataset('Temperature', data=np.asarray(temp), dtype=np.float64)
        grpc_p.create_dataset('Direction', data=np.asarray(direction), dtype=np.int32)

    if calculation.get_conductivity_optical_nonlinear:
        grpc_p = grpc.create_group('conductivity_optical_nonlinear')

        moments, random, point, temp, direction, special = [], [], [], [], [], []
        for single_cond_opt_non in calculation.get_conductivity_optical_nonlinear:
            moments.append(single_cond_opt_non['num_moments'])
            random.append(single_cond_opt_non['num_random'])
            point.append(single_cond_opt_non['num_points'])
            temp.append(single_cond_opt_non['temperature'])
            direction.append(single_cond_opt_non['direction'])
            special.append(single_cond_opt_non['special'])

        grpc_p.create_dataset('NumMoments', data=np.asarray(moments), dtype=np.int32)
        grpc_p.create_dataset('NumRandoms', data=np.asarray(random), dtype=np.int32)
        grpc_p.create_dataset('NumPoints', data=np.asarray(point), dtype=np.int32)
        grpc_p.create_dataset('Temperature', data=np.asarray(temp), dtype=np.float64)
        grpc_p.create_dataset('Direction', data=np.asarray(direction), dtype=np.int32)
        grpc_p.create_dataset('Special', data=np.asarray(special), dtype=np.int32)

    if calculation.get_singleshot_conductivity_dc:
        grpc_p = grpc.create_group('singleshot_conductivity_dc')

        moments, random, energies, gamma, direction = [], [], [], [], []

        for single_singlshot_cond in calculation.get_singleshot_conductivity_dc:
            moments.append(single_singlshot_cond['num_moments'])
            random.append(single_singlshot_cond['num_random'])
            energies.append(single_singlshot_cond['energy'])
            gamma.append(single_singlshot_cond['gamma'])
            direction.append(single_singlshot_cond['direction'])

        grpc_p.create_dataset('NumMoments', data=np.asarray(moments), dtype=np.int32)
        grpc_p.create_dataset('NumRandoms', data=np.asarray(random), dtype=np.int32)
        grpc_p.create_dataset('Energy', data=np.asarray(energies), dtype=np.float64)
        grpc_p.create_dataset('Gamma', data=np.asarray(gamma), dtype=np.float64)
        grpc_p.create_dataset('Direction', data=np.asarray(direction), dtype=np.int32)

    f.close()
