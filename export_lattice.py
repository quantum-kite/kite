import numpy as np
import h5py as hp

from scipy.sparse import coo_matrix


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
    def add_disorder(self, sublattice, type, mean_value, standard_deviation):
        if isinstance(type, list):
            for name in sublattice:
                self.add_local_disorder(name, type, mean_value, standard_deviation)
        else:
            for name in sublattice:
                self.add_local_disorder(name, [type], [mean_value], [standard_deviation])

    def add_local_disorder(self, sublattice_name, type, mean_value, standard_deviation):

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
            orbital_dis_type.append(type[index])
            if type[index] in dis_number:
                orbital_dis_type_id.append(dis_number[type[index]])
                if type[index] == 'Deterministic' or type[index] == 'deterministic':
                    if standard_deviation[index] != 0:
                        raise SystemExit(
                            'Standard deviation of deterministic disorder must be 0.')
            else:
                raise SystemExit(
                    'Disorder not present! Try between Gaussian, Deterministic, and Uniform case insensitive ')

        if not (all(np.asarray(i).shape == size_orb for i in [type, mean_value, standard_deviation])):
            print('Shape of disorder', len(type), len(mean_value), len(standard_deviation),
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
    def __init__(self, fname='DOS', num_moments=1, num_random=1, num_disorder=1):

        num_f = len(fname) if fname else None
        if not (all(len(i) == num_f for i in [num_moments, num_random, num_disorder])):
            print('Number of different functions is different than the entered parameters, num_moments, num_randoms, '
                  'or num_disorder \n')
            raise SystemExit('All parameters should have the same length! ')

        self._number = []
        self._num_moments = []
        self._fname = []
        self._num_random = []
        self._num_disorder = []

        available_functions = {'DOS', 'CondXX', 'CondXY', 'OptCond', 'SpinCond'}
        fun_number = {'DOS': 1, 'CondXX': 2, 'CondXY': 3, 'OptCond': 4, 'SpinCond': 5}

        for f in fname:
            if f not in available_functions:
                print('Available functions are \n', list(available_functions))
                raise SystemExit('Function not available for the calculation! ')

            if len(fname) > 1:
                self._number.append(fun_number[f])
                idx = fname.index(f)
                self._fname.append(f)
                self._num_moments.append(num_moments[idx])
                self._num_random.append(num_random[idx])
                self._num_disorder.append(num_disorder[idx])
            else:
                self._number.append(fun_number[f])
                self._fname.append(f)
                self._num_moments.append(num_moments)
                self._num_random.append(num_random)
                self._num_disorder.append(num_disorder)

    @property
    def fname(self):  # function name:
        """Returns the desired function name."""
        return self._fname

    @property
    def number(self):  # -> function number:
        """Returns the predefined number of desired function."""
        return self._number

    @property
    def name(self):  # -> function name:
        """Returns the desired function name."""
        return self._fname

    @property
    def moments(self):  # -> function name:
        """Returns the number of moments given for the calc."""
        return self._num_moments

    @property
    def randoms(self):  # -> function name:
        """Returns the number of random vectors given for the calc."""
        return self._num_random

    @property
    def disorder(self):  # -> function name:
        """Returns the number of disorder realisations given for the calc."""
        return self._num_disorder


class Configuration:
    def __init__(self, divisions=(1, 1), length=None, boundaries=(False, False), is_complex=False, precision=1):
        self._is_complex = int(is_complex)
        self._precision = precision
        self._divisions = divisions
        self._boundaries = np.asarray(boundaries).astype(int)
        self._length = length
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
    def comp(self):  # -> is_complex:
        """Returns 0 if hamiltonian is real and 1 elsewise."""
        return self._is_complex

    @property
    def prec(self):  # -> precision:
        """Returns 0, 1, 2 if precision if float, double, and long double respectively."""
        return self._precision

    @property
    def div(self):  # -> divisions:
        """Returns the number of decomposed elements of matrix in x and y direction. Product of the two gives the total
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
    # get the lattice vectors and set the size of space (1D, 2D or 3D) as the total number of vectors.
    disorder = kwargs.get('disorder', 0)
    vectors = np.asarray(lattice.vectors)
    space_size = vectors.shape[0]
    vectors = vectors[:, 0:space_size]

    # hamiltonian is complex 1 or real 0
    complx = int(config.comp)

    # get all atom positions to the position array.
    position_atoms = []
    # get number of orbitals at each atom.
    num_orbitals = []

    # iterate through all the sublattices and add onsite energies to hoppings list
    # count num of orbitals and read the atom positions.
    hoppings = []
    num_orbitals = np.zeros(lattice.nsub, dtype=np.int64)
    for name, sub in lattice.sublattices.items():
        # num of orbitals at each sublattice is equal to size of onsite energy
        num_energies = np.asarray(sub.energy).shape[0]
        num_orbitals[sub.alias_id] = num_energies
        # position_atoms is a list of vectors of size space_size
        position_atoms.append(sub.position[0:space_size])
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
    t = np.zeros((matrix.shape[0], max_hop))
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
    # Hamiltonian group
    grp = f.create_group('Hamiltonian')
    # Hamiltonian group
    grp.create_dataset('NHoppings', data=num_hoppings, dtype='u4')
    # distance
    grp.create_dataset('d', data=d, dtype='u4')

    if complx:
        # hoppings
        grp.create_dataset('Hoppings', data=t.astype(config.type))
    else:
        # hoppings
        grp.create_dataset('Hoppings', data=t.real.astype(config.type))

    # magnetic field
    if modification.magnetic_field:
        grp.create_dataset('MagneticField', data=int(modification.magnetic_field), dtype='u4')

    # TODO: Change this comparison
    if disorder != 0:
        grp_dis = grp.create_group('Disorder')
        grp_dis.create_dataset('OnsiteDisorderModelType', data=disorder._type_id)
        grp_dis.create_dataset('OrbitalNum', data=disorder._orbital)
        grp_dis.create_dataset('OnsiteDisorderMeanValue', data=disorder._mean)
        grp_dis.create_dataset('OnsiteDisorderMeanStdv', data=disorder._stdv)

    # Calculation function defined with num_moments, num_random vectors, and num_disorder realisations
    grpc = f.create_group('Calculation')
    grpc.create_dataset('FunctionNum', data=np.asarray(calculation.number), dtype=np.int32)
    grpc.create_dataset('NumMoments', data=np.asarray(calculation.moments), dtype=np.int32)
    grpc.create_dataset('NumRandoms', data=np.asarray(calculation.randoms), dtype=np.int32)
    grpc.create_dataset('NumDisorder', data=np.asarray(calculation.disorder), dtype=np.int32)
    f.close()
