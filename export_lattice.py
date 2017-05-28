import numpy as np
import h5py as hp

from scipy.sparse import coo_matrix

class Modifications:
    def __init__(self, **kwargs):
        self._magnetic_field = kwargs.get('magnetic_field', False)
        self._rec_width = 0
        self._rec_mean_value = 0
        self._gauss_width = 0
        self._gauss_mean_value = 0

        available_distr = ['rectangular', 'gaussian']
        disorder_realisation = kwargs.get('disorder', {})
        if not disorder_realisation:
            self._rectangular = False
            self._gaussian = False
        elif disorder_realisation['name'] == available_distr[0]:
            self._rectangular = True
            self._gaussian = False
            self._rec_width = disorder_realisation['width']
            self._rec_mean_value = disorder_realisation['mean_value']
        elif disorder_realisation['name'] == available_distr[1]:
            self._rectangular = False
            self._gaussian = True
            self._gauss_width = disorder_realisation['width']
            self._gauss_mean_value = disorder_realisation['mean_value']
        else:
            print('Available distributions are \n', list(available_distr))
            raise SystemExit('Selected distribution not available for the calculation! ')

    @property
    def magnetic_field(self):  # -> magnetic_field:
        """Returns true if magnetic field is on, else off."""
        return self._magnetic_field

    @property
    def rectangular(self):  # -> rectangular disorder realisation:
        """Returns True if rectangular disorder realisation, it's width and mean value."""
        return {'bool': self._rectangular, 'width': self._rec_width, 'mean_value': self._rec_mean_value}

    @property
    def gaussian(self):  # -> gausian disorder realisation:
        """Returns True if gaussian disorder realisation, it's width and mean value."""
        return {'bool': self._gaussian, 'width': self._gauss_width, 'mean_value': self._gauss_mean_value}


class Calculation:
    def __init__(self, fname='DOS', num_moments=1, num_random=1, num_disorder=1):
        available_functions = {'DOS', 'CondXX', 'CondXY', 'OptCond', 'SpinCond'}
        fun_number = {'DOS': 1, 'CondXX': 2, 'CondXY': 3, 'OptCond': 4, 'SpinCond': 5}
        if fname not in available_functions:
            print('Available functions are \n', list(available_functions))
            raise SystemExit('Function not available for the calculation! ')
        self._number = fun_number[fname]
        self._fname = fname
        self._num_moments = num_moments
        self._num_random = num_random
        self._num_disorder = num_disorder

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


class Config:
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


def export_lattice(lattice, config, calculation, modification, filename):
    # get the lattice vectors and set the size of space (1D, 2D or 3D) as the total number of vectors.

    vectors = np.asarray(lattice.vectors)
    space_size = vectors.shape[0]
    vectors = vectors[:, 0:space_size]

    # hamiltonian is complex 1 or real 0
    complx = int(config.comp)

    # get all positions to the position array.
    position = []
    # get number of orbitals at each atom.
    num_orbitals = []

    # iterate through all the sublattices and add onsite energies to hoppings list
    # count num of orbitals and read the positions.
    hoppings = []
    for name, sub in lattice.sublattices.items():
        # num of orbitals at each sublattice is equal to size of onsite energy
        num_energies = np.asarray(sub.energy).shape[0]
        num_orbitals.append(num_energies)
        # position is a list of vectors of size space_size
        position.append(sub.position[0:space_size])
        # define hopping dict from relative hopping index from and to id (relative number of sublattice in relative
        # index lattice) and onsite energy
        hopping = {'relative_index': np.zeros(space_size, dtype=np.int32), 'from_id': sub.alias_id,
                   'to_id': sub.alias_id, 'hopping_energy': sub.energy}
        hoppings.append(hopping)

    num_orbitals = np.asarray(num_orbitals)
    position = np.array(position)

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
    # number of orbitals before this sublattice
    orbitals_before = np.cumsum(num_orbitals) - num_orbitals

    # iterate through all hoppings, and define unique orbital hoppings
    # orbital_to in unit cell [i, j] is defined  as [i, j] x [1, 3] + relative_orbital_num*3**2 2D
    # orbital_to in unit cell [i, j, k] is defined  as [i, j, k] x [1, 3, 9] + relative_orbital_num*3**3 3D
    # relative index of orbital_from is unique as only hoppings from the orbitals in the initial unit cell are exported
    for h in hoppings:
        hopping_energy = h['hopping_energy']
        it = np.nditer(hopping_energy, flags=['multi_index'])

        while not it.finished:
            relative_move = np.dot(h['relative_index'] + 1, 3**np.linspace(0, space_size - 1, space_size, dtype=np.int32))
            if hopping_energy.size > 1:
                orbital_from.append(orbitals_before[h['from_id']] + it.multi_index[0])
                orbital_to.append(relative_move + (orbitals_before[h['to_id']] + it.multi_index[1])*3**space_size)
            else:
                orbital_from.append(h['from_id'])
                orbital_to.append(relative_move + h['to_id']*3**space_size)

            orbital_hop.append(it[0] if complx else np.real(it[0]))
            it.iternext()

    # extract t - hoppings where each row coresponds to hopping from row number orbital and d - for each hopping it's
    # unique identifier
    t_list = []
    d_list = []
    # make a sparse matrix from orbital_hop, and (orbital_from, orbital_to) as it's easier to take nonzero hoppings from
    # sparse matrix
    matrix = coo_matrix((orbital_hop, (orbital_from, orbital_to)), shape=(np.max(orbital_from)+1, np.max(orbital_to)+1))
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

    t = np.asarray(t_list)
    d = np.asarray(d_list)

    f = hp.File(filename, 'w')

    f.create_dataset('IS_COMPLEX', data=complx, dtype='u4')
    # precision of hamiltonian float, double, long double
    f.create_dataset('PRECISION', data=config.prec, dtype='u4')
    # number of repetitions in each of the directions
    f.create_dataset('L', data=config.leng, dtype='u4')
    # periodic boundary conditions, 0 - no, 1 - yes.
    bound = config.bound
    f.create_dataset('Boundaries', data=bound, dtype='u4')
    # number of divisions of the in each direction of hamiltonian. nx x ny = num_threads
    f.create_dataset('Divisions', data=config.div, dtype='u4')
    # space dimension of the lattice 1D, 2D, 3D
    f.create_dataset('DIM', data=space_size, dtype='u4')
    # lattice vectors. Size is same as DIM
    f.create_dataset('LattVectors', data=vectors, dtype=np.float64)
    # position for each atom
    f.create_dataset('OrbPositions', data=position, dtype=np.float64)
    # total number of atom
    f.create_dataset('NOrbitals', data=position.shape[0], dtype='u4')
    # Hamiltonian group
    grp = f.create_group('Hamiltonian')
    # Hamiltonian group
    grp.create_dataset('NHoppings', data=num_hoppings, dtype='u4')
    # number of orbitals at each atom
    grp.create_dataset('NumOrbitals', data=num_orbitals, dtype='u4')
    # distance
    grp.create_dataset('d', data=d)

    if complx:
        # hoppings
        grp.create_dataset('Hoppings', data=t.astype(config.type))
        # # onsite
        # grp.create_dataset('Onsite', data=onsite.astype(config.type))
    else:
        # hoppings
        grp.create_dataset('Hoppings', data=t.real.astype(config.type))
        # # onsite potential at each atom
        # grp.create_dataset('Onsite', data=onsite.real.astype(config.type))

    # Disorder group which has 'rectangular' and 'gaussian subgroups'
    # all of them have mean value and width
    grp_dis = f.create_group('Disorder')

    rec = modification.rectangular
    dis_type = grp_dis.create_group('rectangular')
    dis_type.create_dataset('bool', data=int(rec['bool']), dtype='u4')
    dis_type.create_dataset('width', data=rec['width'], dtype=np.float64)
    dis_type.create_dataset('mean_value', data=rec['mean_value'], dtype=np.float64)

    rec = modification.gaussian
    dis_type = grp_dis.create_group('gaussian')
    dis_type.create_dataset('bool', data=int(rec['bool']), dtype='u4')
    dis_type.create_dataset('width', data=rec['width'], dtype=np.float64)
    dis_type.create_dataset('mean_value', data=rec['mean_value'], dtype=np.float64)

    # magnetic field
    if modification.magnetic_field:
        grp.create_dataset('MagneticField', data=int(modification.magnetic_field), dtype='u4')

    # Calculation function defined with num_moments, num_random vectors, and num_disorder realisations
    grpc = f.create_group('Calculation')
    grpc.create_dataset('FunctionNum', data=calculation.number, dtype='u4')
    grpc.create_dataset('NumMoments', data=calculation.moments, dtype='u4')
    grpc.create_dataset('NumRandoms', data=calculation.randoms, dtype='u4')
    grpc.create_dataset('NumDisorder', data=calculation.disorder, dtype='u4')
    f.close()
