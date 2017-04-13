import numpy as np
import h5py as hp


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
        return {'bool': self._gaussian, 'width':  self._gauss_width, 'mean_value':  self._gauss_mean_value}


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
    vectors = lattice.vectors
    space_size = len(vectors)

    # get all the onsite values (as numbers or matrices) to the onsite array.
    onsite = []
    # get all positions to the position array.
    position = []
    # get number of orbitals at each atom.
    num_orbitals = []

    # iterate through all the sublattices and count num of orbitals
    # read onsite potential, and positions.
    for name, sub in lattice.sublattices.items():
        num_orbitals.append(np.asarray(sub.energy).shape[0])
        onsite.append(sub.energy)
        position.append(sub.position[0:space_size])
    onsite = np.array(onsite)
    position = np.array(position)

    hoppings_id = []
    hoppings = []
    # hopping_id matrix is made of the following data:
    # - each row is representing one hopping element,
    # - first space_size elements (1, 2, or 3) represent the relative lattice index with respect
    #   to the [0, 0, 0] lattice,
    # - next 2 numbers show the relative index of the atom in the lattice [0, 0, 0] from which the hopping comes
    #   and the relative index of the atom in the lattice given at the beginning of the row where the hopping goes.
    # - final two numbers are the number of rows and cols of hopping matrix which represent the number of orbitals of
    #   atom "from" and atom "to" between which the hopping happens.

    # hopping matrix with each row represents one hopping element (in case of the matrix the hopping is flatened
    # to fit the row) with properties given in the same row of matrix hopping_id.
    for name, hop in lattice.hoppings.items():
        hopping_energy = hop.energy
        for term in hop.terms:
            hoppings_id.append(np.hstack((term.relative_index[0:space_size].flatten(),
                                          term.from_id, term.to_id, hopping_energy.shape)))
            hoppings.append(hopping_energy.flatten())

    hoppings = np.array(hoppings)
    hoppings_id = np.array(hoppings_id)

    f = hp.File(filename, 'w')
    # hamiltonian is complex 1 or real 0
    complx = int(config.comp)
    f.create_dataset('IsComplex', data=complx)
    # precision of hamiltonian float, double, long double
    f.create_dataset('Precision', data=config.prec)
    # number of repetitions in each of the directions
    f.create_dataset('Length', data=config.leng, dtype='u4')
    # periodic boundary conditions, 0 - no, 1 - yes.
    bound = config.bound
    f.create_dataset('Boundaries', data=bound, dtype='u4')
    # number of divisions of the in each direction of hamiltonian. nx x ny = num_threads
    f.create_dataset('Divisions', data=config.div, dtype='u4')
    # space dimension of the lattice 1D, 2D, 3D
    f.create_dataset('Dim', data=space_size, dtype='u4')
    # lattice vectors. Size is same as DIM
    f.create_dataset('LattVector', data=vectors, dtype=np.float64)
    # position for each atom
    f.create_dataset('AtomPositions', data=position, dtype=np.float64)
    # total number of atom
    f.create_dataset('NumAtoms', data=position.shape[0], dtype='u4')
    # Hamiltonian group
    grp = f.create_group('Hamiltonian')
    # Hamiltonian group
    grp.create_dataset('HoppingsID', data=hoppings_id)
    # number of orbitals at each atom
    grp.create_dataset('NumOrbitals', data=num_orbitals)
    if complx:
        # hoppings
        grp.create_dataset('Hoppings', data=hoppings.astype(config.type))
        # onsite
        grp.create_dataset('Onsite', data=onsite.astype(config.type))
    else:
        # hoppings
        grp.create_dataset('Hoppings', data=hoppings.real.astype(config.type))
        # onsite potential at each atom
        grp.create_dataset('Onsite', data=onsite.real.astype(config.type))
    # labels for the onsite potential

    # Disorder group which has 'rectangular' and 'gaussian subgroups'
    # all of them have mean value and width
    grp_dis = f.create_group('Disorder')
    rec = modification.rectangular

    dis_type = grp_dis.create_group('rectangular')
    int(rec['bool'])

    dis_type.create_dataset('bool', data=int(rec['bool']), dtype='u4')
    dis_type.create_dataset('width', data=rec['width'])
    dis_type.create_dataset('mean_value', data=rec['mean_value'])

    rec = modification.gaussian
    dis_type = grp_dis.create_group('gaussian')
    dis_type.create_dataset('bool', data=int(rec['bool']), dtype='u4')
    dis_type.create_dataset('width', data=rec['width'])
    dis_type.create_dataset('mean_value', data=rec['mean_value'])

    # magnetic field
    if modification.magnetic_field:
        grp.create_dataset('MagneticField', data=int(modification.magnetic_field), dtype='u4')

    # Calculation
    grpc = f.create_group('Calculation')
    grpc.create_dataset('FunctionNum', data=calculation.number)
    grpc.create_dataset('NumMoments', data=calculation.moments)
    grpc.create_dataset('NumRandoms', data=calculation.randoms)
    grpc.create_dataset('NumDisorder', data=calculation.disorder)
    f.close()
