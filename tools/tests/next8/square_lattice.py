import kite
import numpy as np
import sys
import pybinding as pb


L = int(sys.argv[1])
N = int(sys.argv[2])
W = float(sys.argv[3])
# B = int(sys.argv[4])
# W = 1
# B = 1

def square_lattice(onsite=0):
    a1 = np.array([1, 0])
    a2 = np.array([0, 1])

    lat = pb.Lattice(a1=a1, a2=a2)
    lat.add_sublattices(('A', [0, 0], onsite))
    lat.add_hoppings(   ([1, 0], 'A', 'A', -1),
                        ([0, 1], 'A', 'A', -1))
    return lat


std_dev = W*W/12.0
mean    = 0

lattice = square_lattice()
disorder = kite.Disorder(lattice)
disorder.add_disorder('A', 'Uniform', mean, std_dev)

# min_B = 2.019368814847448*(2048.0/L)*B
# mod   = kite.Modification(magnetic_field=min_B)

nx = ny = 8
lx = ly = L
energy = -3.5
# name = 'square_lattice_WL_B' + str(B) + '_W' + str(W) + '_N' + str(N) + '_L' + str(L) + '.h5'
name = 'square_lattice_WL' + '_W' + str(W) + '_N' + str(N) + '_L' + str(L) + '.h5'

b1, b2 = lattice.reciprocal_vectors()
Gamma = [0, 0]
R = b1[0:2] + b2[0:2]
M = b1[0:2]
K = [0.23*np.pi, 0]

points = [K, K]
k_path = pb.results.make_path(*points)
weights = 1.0


SCALE = 4.02 + W/2.0
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=["periodic", "periodic"], is_complex=True, precision=1, spectrum_range=[-SCALE, SCALE])
calculation = kite.Calculation(configuration)
# calculation.singleshot_conductivity_dc(energy=energy, num_moments=N, num_random=1, num_disorder=1, direction='yy', eta=scattering)
# calculation.dos(num_points=1000, num_moments=N, num_random=1, num_disorder=1)
calculation.arpes(k_vector=k_path, num_moments=N, weight=weights,num_disorder=1)
kite.config_system(lattice, configuration, calculation, filename=name, disorder=disorder)
