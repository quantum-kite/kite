import numpy as np
import pybinding as pb
import kite

from pybinding.constants import hbar
from numpy import linalg as LA

t = 2.8  # eV
a_cc = 0.142  # nm
a = a_cc * np.sqrt(3)

a1_2atom = np.array([a / 2 * np.sqrt(3), +a / 2, 0])  #: [nm] unit cell vectors graphene
a2_2atom = np.array([a / 2 * np.sqrt(3), -a / 2, 0])

a1_4atom = a1_2atom + a2_2atom
a2_4atom = a1_2atom - a2_2atom

# std of the gaussian

posA = np.array([-a_cc / 2.0, 0.0])
posB = np.array([+a_cc / 2.0, 0.0])


def Hamiltonian_Graphene(k1, k2, lattice):
    b1, b2 = lattice.reciprocal_vectors()

    kp = k1 * b1[0:2] + k2 * b2[0:2]
    delta1 = posB - posA
    delta2 = delta1 + a2_2atom[0:2]
    delta3 = delta1 + a1_2atom[0:2]

    d1 = np.dot(kp, delta1)
    d2 = np.dot(kp, delta2)
    d3 = np.dot(kp, delta3)

    f = np.exp(1.j * d1) + np.exp(1.j * d2) + np.exp(1.j * d3)
    ham = np.array([[0, t * f], [t * np.conj(f), 0]])
    w, v = LA.eigh(ham)
    return [[w[0], v[:, 0]], [w[1], v[:, 1]]]


def graphene(onsite=(0, 0)):
    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(
        a1=a1_2atom, a2=a2_2atom
    )

    # Add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', posA, onsite[0]),
        ('B', posB, onsite[1]),
    )

    # Add hoppings
    lat.add_hoppings(
        # ([f - i ], i , f )
        # inside the main cell, between which atoms, and the value
        ([0, 0], 'A', 'B', -t),
        # between neighboring cells, between which atoms, and the value
        ([0, -1], 'A', 'B', -t),
        ([-1, 0], 'A', 'B', -t),
    )

    return lat


lx = 256
ly = 256

# STD of the Gaussian
sigma = lx * 3 * a_cc / 120

domain_decompose_1 = 1
domain_decompose_2 = 1

energy = 0.1
delta = 0.2

lattice_gr = graphene((-delta, delta))
b1, b2 = lattice_gr.reciprocal_vectors()

vf = 1.5 * a_cc * np.abs(t)  #: [ev * nm] redefined fermi velocity in natural units
k = np.abs(energy / vf)  #: [nm^-1] find the wavevector length

delta1 = 2 * k / np.linalg.norm(b1)
delta2 = 2 * k / np.linalg.norm(b2)

# you need to pick exact k point fitting the lattice which can be exactly compared with continuum
ks = []
spinors = []

Kprm_point_rel = [+1 / 3, +2 / 3]
Ksec_point_rel = [+2 / 3, +1 / 3]

k1 = Kprm_point_rel[0] * b1[0:2] + Kprm_point_rel[1] * b2[0:2]
k2 = Ksec_point_rel[0] * b1[0:2] + Ksec_point_rel[1] * b2[0:2]

theta = np.pi/2

tolerance_angle = 0.2
tolerance = 0.5

# pick the k vectors around the K point
# for kk in [Kprm_point_rel, Ksec_point_rel]:
for kk in [Kprm_point_rel, Ksec_point_rel]:
    i1min, i1max = int((kk[0] - delta1) * lx), int((kk[0] + delta1) * lx)
    i2min, i2max = int((kk[1] - delta2) * ly), int((kk[1] + delta2) * ly)

    for i1 in range(i1min, i1max):
        for i2 in range(i2min, i2max):
            # x direction of motion ix = +iy
            # y direction of motion ix = -iy
            kx = i1 * 1. / lx
            ky = i2 * 1. / ly

            result = Hamiltonian_Graphene(kx, ky, lattice_gr)
            for i in range(2):
                [e, v] = result[i]
                if np.abs(e - energy) < tolerance * energy:
                    k_vec = kx * b1[0:2] + ky * b2[0:2]
                    k_point = kk[0] * b1[0:2] + kk[1] * b2[0:2]

                    k_rel = k_vec - k_point
                    rel_angle = np.arctan2(k_rel[1], k_rel[0])
                    # print('rel angle', rel_angle)
                    # select only one direction
                    if np.abs(rel_angle - theta) < tolerance_angle:
                        print('rel_angle', rel_angle)
                        ks.append([kx, ky])
                        # v = np.concatenate((v * electron_like))  # for both electrons and holes
                        v /= np.linalg.norm(v)
                        spinors.append(v)

k_vector_rel = np.array(ks)
spinor = np.array(spinors)
print('k vector: ', k_vector_rel, '\nspinor: ', spinor)

emin, emax = -3.03 * np.abs(t) - delta, 3.03 * np.abs(t) + delta
use_disorder = False

configuration = kite.Configuration(divisions=[domain_decompose_1, domain_decompose_2], length=[lx, ly],
                                   boundaries=[True, True], is_complex=True, precision=1, spectrum_range=[emin, emax])

scale_a = (emax - emin) / 2
deltaT = 1e-15  #: [fs] timestep
num_points = 200

# estimate the num moments
num_moments = int(1.5 * 2.0 * scale_a * deltaT / hbar)

calculation = kite.Calculation(configuration)
calculation.gaussian_wave_packet(num_points=num_points, num_moments=num_moments, num_disorder=1,
                                 k_vector=k_vector_rel,
                                 spinor=spinor, width=sigma, timestep=scale_a * deltaT / hbar,
                                 mean_value=[int(lx / 2), int(ly / 2)])


name = 'gaussian_planewave.h5'
kite.config_system(lattice_gr, configuration, calculation, filename=name)
