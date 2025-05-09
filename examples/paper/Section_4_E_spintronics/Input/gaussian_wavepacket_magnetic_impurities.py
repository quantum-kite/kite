""" Time-Evolution of Wave-Packet with magnetic impurities

    ##########################################################################
    #                         Copyright 2020/22, KITE                        #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in eV
    Lattice: Honeycomb (Kane-Mele and Rashba SOC)
    Disorder: Gaussian
    Configuration: Periodic boundary conditions, double precision,
                    given rescaling, size of the system flexible, with domain decomposition (nx=ny=4)
    Calculation type: Time Propagation of a wave-packet
    Last updated: 08/05/2025
"""
import kite
import numpy as np
import pybinding as pb
import matplotlib.pyplot as plt
import subprocess
import sys
import os
import h5py
import shutil
import multiprocessing
from numpy import linalg as LA

from pybinding.repository.graphene import a_cc, a
from pybinding.constants import hbar
from math import sqrt, pi

el_charge = 1.602 * 10 ** -19  #: [C] electron charge
t = -2.507  #: [eV] graphene nearest neighbor hopping
vf = 1.5 * a_cc * np.abs(t) / hbar * 10 ** -9  #: [m/s] Fermi velocity

energy = float(sys.argv[1])  # desired energy
lambda_R = float(sys.argv[2])  #: [eV] Rashba SOC, Phys. Rev. B 93, 155104
lambda_sv = float(sys.argv[3])

lambda_I_A = -lambda_sv
lambda_I_B = +lambda_sv

rashba_so = lambda_R * 2.0 * 1.0j / 3.0  #: [eV] constant and geometrical factors that will define Rashba SOC

km_so_A = lambda_I_A * 1j / (3 * sqrt(3))  #: [eV] constant and geometrical factors that will define intrinsic A SOC
km_so_B = lambda_I_B * 1j / (3 * sqrt(3))  #: [eV] constant and geometrical factors that will define intrinsic B SOC

# num points is the number of iterations in steps of deltaT
num_points = 500
deltaT = 2e-15  #: [fs] timestep
disorder_strength = float(sys.argv[4])
disorder_strength_and = float(sys.argv[9])
lx = ly = int(sys.argv[5])  # number of unit cells
tolerance = 0.001
# number of domains
nx = int(sys.argv[6])
ny = int(sys.argv[7])

calc_spin_x = int(sys.argv[8])

a1 = np.array([+a / 2, a * sqrt(3) / 2, 0])  #: [nm] unit cell vectors graphene
a2 = np.array([-a / 2, a * sqrt(3)/ 2, 0])

# std of the gaussian
sigma = ly * a / 16

posA = np.array([0, 0.0])
posB = np.array([0, -a_cc])


def honeycomb_lattice_with_SO(onsite=(0, 0)):
    """Make a honeycomb lattice with nearest neighbor hopping with SOC

    Parameters
    ----------
    onsite : tuple or list
        Onsite energy at different sublattices.
    """

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(
        a1=a1, a2=a2
    )
    print(a1, a2)
    # Add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('Aup', posA, onsite[0]),
        ('Bup', posB, onsite[1]),
        ('Adown', posA, onsite[0]),
        ('Bdown', posB, onsite[1])
    )

    # Add hoppings
    lat.add_hoppings(
        # ([f - i ], i , f )
        # inside the main cell, between which atoms, and the value
        ([0, 0], 'Adown', 'Bdown', t),
        ([0, 0], 'Aup', 'Bup', t),
        # between neighboring cells, between which atoms, and the value
        ([0, +1], 'Aup', 'Bup', t),
        ([0, +1], 'Adown', 'Bdown', t),

        ([+1, 0], 'Aup', 'Bup', t),
        ([+1, 0], 'Adown', 'Bdown', t)
    )

    if np.abs(lambda_R) > 0:
        lat.add_hoppings(
            # Rashba nearest neighbor, spin flip
            # inside the main cell, between which atoms, and the value
            ([0, 0], 'Aup', 'Bdown', -1.0 * rashba_so),  # delta1
            ([0, +1], 'Aup', 'Bdown', (+0.5 - sqrt(3) / 2 * 1j) * rashba_so),  # delta2
            ([+1, 0], 'Aup', 'Bdown', (+0.5 + sqrt(3) / 2 * 1j) * rashba_so),  # delta3

            ([0, 0], 'Adown', 'Bup', -1.0 * rashba_so),  # delta1
            ([0, +1], 'Adown', 'Bup', (+0.5 + sqrt(3) / 2 * 1j) * rashba_so),  # delta2
            ([+1, 0], 'Adown', 'Bup', (+0.5 - sqrt(3) / 2 * 1j) * rashba_so)  # delta3
        )

    if np.abs(lambda_sv) > 0:
        # Kane-Mele SOC, same spin next-nearest
        # between neighboring cells, between which atoms, and the value
        lat.add_hoppings(
            ([0, +1], 'Aup', 'Aup', -km_so_A),
            ([0, +1], 'Adown', 'Adown', +km_so_A),

            ([+1, 0], 'Aup', 'Aup', +km_so_A),
            ([+1, 0], 'Adown', 'Adown', -km_so_A),

            ([1, -1], 'Aup', 'Aup', -km_so_A),
            ([1, -1], 'Adown', 'Adown', +km_so_A),

            ([0, +1], 'Bup', 'Bup', +km_so_B),
            ([0, +1], 'Bdown', 'Bdown', -km_so_B),

            ([+1, 0], 'Bup', 'Bup', -km_so_B),
            ([+1, 0], 'Bdown', 'Bdown', +km_so_B),

            ([1, -1], 'Bup', 'Bup', +km_so_B),
            ([1, -1], 'Bdown', 'Bdown', -km_so_B),
        )

    return lat


def Hamiltonian_Graphene(k1, k2, lattice):
    b1, b2 = lattice.reciprocal_vectors()

    kp = k1 * b1[0:2] + k2 * b2[0:2]
    delta1 = posB - posA
    delta2 = delta1 + a2[0:2]
    delta3 = delta1 + a1[0:2]

    d1 = np.dot(kp, delta1)
    d2 = np.dot(kp, delta2)
    d3 = np.dot(kp, delta3)

    f = np.exp(1.j * d1) + np.exp(1.j * d2) + np.exp(1.j * d3)
    ham = np.array([[0, t * f], [t * np.conj(f), 0]])
    w, v = LA.eigh(ham)
    return [[w[0], v[:, 0]], [w[1], v[:, 1]]]


if __name__ == '__main__':

    use_disorder = True  # with or without the disorder
    lattice = honeycomb_lattice_with_SO()  # define the lattice
    conc = lambda_R / disorder_strength
    struc_disorder_one = kite.StructuralDisorder(lattice, concentration=conc)
    struc_disorder_two = kite.StructuralDisorder(lattice, concentration=conc)
    struc_disorder_one.add_structural_disorder(
        # in this way we can add onsite disorder in the form [unit cell], 'sublattice', value
        ([+0, +0], 'Aup', +disorder_strength),
        ([+0, +0], 'Adown', -disorder_strength),
    )
    # It is possible to add multiple different disorder type which should be forwarded to the export_lattice function
    # as a list.
    struc_disorder_two.add_structural_disorder(
        # in this way we can add onsite disorder in the form [unit cell], 'sublattice', value
        ([+0, +0], 'Bup', +disorder_strength),
        ([+0, +0], 'Bdown', -disorder_strength),
    )

    disorder = kite.Disorder(lattice)  # Disorder definition
    disorder_deviation_and = disorder_strength_and / (2 * sqrt(3))
    disorder.add_disorder(['Aup', 'Adown'], dis_type='uniform', mean_value=0, standard_deviation=disorder_deviation_and)
    disorder.add_disorder(['Bup', 'Bdown'], dis_type='uniform', mean_value=0, standard_deviation=disorder_deviation_and)

    b1, b2 = lattice.reciprocal_vectors()

    vf = 1.5 * a_cc * np.abs(t)  #: [ev * nm] redefined fermi velocity in natural units
    k = np.abs(energy / vf)  #: [nm^-1] find the wavevector length

    delta1 = 2 * k / np.linalg.norm(b1)
    delta2 = 2 * k / np.linalg.norm(b2)

    # you need to pick exact k point fitting the lattice which can be exactly compared with continuum

    ks = []
    spinors = []
    Kprm_point = np.array([+4 * np.pi/(3 * np.sqrt(3) * a_cc), 0])

    Ksec_point = np.array([-4 * np.pi/(3 * np.sqrt(3) * a_cc), 0])
    Kprm_point_rel = [+1/3, -1/3]
    Ksec_point_rel = [-1/3, +1/3]

    k1 = Kprm_point_rel[0] * b1[0:2] + Kprm_point_rel[1] * b2[0:2]
    k2 = Ksec_point_rel[0] * b1[0:2] + Ksec_point_rel[1] * b2[0:2]

    # pick the k vectors around the K point
    for kk in [Kprm_point_rel, Ksec_point_rel]:
        i1min, i1max = int((kk[0] - delta1) * lx), int((kk[0] + delta1) * lx)
        i2min, i2max = int((kk[1] - delta2) * ly), int((kk[1] + delta2) * ly)

        for i1 in range(i1min, i1max):
            for i2 in range(i2min, i2max):
                # x direction of motion ix = +iy
                # y direction of motion ix = -iy
                kx = i1 * 1. / lx
                ky = i2 * 1. / ly

                result = Hamiltonian_Graphene(kx, ky, lattice)
                for i in range(2):
                    [e, v] = result[i]
                    # print(e)
                    if np.abs(e - energy) < tolerance * energy:
                        k_vec = kx * b1[0:2] + ky * b2[0:2]
                        k_point = kk[0] * b1[0:2] + kk[1] * b2[0:2]

                        ks.append([kx, ky])
                        v = np.concatenate((v, v * int(calc_spin_x==True)))
                        v /= np.linalg.norm(v)
                        spinors.append(v)

    k_vector_rel = np.array(ks)
    spinor = np.array(spinors)
    print('size ', k_vector_rel, spinor)

    emin, emax = -3.03 * np.abs(t), 3.03 * np.abs(t)

    if use_disorder:
        emin -= 1.1 * np.abs(disorder_strength) / 2
        emax += 1.1 * np.abs(disorder_strength) / 2

    scale_a = (emax - emin) / 2
    configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=["periodic", "periodic"],
                                       is_complex=True, precision=1, spectrum_range=[emin, emax])

    # estimate the num moments
    num_moments = int(1.5 * 2.0 * scale_a * deltaT / hbar)

    calculation = kite.Calculation(configuration)
    calculation.gaussian_wave_packet(num_points=num_points, num_moments=num_moments, num_disorder=4,
                                     k_vector=k_vector_rel,
                                     spinor=spinor, width=sigma, timestep=scale_a * deltaT / hbar,
                                     mean_value=[int(lx / 2), int(ly / 2)])

    dirname = 'dis_{:.2f}_dis_res_{:.2f}_k_{:.2f}_sigma_{:.2f}_size_l_{:1d}'.format(disorder_strength_and * int(use_disorder == True), disorder_strength * int(use_disorder == True), k, sigma, lx)

    if os.path.isdir(dirname):
       shutil.rmtree(dirname)
    os.mkdir(dirname)
    os.chdir(dirname)

    # configure the *.h5 file
    name = 'gaussian_wavepacket.h5'

    if use_disorder:
        kite.config_system(lattice, configuration, calculation, filename=name, disorder=disorder,
disorder_structural=[struc_disorder_one, struc_disorder_two])
    else:
        kite.config_system(lattice, configuration, calculation, filename=name)

    pb.utils.tic()  # measure the time

    # # I have an alias for KITEx in zsh, useful for having everything in one script
    sp = subprocess.Popen(['/usr/bin/bash', '-i', '-c', 'KITEwavepacketLeib {}'.format(name)])
    sp.communicate()

    pb.utils.toc('Time for time evolution ')
