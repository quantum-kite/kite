""" Honeycomb lattice with Kane-Mele and Rashba SOC model
"""
import kite
import numpy as np
import pybinding as pb
import sys
import h5py

from pybinding.repository.graphene import a_cc, a
from pybinding.constants import hbar
from math import sqrt, pi
from scipy.spatial import cKDTree

el_charge = 1.602 * 10 ** -19  #: [C] electron charge
t = -2.507  #: [eV] graphene nearest neighbor hopping
delta = 0.54e-3   #: [eV] onsite mass term, Phys. Rev. B 93, 155104
vf = 1.5 * a_cc * np.abs(t) / hbar * 10e-9  #: [m/s] Fermi velocity

lambda_R = 0.56e-3  #: [eV] Rashba SOC, Phys. Rev. B 93, 155104
lambda_I_A = -1.22e-3  #: [eV] intrisic sublattice A SOC, Phys. Rev. B 93, 155104
lambda_I_B = +1.16e-3  #: [eV] intrisic sublattice B SOC, Phys. Rev. B 93, 155104

lambda_PIA_A = -2.69e-3  #: [eV] PIA sublattice A SOC, Phys. Rev. B 93, 155104
lambda_PIA_B = -2.54e-3  #: [eV] PIA sublattice B SOC, Phys. Rev. B 93, 155104

rashba_so = lambda_R * 2.0 * 1.0j / 3.0  #: [eV] constant and geometrical factors that will define Rashba SOC

km_so_A = lambda_I_A * 1j / (3 * sqrt(3))  #: [eV] constant and geometrical factors that will define intrinsic A SOC
km_so_B = lambda_I_B * 1j / (3 * sqrt(3))  #: [eV] constant and geometrical factors that will define intrinsic B SOC

pia_so_A = lambda_PIA_A * 2.0 * 1.0j / 3.0  #: [eV] constant and geometrical factors that will define PIA A SOC
pia_so_B = lambda_PIA_B * 2.0 * 1.0j / 3.0  #: [eV] constant and geometrical factors that will define PIA B SOC

a1 = [a / 2 * sqrt(3), a / 2, 0]  #: [nm] unit cell vectors graphene
a2 = [a / 2 * sqrt(3), -a / 2, 0]


def honeycomb_lattice_with_SO(onsite=(0, 0)):
    lat = pb.Lattice(
        a1=a1, a2=a2
    )

    # Add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('Aup', [a - a_cc / 2.0, 0.0], onsite[0]),
        ('Bup', [a + a_cc / 2.0, 0.0], onsite[1]),
        ('Adown', [a - a_cc / 2.0, 0.0], onsite[0]),
        ('Bdown', [a + a_cc / 2.0, 0.0], onsite[1]))

    # Add hoppings
    lat.add_hoppings(
        # inside the main cell, between which atoms, and the value
        ([0, 0], 'Adown', 'Bdown', t),
        ([0, 0], 'Aup', 'Bup', t),
        # between neighboring cells, between which atoms, and the value
        ([0, -1], 'Aup', 'Bup', t),
        ([0, -1], 'Adown', 'Bdown', t),
        ([-1, 0], 'Aup', 'Bup', t),
        ([-1, 0], 'Adown', 'Bdown', t),

        # Rashba nearest neighbor, spin flip

        # inside the main cell, between which atoms, and the value
        ([0, 0], 'Aup', 'Bdown', -1j * rashba_so),
        ([0, 0], 'Adown', 'Bup', +1j * rashba_so),
        # between neighboring cells, between which atoms, and the value
        ([0, -1], 'Aup', 'Bdown', (+sqrt(3) / 2 + 0.5j) * rashba_so),
        ([0, -1], 'Adown', 'Bup', (+sqrt(3) / 2 - 0.5j) * rashba_so),
        ([-1, 0], 'Aup', 'Bdown', (-sqrt(3) / 2 + 0.5j) * rashba_so),
        ([-1, 0], 'Adown', 'Bup', (-sqrt(3) / 2 - 0.5j) * rashba_so),

        # Kane-Mele SOC, same spin next-nearest

        # between neighboring cells, between which atoms, and the value
        ([0, -1], 'Aup', 'Aup', +km_so_A),
        ([0, -1], 'Adown', 'Adown', -km_so_A),
        ([-1, 0], 'Aup', 'Aup', -km_so_A),
        ([-1, 0], 'Adown', 'Adown', +km_so_A),
        ([1, -1], 'Aup', 'Aup', -km_so_A),
        ([1, -1], 'Adown', 'Adown', +km_so_A),

        ([0, -1], 'Bup', 'Bup', -km_so_B),
        ([0, -1], 'Bdown', 'Bdown', +km_so_B),
        ([-1, 0], 'Bup', 'Bup', +km_so_B),
        ([-1, 0], 'Bdown', 'Bdown', -km_so_B),
        ([1, -1], 'Bup', 'Bup', +km_so_B),
        ([1, -1], 'Bdown', 'Bdown', -km_so_B),

        # PIA SO next-nearest

        ([0, -1], 'Aup', 'Adown', +(+sqrt(3) / 2 + 1.5j) * pia_so_A),
        ([0, -1], 'Adown', 'Aup', +(+sqrt(3) / 2 - 1.5j) * pia_so_A),
        ([-1, 0], 'Aup', 'Adown', +(-sqrt(3) / 2 + 1.5j) * pia_so_A),
        ([-1, 0], 'Adown', 'Aup', +(-sqrt(3) / 2 - 1.5j) * pia_so_A),
        ([1, -1], 'Aup', 'Adown', +(sqrt(3)) * pia_so_A),
        ([1, -1], 'Adown', 'Aup', +(sqrt(3)) * pia_so_A),

        ([0, -1], 'Bup', 'Bdown', +(+sqrt(3) / 2 + 1.5j) * pia_so_B),
        ([0, -1], 'Bdown', 'Bup', +(+sqrt(3) / 2 - 1.5j) * pia_so_B),
        ([-1, 0], 'Bup', 'Bdown', +(-sqrt(3) / 2 + 1.5j) * pia_so_B),
        ([-1, 0], 'Bdown', 'Bup', +(-sqrt(3) / 2 - 1.5j) * pia_so_B),
        ([1, -1], 'Bup', 'Bdown', +(sqrt(3)) * pia_so_B),
        ([1, -1], 'Bdown', 'Bup', +(sqrt(3)) * pia_so_B),
    )

    return lat


lattice = honeycomb_lattice_with_SO((delta / 2, -delta / 2))  # define the lattice

disorder = kite.Disorder(lattice)  # Disorder definition
disorder_strength = 0.5 * t
disorder_deviation = disorder_strength / sqrt(3)
disorder.add_disorder('Aup', dis_type='uniform', mean_value=0, standard_deviation=disorder_deviation)
disorder.add_disorder('Adown', dis_type='uniform', mean_value=0, standard_deviation=disorder_deviation)
disorder.add_disorder('Bup', dis_type='uniform', mean_value=0, standard_deviation=disorder_deviation)
disorder.add_disorder('Bdown', dis_type='uniform', mean_value=0, standard_deviation=disorder_deviation)

lx = 128
ly = 128

b1, b2 = lattice.reciprocal_vectors()

nx = 1
ny = 1

emin, emax = -7.61, 7.61
scale_a = (emax - emin) / 2

emin -= np.abs(disorder_strength)
emax += np.abs(disorder_strength)

configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=["periodic", "periodic"], is_complex=True, precision=1, spectrum_range=[emin, emax])

Gamma = [0, 0]
M = [np.pi/np.sqrt(3), np.pi/3.0]
K = [4*np.pi/np.sqrt(27), 0.0]
moments = 16
points = [Gamma, M, K, Gamma]

k_path = pb.results.make_path(*points, step=0.3)

calculation_arpes = kite.Calculation(configuration)
calculation_arpes.arpes(k_vector=k_path, weight=[1.5, 0.5, 1.3, 0.3], num_moments=moments, num_disorder=1)
kite.config_system(lattice, configuration, calculation_arpes, filename='config.h5')
