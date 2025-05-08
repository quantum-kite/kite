""" Bilayer Graphene with Rashba SOC (ARPES)

    ##########################################################################
    #                         Copyright 2020/2022, KITE                      #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in eV, Length in nm
    Lattice: honeycomb
    Configuration: Periodic boundary conditions, double precision,
                   manually scaling, size of the system 128x128, with domain decomposition (nx=ny=2)
    Disorder: Disorder class Uniform at different sublattices
    Calculation type:  one-particle spectral function of relevance to ARPES
    Last updated: 28/07/2022
"""

import pybinding as pb
from math import sqrt
import kite
import sys
from os import system as terminal


def bilayer_graphene_rashba():

  # Energy values
  t0 = -2.8                                 # [eV] intralayer hopping
  t1 = -0.4                                 # [eV] interlayer hopping
  onsite = 0.05                             # [eV] onsite potential
  lambda_R  = 0.075                         # [eV] Rashba SOC
  rashba_so = lambda_R * 2.0 * 1.0j / 3.0   # [eV] constant and geometrical factors that will define Rashba SOC

  # Length values
  a = 0.24595                               # [nm] unit cell length
  a_cc = 0.142                              # [nm] carbon-carbon distance
  c0 = 0.335                                # [nm] interlayer spacing
  lat = pb.Lattice(
      a1=[a/2,  a/2 * sqrt(3)], 
      a2=[-a/2, a/2 * sqrt(3)])

  lat.add_sublattices(
  ('A1u', [0,  -a_cc/2,   0], onsite),
  ('A1d', [0,  -a_cc/2,   0], onsite),
  ('B1u', [0,   a_cc/2,   0], onsite),
  ('B1d', [0,   a_cc/2,   0], onsite),
  ('B2u', [0, 3*a_cc/2, -c0],-onsite),
  ('B2d', [0, 3*a_cc/2, -c0],-onsite),
  ('A2u', [0,   a_cc/2, -c0],-onsite),
  ('A2d', [0,   a_cc/2, -c0],-onsite)
  )

  lat.register_hopping_energies({
  'gamma0': t0,  # [eV] intralayer
  'gamma1': t1,  # [eV] interlayer
  })

  lat.add_hoppings(
# layer 1, spin up
  ([ 0, 0], 'A1u', 'B1u', 'gamma0'),
  ([ 0, 1], 'A1u', 'B1u', 'gamma0'),
  ([ 1, 0], 'A1u', 'B1u', 'gamma0'),
# layer 1, spin down
  ([ 0, 0], 'A1d', 'B1d', 'gamma0'),
  ([ 0, 1], 'A1d', 'B1d', 'gamma0'),
  ([ 1, 0], 'A1d', 'B1d', 'gamma0'),
# layer 2, spin up
  ([ 0, 0], 'A2u', 'B2u', 'gamma0'),
  ([ 0, 1], 'A2u', 'B2u', 'gamma0'),
  ([ 1, 0], 'A2u', 'B2u', 'gamma0'),
# layer 2, spin down
  ([ 0, 0], 'A2d', 'B2d', 'gamma0'),
  ([ 0, 1], 'A2d', 'B2d', 'gamma0'),
  ([ 1, 0], 'A2d', 'B2d', 'gamma0'),
# interlayer, spin up and down
  ([ 0, 0], 'B1u', 'A2u', 'gamma1'),
  ([ 0, 0], 'B1d', 'A2d', 'gamma1')
  )

#Rashba SO coupling only in the first layer
  lat.add_hoppings(
    ([0, 0], 'A1u', 'B1d', -rashba_so*(-1.0)),  # delta1
    ([0, 1], 'A1u', 'B1d', -rashba_so*(0.5 - 1j*sqrt(3) / 2.0)),  # delta2
    ([1, 0], 'A1u', 'B1d', -rashba_so*(0.5 + 1j*sqrt(3) / 2.0)),  # delta3

    ([0, 0], 'A1d', 'B1u', -rashba_so*(-1.0)),  # delta1
    ([0, 1], 'A1d', 'B1u', -rashba_so*(0.5 + 1j*sqrt(3) / 2.0)),  # delta2
    ([1, 0], 'A1d', 'B1u', -rashba_so*(0.5 - 1j*sqrt(3) / 2.0))  # delta3
  )
#Rashba SO coupling only in the second layer
  lat.add_hoppings(
    ([0, 0], 'A2u', 'B2d', -rashba_so*(-1.0)),  # delta1
    ([0, 1], 'A2u', 'B2d', -rashba_so*(0.5 - 1j*sqrt(3) / 2.0)),  # delta2
    ([1, 0], 'A2u', 'B2d', -rashba_so*(0.5 + 1j*sqrt(3) / 2.0)),  # delta3

    ([0, 0], 'A2d', 'B2u', -rashba_so*(-1.0)),  # delta1
    ([0, 1], 'A2d', 'B2u', -rashba_so*(0.5 + 1j*sqrt(3) / 2.0)),  # delta2
    ([1, 0], 'A2d', 'B2u', -rashba_so*(0.5 - 1j*sqrt(3) / 2.0))  # delta3
  )

  return lat

def main():
    # Disorder parameters
    W = 0.4                                   # [eV] strength of Anderson disorder
    stddev = W/sqrt(12.0)                     # [eV] standard deviation of distribution
    mean = 0.0                                # [eV] average value of the distribution

    # Simulation parameters
    moments = 512
    nx = ny = 2
    lx = ly = 128
    LIM = 10


    lattice = bilayer_graphene_rashba()
    sublattices = ['A1u', 'A2u', 'A1d', 'A2d', 'B1u', 'B2u', 'B1d', 'B2d']

    # Add disorder to each sublattice individually
    disorder = kite.Disorder(lattice)
    for subl in sublattices:
        disorder.add_disorder(subl,'Uniform', mean, stddev)

    # Choose the path in k-space
    b1, b2 = lattice.reciprocal_vectors()
    Gamma = [0, 0] 
    K = 1 / 3 * (b1[0:2] - b2[0:2])
    M = b1[0:2]/2
    dk=0.5
    points = [Gamma, M, K, Gamma]
    k_path = pb.results.make_path(*points, step=dk)


    weights = [1 for i in range(8)]

    configuration = kite.Configuration(
            divisions=[nx, ny], 
            length=[lx, ly], 
            boundaries=["periodic", "periodic"], 
            is_complex=True, 
            precision=1, 
            spectrum_range=[-LIM,LIM])

    calculation = kite.Calculation(configuration)
    calculation.arpes(
            k_vector=k_path, 
            num_moments=moments, 
            weight=weights,
            num_disorder=1)

    output_file = "arpes-output.h5"
    kite.config_system(lattice, configuration, calculation, filename=output_file, disorder=disorder)
    return output_file


if __name__ == "__main__":
    hdf5_file = main() # generate the Configuration file

    if len(sys.argv) > 1 and sys.argv[1] == "complete":
        import run_all_examples as ra
        import process_arpes as pa
        ra.run_calculation(hdf5_file)
        ra.run_tools(hdf5_file, options="--ARPES -K green 0.1 -E -10 10 2048 -F 100")
        pa.process_arpes("arpes.dat")


