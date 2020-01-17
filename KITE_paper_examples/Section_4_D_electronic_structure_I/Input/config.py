import pybinding as pb
import matplotlib.pyplot as plt
from math import sqrt, pi
import kite
import sys

# Energy values

t0 = -2.8                                 # [eV] intralayer hopping
t1 = -0.4                                 # [eV] interlayer hopping
onsite = 0.05                             # [eV] onsite potential
lambda_R  = 0.075                          # [eV] Rashba SOC


#onsite = 0.2*t0                             # [eV] onsite potential
#lambda_R  = 0.15*t0                          # [eV] Rashba SOC
rashba_so = lambda_R * 2.0 * 1.0j / 3.0   # [eV] constant and geometrical factors that will define Rashba SOC
# Length values
a = 0.24595                               # [nm] unit cell length
a_cc = 0.142                              # [nm] carbon-carbon distance
c0 = 0.335                                # [nm] interlayer spacing

def bilayer_graphene_rashba():

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
  ([1, 0], 'A1u', 'B1u', 'gamma0'),
# layer 1, spin down
  ([ 0, 0], 'A1d', 'B1d', 'gamma0'),
  ([ 0, 1], 'A1d', 'B1d', 'gamma0'),
  ([1, 0], 'A1d', 'B1d', 'gamma0'),
# layer 2, spin up
  ([ 0, 0], 'A2u', 'B2u', 'gamma0'),
  ([ 0, 1], 'A2u', 'B2u', 'gamma0'),
  ([1, 0], 'A2u', 'B2u', 'gamma0'),
# layer 2, spin down
  ([ 0, 0], 'A2d', 'B2d', 'gamma0'),
  ([ 0, 1], 'A2d', 'B2d', 'gamma0'),
  ([1, 0], 'A2d', 'B2d', 'gamma0'),
# interlayer, spin up and down
  ([ 0,  0], 'B1u', 'A2u', 'gamma1'),
  ([ 0,  0], 'B1d', 'A2d', 'gamma1')
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


lattice = bilayer_graphene_rashba()


W = float(sys.argv[13])
stddev = W/sqrt(12.0)
disorder = kite.Disorder(lattice)
disorder.add_disorder('A1u','Uniform', 0.0, stddev)
disorder.add_disorder('A2u','Uniform', 0.0, stddev)
disorder.add_disorder('A1d','Uniform', 0.0, stddev)
disorder.add_disorder('A2d','Uniform', 0.0, stddev)

disorder.add_disorder('B1u','Uniform', 0.0, stddev)
disorder.add_disorder('B2u','Uniform', 0.0, stddev)
disorder.add_disorder('B1d','Uniform', 0.0, stddev)
disorder.add_disorder('B2d','Uniform', 0.0, stddev)

b1, b2 = lattice.reciprocal_vectors()
Gamma = [0, 0]
K = 1 / 3 * (b1[0:2] - b2[0:2])
M = b1[0:2]/2

GM = M - Gamma
GM1 = Gamma + 1.0*GM/6.0
GM2 = Gamma + 2.0*GM/6.0
GM3 = Gamma + 3.0*GM/6.0
GM4 = Gamma + 4.0*GM/6.0
GM5 = Gamma + 5.0*GM/6.0

MK = K - M
MK1 = M + 1.0*MK/6.0
MK2 = M + 2.0*MK/6.0
MK3 = M + 3.0*MK/6.0
MK4 = M + 4.0*MK/6.0
MK5 = M + 5.0*MK/6.0

KG  = Gamma - K
KG1 = K + 1.0/6.0*KG
KG2 = K + 2.0/6.0*KG
KG3 = K + 3.0/6.0*KG
KG4 = K + 4.0/6.0*KG
KG5 = K + 5.0/6.0*KG

#points = [Gamma, M, K, Gamma]
points = []

points.append([Gamma,   GM1])
points.append([GM1, GM2])
points.append([GM2, GM3])
points.append([GM3, GM4])
points.append([GM4, GM5])
points.append([GM5, M  ])

points.append([M,   MK1])
points.append([MK1, MK2])
points.append([MK2, MK3])
points.append([MK3, MK4])
points.append([MK4, MK5])
points.append([MK5, K  ])

points.append([K,   KG1])
points.append([KG1, KG2])
points.append([KG2, KG3])
points.append([KG3, KG4])
points.append([KG4, KG5])
points.append([KG5, Gamma])


path_num = int(sys.argv[12])
k_path = pb.results.make_path(*points[path_num], step=float(sys.argv[3]))
moments = int(sys.argv[2])
LIM = 10
nx = ny = 4
lx = ly = int(sys.argv[1])
weights = [float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]), 
           float(sys.argv[8]), float(sys.argv[9]), float(sys.argv[10]), float(sys.argv[11])]

configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=True, precision=1, spectrum_range=[-LIM,LIM])
calculation_arpes = kite.Calculation(configuration)
calculation_arpes.arpes(k_vector=k_path, num_moments=moments, weight=weights,num_disorder=1)
kite.config_system(lattice, configuration, calculation_arpes, filename='config.h5', disorder=disorder)
