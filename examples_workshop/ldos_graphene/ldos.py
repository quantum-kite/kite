import matplotlib.pyplot as plt
import kite
import numpy as np
import pybinding as pb

timp = -0.0
limp = -1.0
EnergyScale = 4.0

def graphene_initial(onsite=(0, 0)):
    theta = np.pi / 3
    a1 = np.array([1 + np.cos(theta), np.sin(theta)])
    a2 = np.array([0, 2 * np.sin(theta)])

    lat = pb.Lattice(
        a1=a1,
        a2=a2
    )
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, 0], onsite[0]),
        ('B', [1, 0], onsite[1])
    )

    lat.add_hoppings(
        # inside the main cell, between which atoms, and the value
        ([0, 0], 'A', 'B', - 1),
        # between neighboring cells, between which atoms, and the value
        ([-1, 0], 'A', 'B', - 1/EnergyScale),
        ([-1, 1], 'A', 'B', - 1/EnergyScale)
    )


    node0 = [[+0,+0],'A']
    node1 = [[+0,+0],'B']
    node2 = [[-1,+0],'B']
    node3 = [[-1,+1],'B']
    
    struc_disorder_one = kite.StructuralDisorder(lat, position=[[32, 32]])
    struc_disorder_one.add_structural_disorder(
        (*node0, *node1, timp),
        (*node0, *node2, timp),
        (*node0, *node3, timp),
        (*node0, limp)
    )
    
    return lat, [struc_disorder_one]

lattice, disorded_structural = graphene_initial()

nx = ny = 1
lx = 128
ly = 128
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                 is_complex=True, precision=1, spectrum_range=[-EnergyScale, EnergyScale])


pos_matrix=[]
sub_matrix=[]
d1 = 32
d2 = 32
N = 64
for i in range(-5,5):
  for j in range(-5,5):
    pos_matrix.append([d1+i,d2+j])
    pos_matrix.append([d1+i,d2+j])
    sub_matrix.append('A')
    sub_matrix.append('B')

calculation = kite.Calculation(configuration)
calculation.ldos(energy=np.linspace(-1, 1, 100), num_moments=N, num_disorder=1,
                     position=pos_matrix, sublattice=sub_matrix)

kite.config_system(lattice, configuration, calculation, filename='ldos.h5', disorder_structural=disorded_structural)
