import matplotlib.pyplot as plt
import kite
import numpy as np
import pybinding as pb


# define lattice of monolayer graphene with 1[nm] interatomic distance and t=1/3[eV] hopping
#  INFO: other examples are defined in define_lattice.py script

rho  = 0.0000+0.0000j
timp = 0.1500+0.0000j
nu   = 0.0000+0.1000j
EnergyScale = (3  +  6*np.abs(np.real(timp))) * 1.01

def graphene_initial(onsite=(0, 0)):
    """Return the basic lattice specification for monolayer graphene with nearest neighbor"""

    theta = np.pi / 3
    a1 = np.array([1 + np.cos(theta), np.sin(theta)])
    a2 = np.array([0, 2 * np.sin(theta)])

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(
        a1=a1,
        a2=a2
    )

    # Add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, 0], onsite[0]),
        ('B', [1, 0], onsite[1])
    )

    # Add hoppings
    lat.add_hoppings(
        # inside the main cell, between which atoms, and the value
        ([0, 0], 'A', 'B', - 1/EnergyScale),
        # between neighboring cells, between which atoms, and the value
        ([-1, 0], 'A', 'B', - 1/EnergyScale),
        ([-1, 1], 'A', 'B', - 1/EnergyScale)
    )

    # Add disorder
    # Each sublattice can have different disorder. If there are multiple orbitals at one sublattice, one needs to add
    # disorder vector of the same size as the number of orbitals. Type of disorder available are Gaussian,
    # Deterministic and Uniform. Each of the needs the have mean value, and standard deviation, where standard deviation
    # of deterministic disorder should be 0.

    disorder = kite.Disorder(lat)
    disorder.add_disorder('A', 'Deterministic', 0.0, 0)
    disorder.add_disorder('B', 'Deterministic', 0.0, 0)


    # Same procedure as adding the hopping + concentration
    node0 = [[+0,+0],'A']
    node1 = [[+0,+0],'B']
    node2 = [[+1,+0],'A']
    node3 = [[+0,+1],'B']
    node4 = [[+0,+1],'A']
    node5 = [[-1,+1],'B']
    
    struc_disorder_one = kite.StructuralDisorder(lat, position=[[64,64], [64,32], [32, 32]])
    struc_disorder_one.add_structural_disorder(
        (*node0, *node1, timp/EnergyScale),
        (*node1, *node2, timp/EnergyScale),
        (*node2, *node3, timp/EnergyScale),
        (*node3, *node4, timp/EnergyScale),
        (*node4, *node5, timp/EnergyScale),
        (*node5, *node0, timp/EnergyScale),
        #
        (*node0, *node2, nu/EnergyScale),
        (*node2, *node4, nu/EnergyScale),
        (*node4, *node0, nu/EnergyScale),
        (*node1, *node3, nu/EnergyScale),
        (*node3, *node5, nu/EnergyScale),
        (*node5, *node1, nu/EnergyScale),
        #
        (*node0, *node3, rho/EnergyScale),
        (*node1, *node4, rho/EnergyScale),
        (*node2, *node5, rho/EnergyScale),
        # in this way we can add onsite disorder again        
        (*node0, 0.)
    )
    
    # if there is disorder it should be returned separately from the lattice
    return lat, disorder, [struc_disorder_one]

lattice, disorder, disorded_structural = graphene_initial()
# number of decomposition parts in each direction of matrix.

nx = ny = 4
# number of unit cells in each direction.
lx = 256
ly = 256
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                 is_complex=True, precision=1, spectrum_range=[-10,10])


pos_matrix=[]
sub_matrix=[]
d1 = 32
d2 = 32
for i in range(-5,5):
  for j in range(-5,5):
    pos_matrix.append([d1+i,d2+j])
    pos_matrix.append([d1+i,d2+j])
    sub_matrix.append('A')
    sub_matrix.append('B')

calculation = kite.Calculation(configuration)
calculation.ldos(energy=np.linspace(-1, 1, 4), num_moments=256, num_disorder=1,
                     position=pos_matrix, sublattice=sub_matrix)

kite.config_system(lattice, configuration, calculation, filename='config.h5', disorder_structural=disorded_structural)
