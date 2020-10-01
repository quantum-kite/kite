import numpy as np
import pybinding as pb
import kite



W=0.1
energy_scale = 3.6 + W/2.0

def graphene_initial(onsite=(1.673819, -1.673819)):

    theta = np.pi / 3
    a1 = np.array([2 * np.sin(theta), 0])
    a2 = np.array([np.sin(theta), 1 + np.cos(theta)])

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(
    a1=a1,
    a2=a2
    )

    # Add sublattices
    lat.add_sublattices(
    # name, position, and onsite potential
    ('A', [0, 0], onsite[0]),
    ('B', [0, 1], onsite[1])
    )

    # Add hoppings
    lat.add_hoppings(
    # inside the main cell, between which atoms, and the value
    ([0, 0], 'A', 'B', - 1),
    # between neighboring cells, between which atoms, and the value
    ([0, -1], 'A', 'B', - 1),
    ([1, -1], 'A', 'B', - 1)
    )

    return lat



nx = ny = 2
N=32
B=0
Lx=256
Ly=256

lattice = graphene_initial()

disorder = kite.Disorder(lattice)
disorder.add_disorder('A', 'Uniform', 0.0, W*0.5/np.sqrt(3))
disorder.add_disorder('B', 'Uniform', 0.0, W*0.5/np.sqrt(3))
configuration = kite.Configuration(divisions=[nx, ny], length=[Lx, Ly], boundaries=[True, True],
                                      is_complex=(B!=0), precision=0, spectrum_range=[-energy_scale, energy_scale])
calculation = kite.Calculation(configuration)
calculation.conductivity_optical_nonlinear(num_points=1000, num_moments=N, num_random=1, num_disorder=1,
                                                  direction='yyy', temperature=0.01)

kite.config_system(lattice, configuration, calculation, disorder=disorder, filename='config.h5')

