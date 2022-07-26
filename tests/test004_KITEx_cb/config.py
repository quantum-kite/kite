import numpy as np
import pybinding as pb
import kite

def cubic():
    a1 = np.array([1,0,0]) 
    a2 = np.array([0,1,0]) 
    a3 = np.array([0,0,1]) 
    lat = pb.Lattice( a1=a1, a2=a2, a3=a3)

    lat.add_sublattices( ('A', [0, 0], 0))

    lat.add_hoppings(
        ([0,0,1], 'A', 'A', 1),
        ([0,1,0], 'A', 'A', 1),
        ([1,0,0], 'A', 'A', 1)
    )
    return lat

lattice = cubic()
nx = ny = nz = 1
lx = ly = lz = 32
mode = "periodic"
configuration = kite.Configuration(divisions=[nx,ny,nz], length=[lx,ly,lz], boundaries=[mode, mode, mode], is_complex=False, precision=1, spectrum_range=[-7,7])
calculation = kite.Calculation(configuration)
calculation.dos(num_points=1000, num_moments=64, num_random=1, num_disorder=1)
kite.config_system(lattice, configuration, calculation, filename='config.h5')
