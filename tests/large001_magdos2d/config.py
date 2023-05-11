"""       
        ##############################################################################      
        #                        KITE | Release  1.0                                 #      
        #                                                                            #      
        #                        Kite home: quantum-kite.com                         #           
        #                                                                            #      
        #  Developed by: Simao M. Joao, Joao V. Lopes, Tatiana G. Rappoport,         #       
        #  Misa Andelkovic, Lucian Covaci, Aires Ferreira, 2018-2021                 #      
        #                                                                            #      
        ##############################################################################      
"""
import kite
import numpy as np
import pybinding as pb



def square_lattice(onsite=(0, 0)):
    a1 = np.array([1, 0])
    a2 = np.array([0, 1])

    lat = pb.Lattice(
        a1=a1, a2=a2
    )

    lat.add_sublattices(
        ('A', [0, 0], onsite[0])
    )

    lat.add_hoppings(
        ([1, 0], 'A', 'A', - 1),
        ([0, 1], 'A', 'A', - 1)
    )

    return lat

lattice = square_lattice()
nx = ny = 2
lx = ly = 512
N = 4192
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=["periodic", "periodic"], is_complex=True, precision=1, spectrum_range=[-4.1,4.1])
calculation = kite.Calculation(configuration)
mod = kite.Modification(magnetic_field = 9)
calculation.dos(num_points=4096, num_moments=N, num_random=1, num_disorder=1)

kite.config_system(lattice, configuration, calculation, modification=mod, filename='config.h5')
