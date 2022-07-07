"""       
        ##############################################################################      
        #                        KITE | Release  1.1                                 #      
        #                                                                            #      
        #                        Kite home: quantum-kite.com                         #           
        #                                                                            #      
        #  Developed by: Simao M. Joao, Joao V. Lopes, Tatiana G. Rappoport,         #       
        #  Misa Andelkovic, Lucian Covaci, Aires Ferreira, 2018-2022                 #      
        #                                                                            #      
        ##############################################################################      
"""
""" Twisted bilayer graphene, loading a predefined lattice

    Lattice :  Twisted bilayer graphene;
    Disorder : None;
    Configuration : size of the system flexible, with domain decomposition, periodic boundary conditions,
                    double precision, automatic scaling;
    Calculation : dos;
   

"""

import pybinding as pb
import matplotlib.pyplot as plt
import numpy as np
import kite


# define the angle
angle = 21.787
# define the name of the pb.Lattice object
name = 'lattice_twisted_bilayer/lattice_tblg_{:.3f}'.format(angle)
# number of decomposition parts in each direction of matrix. This divides the lattice into various sections,
# each of which is calculated in parallel
nx = 1
ny = 1
# number of unit cells in each direction.
l1 = l2 = 64 * nx, 64 * nx
# estimated scale for the Hamiltonian
energy_scale = 10
# load a predefined lattice object, the lattice can be saved with pb.save(lattice, name)
lattice = pb.load(name)
# make config object which caries info about
# - the number of decomposition parts [nx, ny],
# - lengths of structure [lx, ly]
# - boundary conditions, setting True as periodic boundary conditions, and False elsewise,
# - info if the exported hopping and onsite data should be complex,
# - info of the precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
# - scaling, if None it's automatic, if present select spectrum_bound=[e_min, e_max]
configuration = kite.Configuration(divisions=[nx, ny], length=[l1, l2], boundaries=[True, True],
                                   is_complex=False, precision=1, spectrum_range=[-10, 10])
# require the calculation of DOS
num_moments = 1000
calculation = kite.Calculation(configuration)
calculation.dos(num_moments=num_moments, num_random=1, num_disorder=1, num_points=1000)
# make modification object which caries info about

# for a quick check, let's make a Pybinding model and check the DOS
model = pb.Model(lattice, pb.rectangle(100, 100),
                 pb.translational_symmetry(a1=50, a2=50))
# if you would like to specify Disorder, use other function that takes of converting KITE to Pybinding disorder objects
# model = kite.make_pybinding_model(lattice)
# dos = pb.kpm(model).calc_dos(np.linspace(-4, 4, 2000), broadening=1e-2, num_random=1)
# dos.plot()
# plt.show()

# configure the *.h5 file
kite.config_system(lattice, configuration, calculation, filename='tblg_{:.3f}.h5'.format(angle))
