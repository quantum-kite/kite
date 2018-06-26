"""       
        ##############################################################################      
        #                        KITE | Pre-Release version 0.1                      #      
        #                                                                            #      
        #                        Kite home: quantum-kite.com                         #           
        #                                                                            #      
        #  Developed by: Simao M. Joao, Joao V. Lopes, Tatiana G. Rappoport,         #       
        #  Misa Andelkovic, Lucian Covaci, Aires Ferreira, 2018                      #      
        #                                                                            #      
        ##############################################################################      
"""
""" Optical conductivity of disordered graphene lattice

    Lattice : Monolayer graphene;
    Disorder : Disorder class Gaussian at different sublattices;
    Configuration : size of the system 512x512, without domain decomposition (nx=ny=1), periodic boundary conditions,
                    double precision, manual scaling;
    Calculation : dos; optical conductivity (xx)

"""

import kite

from pybinding.repository import graphene


# load a monolayer graphene lattice
lattice = graphene.monolayer()
# add Disorder
disorder = kite.Disorder(lattice)
disorder.add_disorder('B', 'Gaussian', 0.0, 0.6)
disorder.add_disorder('A', 'Gaussian', 0.0, 0.6)
# number of decomposition parts in each direction of matrix. This divides the lattice into various sections,
# each of which is calculated in parallel
nx = ny = 1
# number of unit cells in each direction.
lx = ly = 512
# make config object which caries info about
# - the number of decomposition parts [nx, ny],
# - lengths of structure [lx, ly]
# - boundary conditions, setting True as periodic boundary conditions, and False elsewise,
# - info if the exported hopping and onsite data should be complex,
# - info of the precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
# - scaling, if None it's automatic, if present select spectrum_bound=[e_min, e_max]
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=False, precision=1)
# require the calculation of DOS and conductivity_optical
calculation = kite.Calculation(configuration)
calculation.dos(num_points=4000, num_moments=512, num_random=10, num_disorder=1)
calculation.conductivity_optical(num_points=1000, num_disorder=1, num_random=5, num_moments=256, direction='xx')
# configure the *.h5 file
kite.config_system(lattice, configuration, calculation, filename='optcond_gaussian_disorder.h5', disorder=disorder)
