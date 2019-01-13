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
""" Graphene DOS

    Lattice : Graphene lattice;
    Disorder : None;
    Configuration : size of the system 512x512, without domain decomposition (nx=ny=1), periodic boundary conditions,
                    double precision, automatic scaling;
    Calculation : DOS;
"""
import numpy as np
import pybinding as pb
import kite


def graphene(onsite=(0, 0)):
    """Make a honeycomb lattice with nearest neighbor hopping

    Parameters
    ----------
    onsite : tuple or list
        Onsite energy at different sublattices.
    """

    theta = np.pi / 3
    t = 2.8  # eV
    a1 = np.array([1, 0, 0] )
    a2 = np.array([0, 1, 0])
    a3 = np.array([0, 0, 1])
    lat = pb.Lattice(
        a1=a1, a2=a2, a3=a3
    )
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, 0, 0], onsite[0])
    )
    lat.add_hoppings(
        ([1, 0, 0], 'A', 'A', - t),
        ([0, 1, 0], 'A', 'A', - t),
        ([0, 0, 1], 'A', 'A', - t)
    )

    return lat


# make a graphene lattice
lattice = graphene()
disorder = kite.Disorder(lattice)

#disorder.add_disorder('A', 'Uniform', +0.0, 0.5)
#disorder_struc = kite.StructuralDisorder(lattice, position=[1, 1])
#disorder_struc.add_vacancy('A')


# number of decomposition parts in each direction of matrix.
# This divides the lattice into various sections, each of which is calculated in parallel
nx = 1
ny = 1
nz = 2
# number of unit cells in each direction.
lx = ly = lz = 128
# make config object which caries info about
# - the number of decomposition parts [nx, ny],
# - lengths of structure [lx, ly]
# - boundary conditions, setting True as periodic boundary conditions, and False elsewise,
# - info if the exported hopping and onsite data should be complex,
# - info of the precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
# - scaling, if None it's automatic, if present select spectrum_bound=[e_min, e_max]
configuration = kite.Configuration(divisions=[nx, ny, nz], length=[lx, ly, lz], boundaries=[True, True, True],
                                   is_complex=False, precision=1)
# require the calculation of DOS
calculation = kite.Calculation(configuration)
calculation.dos(num_points=1000, num_moments=128, num_random=1, num_disorder=1)
# configure the *.h5 file
kite.config_system(lattice, configuration, calculation, filename='example_initial.h5')
