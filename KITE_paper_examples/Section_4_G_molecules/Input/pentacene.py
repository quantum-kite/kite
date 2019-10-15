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
""" Onsite disorder

    Lattice : Monolayer graphene;
    Disorder : Disorder class Deterministic and Uniform at different sublattices,
    Configuration : size of the system 512x512, without domain decomposition (nx=ny=1), periodic boundary conditions,
                    double precision, manual scaling;
    Calculation : dos;
"""

import kite
import pybinding as pb
from math import pi, sqrt
#from pybinding.repository import graphene



def pentacene():
    """Return the lattice specification for tetracene (4 rings) """
    N=5 # number of rings
    a = 0.2462  # [nm] site-site distance
    al = (N+1)*a*sqrt(3)  # [nm] unit cell length
    t = -1      # [eV] nearest neighbour hopping
    #t2 = 0.2
    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(
        a1=[al, 0],
        a2=[0, 3*a]
    )

    lat.add_sublattices(
        # name and position
        ('C1', [-0.5*a*sqrt(3),0.5*a],0),
        ('C2', [-0.5*a*sqrt(3),-0.5*a],0),
        ('C3', [0, a], 0),
        ('C4', [0, -a], 0),
        ('C5', [0.5 * a * sqrt(3), 0.5 * a], 0),
        ('C6', [0.5 * a * sqrt(3), -0.5 * a], 0),
        ('C7', [a * sqrt(3), a], 0),
        ('C8', [a * sqrt(3), -a], 0),
        ('C9', [1.5 * a * sqrt(3), 0.5 * a], 0),
        ('C10',[1.5 * a * sqrt(3), -0.5 * a], 0),
        ('C11', [2*a * sqrt(3), a], 0),
        ('C12', [2*a * sqrt(3), -a], 0),
        ('C13', [2.5 * a * sqrt(3), 0.5 * a], 0),
        ('C14', [2.5 * a * sqrt(3), -0.5 * a], 0),
        ('C15', [3 * a * sqrt(3), a], 0),
        ('C16', [3 * a * sqrt(3), -a], 0),
        ('C17', [3.5 * a * sqrt(3), 0.5 * a], 0),
        ('C18', [3.5 * a * sqrt(3), -0.5 * a], 0),
        ('C19', [4 * a * sqrt(3), a], 0),
        ('C20', [4 * a * sqrt(3), -a], 0),
        ('C21', [4.5 * a * sqrt(3), 0.5 * a], 0),
        ('C22', [4.5 * a * sqrt(3), -0.5 * a], 0)
    )

    lat.add_hoppings(
        # inside the main cell
        ([0,  0], 'C1', 'C2', t),
        ([0,  0], 'C1', 'C3', t),
        ([0,  0], 'C3', 'C5', t),
        ([0, 0], 'C2', 'C4', t),
        ([0, 0], 'C4', 'C6', t),
        ([0, 0], 'C5', 'C6', t),
        ([0,  0], 'C5', 'C7', t),
        ([0,  0], 'C7', 'C9', t),
        ([0,  0], 'C6', 'C8', t),
        ([0, 0], 'C8', 'C10', t),
        ([0, 0], 'C9', 'C10', t),
        ([0, 0], 'C9', 'C11', t),
        ([0, 0], 'C11', 'C13', t),
        ([0, 0], 'C10', 'C12', t),
        ([0, 0], 'C12', 'C14', t),
        ([0, 0], 'C13', 'C14', t),
        ([0, 0], 'C13', 'C15', t),
        ([0, 0], 'C15', 'C17', t),
        ([0, 0], 'C14', 'C16', t),
        ([0, 0], 'C16', 'C18', t),
        ([0, 0], 'C17', 'C18', t),
        ([0, 0], 'C17', 'C19', t),
        ([0, 0], 'C19', 'C21', t),
        ([0, 0], 'C18', 'C20', t),
        ([0, 0], 'C20', 'C22', t),
        ([0, 0], 'C21', 'C22', t)

        # between neighboring cells
#        ([1, -1], 'A', 'B', t),

    )

    return lat
# load graphene lattice and structural_disorder
lattice = pentacene()
# add Disorder
delta=0.25
disorder = kite.Disorder(lattice)
#disorder.add_disorder('B', 'Deterministic', -1.0)
disorder.add_disorder('C1', 'Uniform', 0, delta)
disorder.add_disorder('C2', 'Uniform', 0, delta)
disorder.add_disorder('C3', 'Uniform', 0, delta)
disorder.add_disorder('C4', 'Uniform', 0, delta)
disorder.add_disorder('C5', 'Uniform', 0, delta)
disorder.add_disorder('C6', 'Uniform', 0, delta)
disorder.add_disorder('C7', 'Uniform', 0, delta)
disorder.add_disorder('C8', 'Uniform', 0, delta)
disorder.add_disorder('C9', 'Uniform', 0, delta)
disorder.add_disorder('C10', 'Uniform', 0, delta)
disorder.add_disorder('C11', 'Uniform', 0, delta)
disorder.add_disorder('C12', 'Uniform', 0, delta)
disorder.add_disorder('C13', 'Uniform', 0, delta)
disorder.add_disorder('C14', 'Uniform', 0, delta)
disorder.add_disorder('C16', 'Uniform', 0, delta)
disorder.add_disorder('C17', 'Uniform', 0, delta)
disorder.add_disorder('C18', 'Uniform', 0, delta)
disorder.add_disorder('C19', 'Uniform', 0, delta)
disorder.add_disorder('C20', 'Uniform', 0, delta)
disorder.add_disorder('C21', 'Uniform', 0, delta)

# number of decomposition parts in each direction of matrix.
# This divides the lattice into various sections, each of which is calculated in parallel
nx = ny = 1
# number of unit cells in each direction.
lx = ly = 128
# make config object which caries info about
# - the number of decomposition parts [nx, ny],
# - lengths of structure [lx, ly]
# - boundary conditions, setting True as periodic boundary conditions, and False elsewise,
# - info if the exported hopping and onsite data should be complex,
# - info of the precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
# - scaling, if None it's automatic, if present select spectrum_bound=[e_min, e_max]
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True],
                                   is_complex=False, precision=1)
# require the calculation of DOS
calculation = kite.Calculation(configuration)
calculation.dos(num_points=4000, num_moments=512, num_random=5, num_disorder=1)
calculation.conductivity_optical(num_points=1000, num_disorder=1, num_random=5, num_moments=256, direction='xx')
# configure the *.h5 file
kite.config_system(lattice, configuration, calculation, filename='pentacene.h5', disorder=disorder)
