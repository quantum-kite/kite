""" DoS for Benzene

    ##########################################################################
    #                         Copyright 2020/22, KITE                        #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in eV
    Lattice: Graphene Monolayer
    Disorder: Deterministic and Uniform at different sublattices
    Configuration : size of the system 512x512, without domain decomposition (nx=ny=1), periodic boundary conditions,
                    double precision, manual scaling
    Calculation type: DoS
    Last updated: 08/05/2025
"""
import kite
import pybinding as pb
from math import pi, sqrt
#from pybinding.repository import graphene

def naphthalene():
    """Return the lattice specification for naphthalene (2 rings) """
    N=2 # number of rings
    a = 0.2462  # [nm] site-site distance
    al = (N+1)*a*sqrt(3)  # [nm] unit cell length
    t = -1      # [eV] nearest neighbour hopping
    #t2 = 0.2
    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(
        a1=[al, 0],
        a2=[0, 3 * a]
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
        ('C10',[1.5 * a * sqrt(3), -0.5 * a], 0)
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
        ([0, 0], 'C9', 'C10', t)

        # between neighboring cells
#        ([1, -1], 'A', 'B', t),

    )

    return lat
# load graphene lattice and structural_disorder
lattice = naphthalene()
# add Disorder
delta=0.4
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
# - scaling, if None it's automatic, if present select spectrum_range=[e_min, e_max]
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=["periodic", "periodic"],
                                   is_complex=False, precision=1)
# require the calculation of DOS
calculation = kite.Calculation(configuration)
calculation.dos(num_points=4000, num_moments=512, num_random=5, num_disorder=1)
calculation.conductivity_optical(num_points=1000, num_disorder=1, num_random=5, num_moments=256, direction='xx')

# configure the *.h5 file
kite.config_system(lattice, configuration, calculation, filename='naphthalene.h5', disorder=disorder)
