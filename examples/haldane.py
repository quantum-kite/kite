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
""" DOS and DC conductivity of the Haldane model

    Lattice : Honeycomb lattice;
    Disorder : Disorder class Uniform at different sublattices;
    Configuration : size of the system 256x256, with domain decomposition (nx=ny=1), periodic boundary conditions,
                    double precision, automatic scaling;
    Calculation : DOS and conductivity_dc (xy);
"""

import kite
import pybinding as pb
from math import sqrt


def haldane():
    """Return the lattice specification for Haldane model"""

    a = 0.24595   # [nm] unit cell length
    a_cc = 0.142  # [nm] carbon-carbon distance

    t=-1
    t2 = t/10
    m=0
    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(
        a1=[a, 0],
        a2=[a/2, a/2 * sqrt(3)]
    )

    lat.add_sublattices(
        # name and position
        ('A', [0, -a_cc/2],-m),
        ('B', [0,  a_cc/2],m)
    )

    lat.add_hoppings(
        # inside the main cell
        ([0,  0], 'A', 'B', t),
        # between neighboring cells
        ([1, -1], 'A', 'B', t),
        ([0, -1], 'A', 'B', t),
        ([1, 0], 'A', 'A', t2 * 1j),
        ([0, -1], 'A', 'A', t2 * 1j),
        ([-1, 1], 'A', 'A', t2 * 1j),
        ([1, 0], 'B', 'B', t2 * -1j),
        ([0, -1], 'B', 'B', t2 * -1j),
        ([-1, 1], 'B', 'B', t2 * -1j)
    )

    return lat


# make a haldane lattice
lattice = haldane()
# add Uniformly distributed disorder
disorder = kite.Disorder(lattice)
disorder.add_disorder('A', 'Uniform', +0.0, 0.4)
disorder.add_disorder('B', 'Uniform', +0.0, 0.4)
# number of decomposition parts in each direction of matrix.
# This divides the lattice into various sections, each of which is calculated in parallel
nx = ny = 1
# number of unit cells in each direction.
lx = ly = 256
# make config object which caries info about
# - the number of decomposition parts [nx, ny],
# - lengths of structure [lx, ly]
# - boundary conditions [mode,mode, ... ] with modes:
#   . "periodic"
#   . "open"
#   . "twist_fixed"     this option needs the extra argument ths=[phi_1,..,phi_DIM] where phi_i \in [0, 2*M_PI]  
#   . "twist_random"
#
# - info if the exported hopping and onsite data should be complex,
# - info of the precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
# - scaling, if None it's automatic, if present select spectrum_range=[e_min, e_max]
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=["periodic", "periodic"],
                                   is_complex=True, precision=0)
calculation = kite.Calculation(configuration)
# require the calculation of DOS and conductivity_dc
calculation.dos(num_points=1000, num_moments=512, num_random=10, num_disorder=1)
calculation.conductivity_dc(num_points=1000, num_moments=256, num_random=50, num_disorder=1,
                            direction='xy', temperature=100)
# configure the *.h5 file
kite.config_system(lattice, configuration, calculation, filename='haldane.h5',
                   disorder=disorder)
# ATTENTION: to generate  the conductivity data file for a desired window of Fermi energies, please use 
#./KITE-tools h5_file.h --CondDC -F Emin Emax NumPoints 
#See ./KITE-tools --help for more options
