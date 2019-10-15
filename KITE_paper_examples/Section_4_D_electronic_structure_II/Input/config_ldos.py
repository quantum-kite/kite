import pybinding as pb
import sys
import kite
import numpy as np
# import matplotlib.pyplot as plt
from math import pi, sqrt
"""Calculate the band structure of Haldane model with sublattice symmetry breaking"""
# the gaps have different sizes in K and K'!

def tmd(signsoc):
    """Return the lattice specification for monolayer graphene"""
    a = 0.32   # [nm] unit cell length
    a_cc = a/sqrt(3)  # [nm] site-site distance
    t = -1      # [eV] nearest neighbour hopping
    t2 = 0.2  # nest-nearest neighbour hopping
    m=0.1  # sub-lattice symmetry breaking term
    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(
        a1=[a, 0],
        a2=[a/2, a/2 * sqrt(3)]
    )
    lambdaM =signsoc*0.251
    lambdaX =signsoc*0.057

    Delta0=-0.872
    Delta1=0.067
    Delta2=-1.511
    Deltap=-3.468
    Deltaz=-3.913

    Vpdsigma=3.603
    Vpdpi=-0.942

    Vddsigma = -1.216
    Vddpi = 0.177
    Vdddelta = 0.243

    Vppsigma=0.749
    Vpppi=0.236




    txx =Vpppi
    tyy =Vpppi
    tzz =Vppsigma
    t2m3x1=sqrt(2)/(7*sqrt(7))*14*Vpdpi
    t2m1x2=sqrt(2)/(7*sqrt(7))*(-6*sqrt(3)*Vpdpi+2*Vpdsigma)
    t2m1x3 = sqrt(2) / (7 * sqrt(7)) * (12 * Vpdpi + sqrt(3) * Vpdsigma)
    t2m2x2 = sqrt(2) / (7 * sqrt(7)) * (-6 * Vpdpi - 4 * sqrt(3)  * Vpdsigma)
    t2m2x3 = sqrt(2) / (7 * sqrt(7)) * (4 * sqrt(3) * Vpdpi - 6 * Vpdsigma)

    t1m1x1 = sqrt(2) / (7 * sqrt(7)) * (-9 * Vpdpi + sqrt(3) * Vpdsigma)
    t1m1x2 = sqrt(2) / (7 * sqrt(7)) * (3 *  sqrt(3) * Vpdpi - 1 * Vpdsigma)
    t1m1x3 = sqrt(2) / (7 * sqrt(7)) * (12 *  Vpdpi +  sqrt(3) * Vpdsigma)

    t1m2x1 = sqrt(2) / (7 * sqrt(7)) * (5 * sqrt(3) * Vpdpi + 3 * Vpdsigma)
    t1m2x2 = sqrt(2) / (7 * sqrt(7)) * (9 * Vpdpi - sqrt(3) * Vpdsigma)
    t1m2x3 = sqrt(2) / (7 * sqrt(7)) * (-2 * sqrt(3) * Vpdpi + 3 * Vpdsigma)


    t1m3x1 = sqrt(2) / (7 * sqrt(7)) * (-1 * Vpdpi -3 * sqrt(3) * Vpdsigma)
    t1m3x2 = sqrt(2) / (7 * sqrt(7)) * (5 * sqrt(3) * Vpdpi + 3 * Vpdsigma)
    t1m3x3 = sqrt(2) / (7 * sqrt(7)) * (6 * Vpdpi -3 * sqrt(3) * Vpdsigma)

    t1m1m1=0.25*(3*Vdddelta+Vddsigma)
    t1m1m2=0.25*0.5*sqrt(3)*(-Vdddelta+Vddsigma)
    t1m1m3=0.25*(-1.5)*(Vdddelta-Vddsigma)
    t1m2m2=0.25*0.25*(Vdddelta+12*Vddpi+3*Vddsigma)
    t1m2m3=0.25*0.25*sqrt(3)*(Vdddelta-4*Vddpi+3*Vddsigma)
    t1m3m3=0.25*0.25*(3*Vdddelta+4*Vddpi+9*Vddsigma)

    t2m1m1 = 0.25 * (3 * Vdddelta + Vddsigma)
    t2m1m2 = 0.25 * sqrt(3)*(Vdddelta - Vddsigma)
    t2m2m2 = 0.25 * ( Vdddelta + 3*Vddsigma)
    t2m3m3 = 0.25 * (4 * Vddpi)

    t1x1x1 = 0.25 * (3 * Vpppi + Vppsigma)
    t1x1x2 = 0.25 * sqrt(3) * (Vpppi - Vppsigma)
    t1x2x2 = 0.25 * (Vpppi + 3 * Vppsigma)
    t1x3x3 = 0.25 * (4 * Vpppi)

    t2x1x1=Vppsigma
    t2x2x2=Vpppi
    t2x3x3=Vpppi


    lat.add_sublattices(
        # name and position
        ('M1', [0, -a_cc/2], Delta0),
        ('M2', [0, -a_cc/2], Delta2),
        ('M3', [0, -a_cc/2], Delta2),
        ('X1', [0,  a_cc/2], Deltap+txx),
        ('X2', [0,  a_cc/2], Deltap+tyy),
        ('X3', [0,  a_cc/2], Deltaz-tzz)
    )

    lat.add_hoppings(
        # inside the main cell
        #([0, 0], 'M1', 'M2', t),
        #([0, 0], 'M1', 'M3', t),
        ([0, 0], 'M2', 'M3', -1j*lambdaM),
        #([0, 0], 'M3', 'M2', 1j * lambdaM),
        ([0, 0], 'X1', 'X2', -1j*lambdaX),
        #([0, 0], 'X2', 'X1', 1j * lambdaX),

        #([0, 0], 'X1', 'X3', t),
        #([0, 0], 'X2', 'X3', t),

         ([0, 0], 'M1', 'X2', t2m1x2),
         ([0, 0], 'M1', 'X3', t2m1x3),
         ([0, 0], 'M2', 'X2', t2m2x2),
         ([0, 0], 'M2', 'X3', t2m2x3),
         ([0, 0], 'M3', 'X1', t2m3x1),
        #
        # # between neighboring cells
         ([1, -1], 'M1', 'X1', t1m1x1),
         ([1, -1], 'M1', 'X2', t1m1x2),
         ([1, -1], 'M1', 'X3', t1m1x3),
         ([1, -1], 'M2', 'X1', t1m2x1),
         ([1, -1], 'M2', 'X2', t1m2x2),
         ([1, -1], 'M2', 'X3', t1m2x3),
         ([1, -1], 'M3', 'X1', t1m3x1),
         ([1, -1], 'M3', 'X2', t1m3x2),
         ([1, -1], 'M3', 'X3', t1m3x3),

         ([0, -1], 'M1', 'X1', -t1m1x1),
         ([0, -1], 'M1', 'X2', t1m1x2),
         ([0, -1], 'M1', 'X3', t1m1x3),
         ([0, -1], 'M2', 'X1', -t1m2x1),
         ([0, -1], 'M2', 'X2', t1m2x2),
         ([0, -1], 'M2', 'X3', t1m2x3),
         ([0, -1], 'M3', 'X1', t1m3x1),
         ([0, -1], 'M3', 'X2', -t1m3x2),
         ([0, -1], 'M3', 'X3', -t1m3x3),

        ([1, 0], 'M1', 'M1', t2m1m1),
        ([1, 0], 'M1', 'M2', t2m1m2),
        ([1, 0], 'M2', 'M1', t2m1m2),
        ([1, 0], 'M2', 'M2', t2m2m2),
        ([1, 0], 'M3', 'M3', t2m3m3),

        ([-1, 1], 'M1', 'M1', t1m1m1),
        ([-1, 1], 'M1', 'M2', t1m1m2),
        ([-1, 1], 'M1', 'M3', t1m1m3),
        ([-1, 1], 'M2', 'M1', t1m1m2),
        ([-1, 1], 'M2', 'M2', t1m2m2),
        ([-1, 1], 'M2', 'M3', t1m2m3),
        ([-1, 1], 'M3', 'M1', t1m1m3),
        ([-1, 1], 'M3', 'M2', t1m2m3),
        ([-1, 1], 'M3', 'M3', t1m3m3),

        ([0, -1], 'M1', 'M1', t1m1m1),
        ([0, -1], 'M1', 'M2', t1m1m2),
        ([0, -1], 'M1', 'M3', -t1m1m3),
        ([0, -1], 'M2', 'M1', t1m1m2),
        ([0, -1], 'M2', 'M2', t1m2m2),
        ([0, -1], 'M2', 'M3', -t1m2m3),
        ([0, -1], 'M3', 'M1', -t1m1m3),
        ([0, -1], 'M3', 'M2', -t1m2m3),
        ([0, -1], 'M3', 'M3', t1m3m3),

         ([1, 0], 'X1', 'X1', t2x1x1),
         ([1, 0], 'X2', 'X2', t2x2x2),
         ([1, 0], 'X3', 'X3', t2x3x3),

         ([-1, 1], 'X1', 'X1', t1x1x1),
         ([-1, 1], 'X1', 'X2', t1x1x2),
         ([-1, 1], 'X2', 'X1', t1x1x2),
         ([-1, 1], 'X2', 'X2', t1x2x2),
         ([-1, 1], 'X3', 'X3', t1x1x1),

         ([0, -1], 'X1', 'X1', t1x1x1),
         ([0, -1], 'X1', 'X2', -t1x1x2),
         ([0, -1], 'X2', 'X1', -t1x1x2),
         ([0, -1], 'X2', 'X2', t1x2x2),
         ([0, -1], 'X3', 'X3', t1x1x1)

    )

    return lat

lattice = tmd(1.0)
delta=0.1


N = int(sys.argv[2])
lim = 10
d1 = d2 = 32

# List of points where to calculate the local density of states
pos_matrix=[]
sub_matrix=[]
for i in range(-lim,lim):
  for j in range(-lim,lim):
    pos_matrix.append([d1+i,d2+j])
    pos_matrix.append([d1+i,d2+j])
    pos_matrix.append([d1+i,d2+j])
    sub_matrix.append('M1')
    sub_matrix.append('M2')
    sub_matrix.append('M3')

# Insert a vacancy in position 32,32
struc_disorder_1 = kite.StructuralDisorder(lattice, position=[[d1,d2]])
struc_disorder_2 = kite.StructuralDisorder(lattice, position=[[d1,d2]])
struc_disorder_3 = kite.StructuralDisorder(lattice, position=[[d1,d2]])
struc_disorder_1.add_vacancy('M1')
struc_disorder_2.add_vacancy('M2')
struc_disorder_3.add_vacancy('M3')
disorder_structural = [struc_disorder_1, struc_disorder_2, struc_disorder_3]  


nx = ny = 4
lx = ly = int(sys.argv[1])
configuration = kite.Configuration(divisions=[nx, ny], length=[lx, ly], boundaries=[True, True], is_complex=True, precision=1, spectrum_range = [-10,10])
calculation = kite.Calculation(configuration)
calculation.ldos(energy=np.linspace(-0.1, 0.1, 3), num_moments=N, num_disorder=1, position=pos_matrix, sublattice=sub_matrix)
calculation.dos(num_points=512, num_moments=512, num_random=1, num_disorder=1)
kite.config_system(lattice, configuration, calculation, filename='config.h5', disorder_structural=disorder_structural)

