""" Dirac Semimetal DOS (JPSPires) 
    script edited by SMJ

    Lattice : Double Node Dirac/Weyl Semimetal Model (Cubic Lattice);
    Disorder : Yes
    Configuration : size of the system 512x512x512, with domain decomposition (nx=ny=4,bz=1), periodic boundary conditions, double precision, fixed scaling)
    Calculation : DOS;
"""
import numpy as np
import matplotlib.pyplot as plt
import pybinding as pb
import kite
import sys

L1 = L2 = L3 = int(sys.argv[1]) #Linear Size of the Sample x-Direction
M = int(sys.argv[2])   #Number of Chebyshev Moments
C = float(sys.argv[3])
ths = int(sys.argv[4]) #Number of Physical Threads
scattering = float(sys.argv[5])
energy = float(sys.argv[6])
m = 0.0

def weyl():
    t = 1.0

    a1 = np.array([1, 0, 0])
    a2 = np.array([0, 1, 0])
    a3 = np.array([0, 0, 1])
    lat = pb.Lattice(a1=a1, a2=a2, a3=a3)
    
    lat.add_sublattices(
        ('A', [0, 0, 0], m),
        ('B', [0, 0, 0], -m)
    )

    lat.add_hoppings(
        ([1,0,0],'A','B',t),
        ([1,0,0],'B','A',t),

        ([0,1,0],'A','B', t*1j),
        ([0,1,0],'B','A',-t*1j),

        ([0,0,1],'A','A',t),
        ([0,0,1],'B','B',-t)
    )

    return lat

energy_scale = 4.0
norm_range = list(energy_scale*np.array([-1.0,1.0]))

lattice = weyl() # Build the Cubic Lattice with the Hoppings


# model = pb.Model(lattice, pb.translational_symmetry())
# solver = pb.solver.lapack(model)
# a = np.pi
# points = []
# points.append([0,0,0])
# points.append([a,0,0])
# points.append([a,a,a])
# points.append([-a,a/2,a])

# bands = solver.calc_bands(*points)
# bands.plot()
# plt.show()
# exit(1)


# Number of Subdomains (Thread Dependent Scheme)
if ths==1:
    nx,ny,nz = 1,1,1
elif ths==2:
    nx,ny,nz = 2,1,1
elif ths==4:
    nx,ny,nz = 2,2,1
elif ths==8:
    nx,ny,nz = 2,2,2
elif ths==16:
    nx,ny,nz = 4,2,2
elif ths==32:
    nx,ny,nz = 4,4,2
elif ths==64: nx,ny,nz = 4,4,4
else:
    sys.exit("INVALID NUMBER OF THREADS\n\nSCRIPT ABORTED!")
    
configuration = kite.Configuration(divisions=[nx, ny, nz], length=[L1, L2, L3], boundaries=[True, True, True], is_complex=True, precision=1,spectrum_range=norm_range, mss=1)

struc_disorder = kite.StructuralDisorder(lattice, concentration=C)
struc_disorder.add_vacancy('A')
struc_disorder.add_vacancy('B')

# Requires the calculation of DOS
calculation = kite.Calculation(configuration)
#calculation.dos(num_points=100, num_moments=M, num_random=1, num_disorder=1)
calculation.singleshot_conductivity_dc(energy=energy, num_moments=M, num_random=1, num_disorder=1, direction='xx', eta=scattering, preserve_disorder=True)

Center = [int(L1/2),int(L2/2),int(L3/2)]
Radius = 4
Potential = 0

kite.config_system(lattice, configuration, calculation, filename='config.h5', disorder_structural=struc_disorder, impurity=np.append(Center,[Radius,Potential]))
