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
""" Bilayer Graphene with Rashba SOC """

import pybinding as pb
import kite
import sys
import numpy as np
from os import system as terminal

def cube(onsite=(0, 0)):
    """Return lattice specification for a cube lattice with nearest neighbor hoppings"""

    # parameters
    t = 1

    # define lattice vectors
    a1 = np.array([1, 0, 0])
    a2 = np.array([0, 1, 0])
    a3 = np.array([0, 0, 1])

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(a1=a1, a2=a2, a3=a3)

    # add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [0, 0, 0], onsite[0])
    )

    # Add hoppings
    lat.add_hoppings(
        # between neighboring cells, between which atoms, and the value
        ([1, 0, 0], 'A', 'A', -t),
        ([0, 1, 0], 'A', 'A', -t),
        ([0, 0, 1], 'A', 'A', -t)
    )
    return lat

def main():
    # Disorder parameters
    W = 0.4                                   # [eV] strength of Anderson disorder
    stddev = W/np.sqrt(12.0)                     # [eV] standard deviation of distribution
    mean = 0.0                                # [eV] average value of the distribution

    # Simulation parameters
    moments = 256
    nx = ny = nz = 2
    lx = ly = lz = 64
    LIM = 7


    lattice = cube()

    # Add disorder 
    disorder = kite.Disorder(lattice)
    disorder.add_disorder("A",'Uniform', mean, stddev)

    # Choose the path in k-space
    Gamma = [0,     0,     0] 
    X     = [np.pi, 0,     0]
    Y     = [0,     np.pi, 0]
    M     = [np.pi, np.pi, 0]
    R     = [np.pi, np.pi, np.pi]
    dk = 0.3
    points = [Gamma, X, Y, M, R, Gamma]
    k_path = pb.results.make_path(*points, step=dk)


    weights = [1]

    configuration = kite.Configuration(
            divisions=[nx, ny, nz], 
            length=[lx, ly, lz], 
            boundaries=["periodic", "periodic", "periodic"], 
            is_complex=True, 
            precision=1, 
            spectrum_range=[-LIM,LIM])

    calculation = kite.Calculation(configuration)
    calculation.arpes(
            k_vector=k_path, 
            num_moments=moments, 
            weight=weights,
            num_disorder=1)

    output_file = "arpes-output.h5"
    kite.config_system(lattice, configuration, calculation, filename=output_file, disorder=disorder)
    return output_file


if __name__ == "__main__":
    hdf5_file = main() # generate the Configuration file

    if len(sys.argv) > 1 and sys.argv[1] == "complete":
        import run_all_examples as ra
        import process_arpes as pa
        ra.run_calculation(hdf5_file)
        ra.run_tools(hdf5_file, options="--ARPES -K green 0.1 -E -10 10 2048 -F 100")
        pa.process_arpes("arpes.dat")


