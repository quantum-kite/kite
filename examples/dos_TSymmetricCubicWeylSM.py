""" Density of states of a T-Symmetric Cubic Weyl Semimetal with Anderson Disorder
    [J. H. Pixley, P. Goswami, and S. das Sarma PRB 93, 085103 (2016)]

    ##########################################################################
    #                         Copyright 2022, KITE                                                                      #
    #                         Home page: quantum-kite.com                                                       #
    ##########################################################################

    Units: Energy in units of nearest neighbour hopping, |t| = 1
    Lattice: simple cubic lattice (two orbitals per site)
    Configuration: random twisted boundary conditions, double precision, manual rescaling
    Calculation type: Average DOS
    Last updated: 26/07/2022
"""

__all__ = ["main"]

import numpy as np
import pybinding as pb
import kite

def Weyl_Semimetal(t = 1.0):
    # define lattice vectors (Simple Cubic lattice)
    a1 = np.array([1, 0, 0])
    a2 = np.array([0, 1, 0])
    a3 = np.array([0, 0, 1])

    # create a lattice with 3 primitive vectors
    lat = pb.Lattice(a1=a1, a2=a2, a3=a3)

     # Add sublattices
    lat.add_sublattices(('A',[0,0,0],0.0),
                                 ('B',[0,0,0],0.0))

    # Add Nearest-Neighbour hoppings
    lat.add_hoppings(([1,0,0],'A','B', 0.5*t),
                              ([1,0,0],'B','A', 0.5*t),
                              ([0,1,0],'A','B', 0.5*t*1j),
                              ([0,1,0],'B','A',-0.5*t*1j),
                              ([0,0,1],'A','A', 0.5*t),
                              ([0,0,1],'B','B',-0.5*t))
    return lat

def main(Anderson_W = 3.0):
    """Prepare the input file for KITEx"""
    # load lattice
    lattice =  Weyl_Semimetal()
    
    # add scalar on-site Disorder (Box distribution in [W/2,W/2])
    # Number of Samples controlled by num_disorder
    disorder = kite.Disorder(lattice)
    disorder.add_disorder('A', 'Uniform', 0.0, Anderson_W/np.sqrt(12))
    disorder.add_disorder('B', 'Uniform', 0.0, Anderson_W/np.sqrt(12))

    # Phase Diagram
    # Anderson_W < 2.55 --> Semimetal Phase (AvDoS Zero at Zero Energy)
    # Anderson_W > 2.55 --> Diffusive Metal Phase (Finite AvDoS at Zero Energy)
    
    # number of decomposition parts [nx,ny,nz] in each direction of matrix.
    # This divides the lattice into various sections, each of which is calculated in parallel
    nx = 2; ny = nz = 1
    # number of unit cells in each direction.
    lx = ly = lz = 64
    # Boundary Mode
    mode = "random"
    
    # make config object which caries info about
    # - the number of decomposition parts [nx, ny,nz],
    # - lengths of structure [lx, ly,lz]
    # - boundary conditions [mode,mode,mode] with modes:
    #   . "periodic"
    #   . "open"
    #   . "twisted" -- this option needs the extra argument angles=[phi_1,..,phi_DIM] where phi_i \in [0, 2*M_PI]
    #   . "random" -- The number of random twists is fixed by num_disorder
    
    configuration = kite.Configuration(divisions = [nx, ny, nz],
                                           length = [lx, ly, lz],
                                           boundaries=[mode,mode,mode],
                                           is_complex=True,
                                           precision=1,
                                           spectrum_range=[-1.8 - Anderson_W,1.8 + Anderson_W])
    
    # require the calculation of DOS
    calculation = kite.Calculation(configuration)
    calculation.dos(num_points  = 1001,
                        num_moments = 512,
                        num_random = 1,
                        num_disorder = 10)
    
    # configure the *.h5 file
    output_file = "Weyl_Semimetal-output.h5"
    kite.config_system(lattice, configuration, calculation, filename=output_file, disorder=disorder)
    
    # for generating the desired output from the generated HDF5-file, run
    # ../build/KITEx Weyl_Semimetal-output.h5
    # ../tools/build/KITE-tools Weyl_Semimetal-output.h5
    
    # returning the name of the created HDF5-file
    return output_file

if __name__ == "__main__":
    main()

