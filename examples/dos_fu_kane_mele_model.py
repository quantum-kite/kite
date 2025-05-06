""" Density of states of the Fu-Kane-Mele Topological Insulator in the Diamond Lattice
    [L. Fu, C.L. Kane and E. J. Mele, PRL 98, 106803 (2007)]

    ##########################################################################
    #                         Copyright 2022, KITE                           #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in units of fixed nearest neighbour hopping, |t| = |lambda_SO| = 1
    Lattice: face-centered cubic lattice (Diamond Structure)
    Configuration: random twisted boundary conditions, double precision, automatic rescaling
    Calculation type: Average DOS
    Last updated: 28/07/2022
"""

__all__ = ["main"]

import kite
import numpy as np
import pybinding as pb


def fu_kane_mele(dt_plus, dt_minus, lambda_so=1, t=1):
    """Retrun FuKaneMele Hamiltonian"""

    # parameters
    lso = -t * lambda_so * 1j / 8.0
    delta_t1 = 0.5 * (dt_plus + dt_minus)  # delta1 Hopping Deformation
    delta_t2 = 0.5 * (dt_plus - dt_minus)  # delta2 Hopping Deformation

    # define lattice vectors (FCC lattice)
    a1 = np.array([0.5, 0.5, 0.0])
    a2 = np.array([0.0, 0.5, 0.5])
    a3 = np.array([0.5, 0.0, 0.5])

    # create a lattice with 3 primitive vectors
    lat = pb.Lattice(a1=a1, a2=a2, a3=a3)

    # add sublattices (Includes Spin-1/2 Degree of Freedom)
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A1', [0.00, 0.00, 0.00], 0.0),
        ('A2', [0.00, 0.00, 0.00], 0.0),
        ('B1', [0.25, 0.25, 0.25], 0.0),
        ('B2', [0.25, 0.25, 0.25], 0.0)
    )

    # Add Nearest-Neighbour hoppings (2 Deformed hoppings + 2 Normal Hoppings)
    lat.add_hoppings(
        # between neighboring cells, between which atoms, and the value
        ([0, 0, 0], 'A1', 'B1', t + delta_t1),
        ([0, 0, 0], 'A2', 'B2', t + delta_t1),
        ([1, 0, 0], 'A1', 'B1', t + delta_t2),
        ([1, 0, 0], 'A2', 'B2', t + delta_t2),
        ([0, 1, 0], 'A1', 'B1', t),
        ([0, 1, 0], 'A2', 'B2', t),
        ([0, 0, 1], 'A1', 'B1', t),
        ([0, 0, 1], 'A2', 'B2', t)
    )
    
    # Add Next-Nearest-Neighbour hoppings (Spin-Orbit Coupling)
    deltas = np.array(
        [[1, 0, 0],
         [0, 1, 0],
         [0, 0, 1],
         [-1, 1, 0],
         [0, -1, 1],
         [1, 0, -1]],
        dtype=int
    )
    
    gg = np.array(
        [[-1, 1, 0],
         [0, -1, 1],
         [1, 0, -1],
         [-1, 0, -1],
         [-1, -1, 0],
         [0, -1, -1]],
        dtype=int
    )
    
    for i in range(len(deltas[:, 0])):
        lat.add_hoppings(
            # between neighboring cells, between which atoms, and the value
            (deltas[i, :], 'A1', 'A1',  lso * gg[i, 2]),
            (deltas[i, :], 'A2', 'A2', -lso * gg[i, 2]),
            (deltas[i, :], 'A1', 'A2',  lso * (gg[i, 0] - 1j * gg[i, 1])),
            (deltas[i, :], 'A2', 'A1',  lso * (gg[i, 0] + 1j * gg[i, 1])),
            (deltas[i, :], 'B1', 'B1', -lso * gg[i, 2]),
            (deltas[i, :], 'B2', 'B2',  lso * gg[i, 2]),
            (deltas[i, :], 'B1', 'B2', -lso * (gg[i, 0] - 1j * gg[i, 1])),
            (deltas[i, :], 'B2', 'B1', -lso * (gg[i, 0] + 1j * gg[i, 1]))
        )
    return lat


def main(dt_plus=0.1, dt_minus=0.1):
    """Prepare the input file for KITEx"""
    # load lattice
    lattice = fu_kane_mele(dt_plus, dt_minus)

    # Phase Diagram
    #    - dt_plus > 0 and dt_minus < 0 --> Strong TI 1;(1,1b,1b)
    #    - dt_plus > 0 and dt_minus > 0 --> Strong TI 1;(1,1,1)
    #    - dt_plus < 0 and dt_minus < 0 --> Weak TI 0;(1,1b,1b)
    #    - dt_plus < 0 and dt_minus > 0 --> Weak TI 0;(1,1,1)
    #    - Transitions --> Dirac Semimetal
    
    # number of decomposition parts [nx,ny,nz] in each direction of matrix.
    # This divides the lattice into various sections, each of which is calculated in parallel
    nx = 2
    ny = nz = 1
    # number of unit cells in each direction.
    lx = ly = lz = 64

    # make config object which caries info about
    # - the number of decomposition parts [nx, ny,nz],
    # - lengths of structure [lx, ly,lz]
    # - boundary conditions [mode,mode,mode] with modes:
    #   . "periodic"
    #   . "open"
    #   . "twisted" -- this option needs the extra argument angles=[phi_1,..,phi_DIM] where phi_i \in [0, 2*M_PI]
    #   . "random" -- The number of random twists is fixed by num_disorder

    # Boundary Mode
    mode = "random"

    # - specify precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
    # - scaling, if None it's automatic, if present select spectrum_range=[e_min, e_max]
    configuration = kite.Configuration(
        divisions=[nx, ny, nz],
        length=[lx, ly, lz],
        boundaries=[mode, mode, mode],
        is_complex=True,
        precision=1
    )

    # configure the *.h5 file
    calculation = kite.Calculation(configuration)
    calculation.dos(
        num_points=1001,
        num_moments=1024,
        num_random=1,
        num_disorder=5
    )

    # configure the *.h5 file
    output_file = "FKM_Model-output.h5"
    kite.config_system(lattice, configuration, calculation, filename=output_file)

    # for generating the desired output from the generated HDF5-file, run
    # ../build/KITEx FKM_Model-output.h5
    # ../build/KITE-tools FKM_Model-output.h5

    # returning the name of the created HDF5-file
    return output_file


if __name__ == "__main__":
    main()
