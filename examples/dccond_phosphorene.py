""" Phosphorene conductivity 'xx'/'yy'

    ##########################################################################
    #                         Copyright 2022, KITE                           #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in eV
    Lattice: Bilayer phosphorene
    Configuration: Periodic boundary conditions, double precision,
                    automatic scaling, size of the system 512x512, with domain decomposition (nx=ny=2)
    Calculation: singleshot_conductivity_dc xx/yy
    Last updated: 18/07/2022
"""

__all__ = ["main"]

import kite
import numpy as np
import pybinding as pb

def phosphorene_lattice(num_hoppings=4):
    """Return lattice specification for a bilayer phosphorene lattice"""

    # parameters
    a = 0.222  # nm
    ax = 0.438  # nm
    ay = 0.332  # nm
    theta = 96.79 * (np.pi / 180)
    phi = 103.69 * (np.pi / 180)
    h = a * np.sin(phi - np.pi / 2)
    s = 0.5 * ax - a * np.cos(theta / 2)
    
    # define lattice vectors
    a1 = np.array([ax, 0])
    a2 = np.array([0, ay])

    # create a lattice with 2 primitive vectors
    lat = pb.Lattice(a1=a1, a2=a2)
    
    # add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('A', [-s / 2, -ay / 2, h], 0),
        ('B', [s / 2, -ay / 2, 0], 0),
        ('C', [-s / 2 + ax / 2, 0, 0], 0),
        ('D', [s / 2 + ax / 2, 0, h], 0)
    )
    
    # Add hopping energies
    lat.register_hopping_energies(
        {'t1': -1.22,
         't2': 3.665,
         't3': -0.205,
         't4': -0.105,
         't5': -0.055}
    )
    
    # Add hoppings
    if num_hoppings < 2:
        raise RuntimeError("t1 and t2 must be included")
    elif num_hoppings > 5:
        raise RuntimeError("t5 is the last one")

    if num_hoppings >= 2:
        lat.add_hoppings(
            # inside the main cell, between which atoms, and the value
            ([0, 0], 'A', 'B', 't2'),
            ([0, 0], 'C', 'D', 't2'),
            # between neighboring cells, between which atoms, and the value
            ([-1, 0], 'A', 'D', 't1'),
            ([-1, -1], 'A', 'D', 't1'),
            ([0, 0], 'B', 'C', 't1'),
            ([0, -1], 'B', 'C', 't1')
        )
    if num_hoppings >= 3:
        lat.add_hoppings(
            # inside the main cell, between which atoms, and the value
            ([0, 0], 'A', 'D', 't3'),
            # between neighboring cells, between which atoms, and the value
            ([0, -1], 'A', 'D', 't3'),
            ([1, 1], 'C', 'B', 't3'),
            ([1, 0], 'C', 'B', 't3')
        )
    if num_hoppings >= 4:
        lat.add_hoppings(
            # inside the main cell, between which atoms, and the value
            ([0, 0], 'A', 'C', 't4'),
            ([0, 0], 'B', 'D', 't4'),
            ([0, -1], 'B', 'D', 't4'),
            # between neighboring cells, between which atoms, and the value
            ([0, -1], 'A', 'C', 't4'),
            ([-1, 0], 'A', 'C', 't4'),
            ([-1, -1], 'A', 'C', 't4'),
            ([-1, 0], 'B', 'D', 't4'),
            ([-1, -1], 'B', 'D', 't4')
        )
    if num_hoppings >= 5:
        lat.add_hoppings(
            # between neighboring cells, between which atoms, and the value
            ([-1, 0], 'A', 'B', 't5'),
            ([-1, 0], 'C', 'D', 't5')
        )
    lat.min_neighbors = 2
    return lat


def main(direction='xx', num_hoppings=4):
    """Prepare the input file for KITEx"""
    # load lattice
    lattice = phosphorene_lattice(num_hoppings=num_hoppings)

    # number of decomposition parts [nx,ny,nz] in each direction of matrix.
    # This divides the lattice into various sections, each of which is calculated in parallel
    nx = ny = 2
    # number of unit cells in each direction.
    lx = ly = 256

    # make config object which caries info about
    # - the number of decomposition parts [nx, ny],
    # - lengths of structure [lx, ly]
    # - boundary conditions [mode,mode, ... ] with modes:
    #   . "periodic"
    #   . "open"
    #   . "twisted" -- this option needs the extra argument angles=[phi_1,..,phi_DIM] where phi_i \in [0, 2*M_PI]
    #   . "random"

    # Boundary Mode
    mode = "periodic"

    # - specify precision of the exported hopping and onsite data, 0 - float, 1 - double, and 2 - long double.
    # - scaling, if None it's automatic, if present select spectrum_range=[e_min, e_max]
    configuration = kite.Configuration(divisions=[nx, ny],
                                       length=[lx, ly],
                                       boundaries=[mode, mode],
                                       is_complex=False,
                                       spectrum_range=[-10, 10],
                                       precision=1)

    # define energy grid
    num_points = 15
    energy = np.linspace(0, 3, num_points)

    # specify calculation type
    calculation = kite.Calculation(configuration)
    # require the calculation of singleshot_conductivity_dc
    calculation.singleshot_conductivity_dc(energy=energy,
                                           num_moments=512,
                                           num_random=1,
                                           num_disorder=1,
                                           direction=direction,
                                           eta=0.02)

    # configure the *.h5 file
    output_file = "ph{0}-output.h5".format(direction)
    kite.config_system(lattice, configuration, calculation, filename=output_file)

    # for generating the desired output from the generated HDF5-file, run
    # ../build/KITEx phxx-output.h5
    # ../tools/build/KITE-tools phxx-output.h5

    # returning the name of the created HDF5-file
    return output_file


def post_process(file_name="phxx-output.h5", out_file_name="condDC.dat"):
    from h5py import File
    np.savetxt(out_file_name, File(file_name, "r+")['Calculation']['singleshot_conductivity_dc']['SingleShot'])


if __name__ == "__main__":
    main()
