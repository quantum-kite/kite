""" Optical Conductivity of Hexagonal Boron Nitride

    ##########################################################################
    #                         Copyright 2020/22, KITE                        #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Units: Energy in eV
    Lattice: Honeycomb
    Disorder: Vacancies
    Configuration: Periodic boundary conditions, double precision,
                    given rescaling, size of the system flexible, with domain decomposition (nx=ny=2)
    Calculation type: Optical Conductivity
    Last updated: 08/05/2025
"""
import pybinding as pb
import numpy as np
import kite


def hbn():
    t = -1.0
    a_cc = 1                              # [nm] carbon-carbon distance
    a = a_cc*np.sqrt(3)                   # [nm] unit cell length
    gap = 0.2                             # [t] gap

    lat = pb.Lattice(
        a1=[a/2,  a/2 * np.sqrt(3)],
        a2=[a/2, -a/2 * np.sqrt(3)])

    lat.add_sublattices(
    ('A', [0,    0], -gap/2.0),
    ('B', [0, a_cc],  gap/2.0)
    )

    lat.add_hoppings(
    ([ 0, 0], 'A', 'B', t),
    ([ 0, 1], 'A', 'B', t),
    ([-1, 0], 'A', 'B', t)
    )


    return lat


def main():
    lattice = hbn()

    N = 256         # number of polynomials
    C = 0.02        # concentration of vacancies
    lx = ly = 128   # system dimensions
    nx, ny = 2, 2   # number of threads in each direction
    LIM = 3.2       # overestimated bounds of the spectrum
    output_file = "hbn-output.h5"

    disorder = kite.Disorder(lattice)
    struct_A = kite.StructuralDisorder(lattice, concentration=C)
    struct_B = kite.StructuralDisorder(lattice, concentration=C)
    struct_A.add_vacancy('A')
    struct_B.add_vacancy('B')
    disorder_structural = [struct_A, struct_B]


    configuration = kite.Configuration(
            divisions=[nx, ny],
            length=[lx, ly],
            boundaries=["periodic", "periodic"],
            is_complex=True,
            precision=1,
            spectrum_range=[-LIM,LIM])

    calculation = kite.Calculation(configuration)
    calculation.conductivity_optical_nonlinear(
            "xxy",
            num_points=512,
            num_moments=N,
            num_random=1,
            num_disorder=1,
            special=1
            )

    kite.config_system(lattice, configuration, calculation, filename=output_file, disorder_structural=disorder_structural)
    return output_file

if __name__ == "__main__":
    main()
