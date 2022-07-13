""" Run all the examples in the kite/examples-folder

    ##############################################################################
    #                        Copyright 2022, KITE                                #
    #                        Home page: quantum-kite.com                         #
    ##############################################################################

    Last updated: 13/07/2022
"""

import os
import matplotlib.pyplot as plt
import numpy as np


def main():
    # Header
    print("    ------------------------------------------------------------------------------")
    print("       Run all the examples in the kite/examples-folder                           ")
    print("    ------------------------------------------------------------------------------")
    print("                                                                                  ")
    print("    ##############################################################################")
    print("    #                        Copyright 2022, KITE                                #")
    print("    #                        Home page: quantum-kite.com                         #")
    print("    ##############################################################################")
    print("                                                                                  ")

    # Example 1: dos_square_lattice.py
    print("    ======= Example 1: DOS for a square lattice                           ========")
    import dos_square_lattice
    print("          - Making configuration file                                             ")
    dos_square_lattice.main()
    print("          - Running KITEx                                                         ")
    os.system("../build/KITEx square_lattice-data.h5")
    print("          - Running KITE-tools                                                    ")
    os.system("../tools/build/KITE-tools square_lattice-data.h5")
    os.system("mv dos.dat square_lattice-dos.dat")
    print("          - Visualizing the results                                               ")
    square_lattice_dos = np.loadtxt("square_lattice-dos.dat")
    fig = plt.figure()
    ax = fig.subplots()
    ax.plot(square_lattice_dos[:, 0], square_lattice_dos[:, 1])
    ax.set_title("DOS Square Lattice")
    ax.set_xlabel("Energy (|t|)")
    ax.set_ylabel("DOS (A.U.)")
    fig.savefig("square_lattice-dos.pdf")
    plt.close(fig)

    # Example 2: dos_square_lattice_twisted_bc.py
    print("    ======= Example 2: DOS for a square lattice with twisted BC             ======")
    import dos_square_lattice_twisted_bc
    print("          - Making configuration file                                             ")
    dos_square_lattice_twisted_bc.main()
    print("          - Running KITEx                                                         ")
    os.system("../build/KITEx square_lattice_twisted_bc-data.h5")
    print("          - Running KITE-tools                                                    ")
    os.system("../tools/build/KITE-tools square_lattice_twisted_bc-data.h5")
    os.system("mv dos.dat square_lattice_twisted_bc-dos.dat")
    print("          - Visualizing the results                                               ")
    square_lattice_twisted_bc_dos = np.loadtxt("square_lattice_twisted_bc-dos.dat")
    fig = plt.figure()
    ax = fig.subplots()
    ax.plot(square_lattice_twisted_bc_dos[:, 0], square_lattice_twisted_bc_dos[:, 1])
    ax.set_title("DOS Square Lattice & twisted BC")
    ax.set_xlabel("Energy (|t|)")
    ax.set_ylabel("DOS (A.U.)")
    fig.savefig("square_lattice_twisted_bc-dos.pdf")
    plt.close(fig)

    # Example 3: dos_checkerboard_lattice.py
    print("    ======= Example 3: DOS for a checkerboard lattice                       ======")
    import dos_checkerboard_lattice
    print("          - Making configuration file                                             ")
    dos_checkerboard_lattice.main()
    print("          - Running KITEx                                                         ")
    os.system("../build/KITEx checkboard_lattice-data.h5")
    print("          - Running KITE-tools                                                    ")
    os.system("../tools/build/KITE-tools checkboard_lattice-data.h5")
    os.system("mv dos.dat checkboard_lattice-dos.dat")
    print("          - Visualizing the results                                               ")
    checkerboard_lattice_dos = np.loadtxt("checkboard_lattice-dos.dat")
    fig = plt.figure()
    ax = fig.subplots()
    ax.plot(checkerboard_lattice_dos[:, 0], checkerboard_lattice_dos[:, 1])
    ax.set_title("DOS Checkerboard")
    ax.set_xlabel("Energy (|t|)")
    ax.set_ylabel("DOS (A.U.)")
    fig.savefig("checkboard_lattice-dos.pdf")
    plt.close(fig)

    # Example 4: dos_graphene.py
    print("    ======= Example 4: DOS for graphene                                     ======")
    import dos_graphene
    print("          - Making configuration file                                             ")
    dos_graphene.main()
    print("          - Running KITEx                                                         ")
    os.system("../build/KITEx graphene_lattice-data.h5")
    print("          - Running KITE-tools                                                    ")
    os.system("../tools/build/KITE-tools graphene_lattice-data.h5")
    os.system("mv dos.dat graphene_lattice-dos.dat")
    print("          - Visualizing the results                                               ")
    dos_graphene = np.loadtxt("graphene_lattice-dos.dat")
    fig = plt.figure()
    ax = fig.subplots()
    ax.plot(dos_graphene[:, 0], dos_graphene[:, 1])
    ax.set_title("DOS Graphene")
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("DOS (A.U.)")
    fig.savefig("graphene-dos.pdf")
    plt.close(fig)

    # Example 5: dos_cubic_lattice_twisted_bc.py
    print("    ======= Example 5: DOS for a cubic lattice                              ======")
    import dos_cubic_lattice_twisted_bc
    print("          - Making configuration file                                             ")
    dos_cubic_lattice_twisted_bc.main()
    print("          - Running KITEx                                                         ")
    os.system("../build/KITEx cubic_lattice-data.h5")
    print("          - Running KITE-tools                                                    ")
    os.system("../tools/build/KITE-tools cubic_lattice-data.h5")
    os.system("mv dos.dat cubic_lattice-dos.dat")
    print("          - Visualizing the results                                               ")
    cubic_lattice_twisted_bc_dos = np.loadtxt("cubic_lattice-dos.dat")
    fig = plt.figure()
    ax = fig.subplots()
    ax.plot(cubic_lattice_twisted_bc_dos[:, 0], cubic_lattice_twisted_bc_dos[:, 1])
    ax.set_title("DOS Cubic lattice")
    ax.set_xlabel("Energy (|t|)")
    ax.set_ylabel("DOS (A.U.)")
    fig.savefig("cubic_lattice-dos.pdf")
    plt.close(fig)

    # Example 6: basic_on_site_disorder.py
    print("    ======= Example 6: DOS for lattice with onsite disorder                ======")
    import basic_on_site_disorder
    print("          - Making configuration file                                             ")
    basic_on_site_disorder.main()
    print("          - Running KITEx                                                         ")
    os.system("../build/KITEx on_site_disorder-data.h5")
    print("          - Running KITE-tools                                                    ")
    os.system("../tools/build/KITE-tools on_site_disorder-data.h5")
    os.system("mv dos.dat on_site_disorder-dos.dat")
    print("          - Visualizing the results                                               ")
    basic_on_site_disorder_dos = np.loadtxt("on_site_disorder-dos.dat")
    fig = plt.figure()
    ax = fig.subplots()
    ax.plot(basic_on_site_disorder_dos[:, 0], basic_on_site_disorder_dos[:, 1])
    ax.set_title("DOS On-site disorder")
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("DOS (A.U.)")
    fig.savefig("on_site_disorder-dos.pdf")
    plt.close(fig)

    # Example 7: basic_vacancies.py
    print("    ======= Example 7: DOS for lattice with vacancies                       ======")
    import basic_vacancies
    print("          - Making configuration file                                             ")
    basic_vacancies.main()
    print("          - Running KITEx                                                         ")
    os.system("../build/KITEx vacancies-data.h5")
    print("          - Running KITE-tools                                                    ")
    os.system("../tools/build/KITE-tools vacancies-data.h5")
    os.system("mv dos.dat vacancies-dos.dat")
    print("          - Visualizing the results                                               ")
    vacancies_dos = np.loadtxt("vacancies-dos.dat")
    fig = plt.figure()
    ax = fig.subplots()
    ax.plot(vacancies_dos[:, 0], vacancies_dos[:, 1])
    ax.set_title("DOS On-site vacancies")
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("DOS (A.U.)")
    fig.savefig("vacancies-dos.pdf")
    plt.close(fig)

    # Example 8: basic_mixed_disorder.py
    print("    ======= Example 8: DOS for lattice with vacancies and onsite disorder   ======")
    import basic_mixed_disorder
    print("          - Making configuration file                                             ")
    basic_mixed_disorder.main()
    print("          - Running KITEx                                                         ")
    os.system("../build/KITEx mixed_disorder-data.h5")
    print("          - Running KITE-tools                                                    ")
    os.system("../tools/build/KITE-tools mixed_disorder-data.h5")
    os.system("mv dos.dat mixed_disorder-dos.dat")
    print("          - Visualizing the results                                               ")
    mixed_disorder_dos = np.loadtxt("mixed_disorder-dos.dat")
    fig = plt.figure()
    ax = fig.subplots()
    ax.plot(mixed_disorder_dos[:, 0], mixed_disorder_dos[:, 1])
    ax.set_title("DOS On-site vacancies and disorder")
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("DOS (A.U.)")
    fig.savefig("mixed_disorder-dos.pdf")
    plt.close(fig)

    # Example 9: optcond_gaussian_disorder.py
    print("    ======= Example 9: Optical conductivity with onsite disorder            ======")
    import optcond_gaussian_disorder
    print("          - Making configuration file                                             ")
    optcond_gaussian_disorder.main()
    print("          - Running KITEx                                                         ")
    os.system("../build/KITEx optcond_gaussian_disorder-data.h5")
    print("          - Running KITE-tools                                                    ")
    os.system("../tools/build/KITE-tools optcond_gaussian_disorder-data.h5")
    os.system("mv dos.dat optcond_gaussian_disorder-dos.dat")
    os.system("mv optcond.dat optcond_gaussian_disorder-optcond.dat")
    print("          - Visualizing the results                                               ")
    optcond_gaussian_disorder_dos = np.loadtxt("optcond_gaussian_disorder-dos.dat")
    fig = plt.figure()
    ax = fig.subplots()
    ax.plot(optcond_gaussian_disorder_dos[:, 0], optcond_gaussian_disorder_dos[:, 1])
    ax.set_title("DOS On-site disorder")
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("DOS (A.U.)")
    fig.savefig("optcond_gaussian_disorder-dos.pdf")
    plt.close(fig)
    optcond_gaussian_disorder_optcond = np.loadtxt("optcond_gaussian_disorder-optcond.dat")
    fig = plt.figure()
    ax = fig.subplots()
    optr = ax.plot(optcond_gaussian_disorder_optcond[:, 0], optcond_gaussian_disorder_optcond[:, 1])
    opti = ax.plot(optcond_gaussian_disorder_optcond[:, 0], optcond_gaussian_disorder_optcond[:, 2])
    ax.legend([optr[0], opti[0]], [r"$\mathcal{R}[\simga_{xx}]$", r"$\mathcal{I}[\simga_{xx}]$"])
    ax.set_title("Optical conductivity On-site disorder")
    ax.set_xlabel(r"$\hbar \omega$ (eV)")
    ax.set_ylabel("Optical conductivity")
    fig.savefig("optcond_gaussian_disorder-optcond.pdf")
    plt.close(fig)

    # Example 10: haldane.py
    print("    ======= Example 10: DOS & optical conductivity for the Haldene model    ======")
    import haldane
    print("          - Making configuration file                                             ")
    haldane.main()
    print("          - Running KITEx                                                         ")
    os.system("../build/KITEx haldane-data.h5")
    print("          - Running KITE-tools                                                    ")
    os.system("../tools/build/KITE-tools haldane-data.h5")
    os.system("mv dos.dat haldane-dos.dat")
    os.system("mv optcond.dat haldane-optcond.dat")
    print("          - Visualizing the results                                               ")
    haldane_dos = np.loadtxt("haldane-dos.dat")
    fig = plt.figure()
    ax = fig.subplots()
    ax.plot(haldane_dos[:, 0], haldane_dos[:, 1])
    ax.set_title("DOS Haldane")
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("DOS (A.U.)")
    fig.savefig("haldane-dos.pdf")
    plt.close(fig)
    haldane_optcond = np.loadtxt("haldane-optcond.dat")
    fig = plt.figure()
    ax = fig.subplots()
    optr = ax.plot(haldane_optcond[:, 0], haldane_optcond[:, 1])
    opti = ax.plot(haldane_optcond[:, 0], haldane_optcond[:, 2])
    ax.legend([optr[0], opti[0]], [r"$\mathcal{R}[\simga_{xx}]$", r"$\mathcal{I}[\simga_{xx}]$"])
    ax.set_title("Optical conductivity haldane")
    ax.set_xlabel(r"$\hbar \omega$ (eV)")
    ax.set_ylabel("Optical conductivity")
    fig.savefig("haldane-optcond.pdf")
    plt.close(fig)

    # Example 11: phosphorene.py
    print("    ======= Example 11: DOS & DC coductivity for phosphorene in XX          ======")
    import phosphorene
    print("          - Making configuration file                                             ")
    phosphorene.main()
    print("          - Running KITEx                                                         ")
    os.system("../build/KITEx phxx-data.h5")
    print("          - Running KITE-tools                                                    ")
    os.system("../tools/build/KITE-tools phxx-data.h5")
    os.system("mv condDC.dat phxx-condDC.dat")
    print("          - Visualizing the results                                               ")
    phxx_conddc = np.loadtxt("phxx-condDC.dat")
    fig = plt.figure()
    ax = fig.subplots()
    ax.plot(phxx_conddc[:, 0], phxx_conddc[:, 1])
    ax.set_title("DC Conductivity XX phosphorene")
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("Optical conductivity")
    fig.savefig("phxx-condDC.pdf")
    plt.close(fig)

    # Example 12: phosphorene.py
    print("    ======= Example 12: DOS & DC coductivity for phosphorene in YY          ======")
    print("          - Making configuration file                                             ")
    phosphorene.main(direction='yy')
    print("          - Running KITEx                                                         ")
    os.system("../build/KITEx phyy.h5")
    print("          - Running KITE-tools                                                    ")
    os.system("../tools/build/KITE-tools phyy.h5")
    os.system("mv condDC.dat phyy-condDC.dat")
    print("          - Visualizing the results                                               ")
    phyy_conddc = np.loadtxt("phyy-condDC.dat")
    fig = plt.figure()
    ax = fig.subplots()
    ax.plot(phyy_conddc[:, 0], phyy_conddc[:, 1])
    ax.set_title("DC Conductivity YY phosphorene")
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("Optical conductivity")
    fig.savefig("phyy-condDC.pdf")
    plt.close(fig)

    # Example 13: twisted_bilayer.py
    print("    ======= Example 12: DOS for twisted bilayer graphene at 21.787 degrees  ======")
    import twisted_bilayer
    print("          - Making configuration file                                             ")
    twisted_bilayer.main(angle_index=3)
    print("          - Running KITEx                                                         ")
    os.system("../build/KITEx tblg_3-data.h5")
    print("          - Running KITE-tools                                                    ")
    os.system("../tools/build/KITE-tools tblg_3-data.h5")
    os.system("mv dos.dat tblg_3-dos.dat")
    print("          - Visualizing the results                                               ")
    tblg_3_dos = np.loadtxt("tblg_3-dos.dat")
    fig = plt.figure()
    ax = fig.subplots()
    ax.plot(mixed_disorder_dos[:, 0], mixed_disorder_dos[:, 1])
    ax.set_title("DOS Twisted Bilayer graphene at 21.787 degrees")
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("DOS (A.U.)")
    fig.savefig("tblg_3-dos.pdf")
    plt.close(fig)


def clean():
    # Header
    print("    ------------------------------------------------------------------------------")
    print("       Remove the previous calculations from the examples folder                  ")
    print("    ------------------------------------------------------------------------------")

    # Example 1
    os.system("rm square_lattice-data.h5 square_lattice-dos.dat square_lattice-dos.pdf")

    # Example 2
    os.system("rm square_lattice_twisted_bc-data.h5 square_lattice_twisted_bc-dos.dat")
    os.system("rm square_lattice_twisted_bc-dos.pdf")

    # Example3
    os.system("rm checkboard_lattice-data.h5 checkboard_lattice-dos.dat checkboard_lattice-dos.pdf")

    # Example 4
    os.system("rm graphene_lattice-data.h5 graphene_lattice-dos.dat graphene-dos.pdf")

    # Example 5
    os.system("rm cubic_lattice-data.h5 cubic_lattice-dos.dat cubic_lattice-dos.pdf")

    # Example 6
    os.system("rm on_site_disorder-data.h5 on_site_disorder-dos.dat on_site_disorder-dos.pdf")

    # Example 7
    os.system("rm vacancies-data.h5 vacancies-dos.dat vacancies-dos.pdf")

    # Example 8
    os.system("rm mixed_disorder-data.h5 mixed_disorder-dos.dat mixed_disorder-dos.pdf")

    # Example 9
    os.system("rm optcond_gaussian_disorder-data.h5 optcond_gaussian_disorder-dos.pdf")
    os.system("rm optcond_gaussian_disorder-optcond.dat optcond_gaussian_disorder-optcond.pdf")

    # Example 10
    os.system("rm haldane-data.h5 haldane-dos.dat haldane-optcond.dat haldane-dos.pdf haldane-optcond.pdf")

    # Example 11
    os.system("rm phxx-data.h5 phxx-condDC.dat phxx-condDC.pdf")

    # Example 12
    os.system("rm phyy-data.h5 phyy-condDC.dat phyy-condDC.pdf")

    # Example 13
    os.system("rm tblg_3-data.h5 tblg_3-dos.dat tblg_3-dos.pdf")


if __name__ == "__main__":
    main()
