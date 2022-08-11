""" Run all the examples in the kite/examples-folder

    ##########################################################################
    #                         Copyright 2022, KITE                           #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Last updated: 30/07/2022
"""

import numpy as np
import matplotlib.pyplot as plt
import process_arpes as pa
from os import system as terminal
from os.path import exists


_kitex_dir = "../build"
_kite_tools_dir = "../tools/build"

KITEx_exists = exists(_kitex_dir)
tools_exists = exists(_kite_tools_dir)

if not KITEx_exists:
    print("Please make sure the KITEx executable is in the correct place relative to this script: ../build/KITEx")
    exit(1)
if not tools_exists:
    print("Please make sure the KITE-tools executable is in the correct place relative to this script: ../tools/build/KITE-tools")
    exit(1)


def run_calculation(input_file="output.h5"):
    """Run KITEx"""
    print_command("- - - -             Doing the KITEx-calculation                - - - - -")
    terminal("{0}/KITEx {1}".format(_kitex_dir, input_file))


def run_tools(input_file="output.h5", options=""):
    """Run KITE-tools to obtain the DOS"""
    print_command("- - - -             Doing the KITE-tools postprocessing        - - - - -")
    terminal("{0}/KITE-tools {1} {2}".format(_kite_tools_dir, input_file, options))


def make_figure_dos(file_data="dos.dat", title="DOS", xlabel="Energy (ev)", ylabel="DOS (1/eV)", file_out="dos.pdf"):
    """Make a figure for the DOS"""
    dos = np.loadtxt(file_data)
    fig = plt.figure()
    ax = fig.subplots()
    ax.plot(dos[:, 0], dos[:, 1])
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.savefig(file_out)
    plt.close(fig)


def make_figure_opt_cond(file_data="optcond.dat", title="Optical conductivity", xlabel=r"$\hbar\omega$ (ev)",
                         ylabel=r"$\sigma_{xx}(\omega)$", file_out="optcond.pdf"):
    """Make a figure for the DOS"""
    optcond = np.loadtxt(file_data)
    fig = plt.figure()
    ax = fig.subplots()
    lines = [ax.plot(optcond[:, 0], optcond[:, 1])[0], ax.plot(optcond[:, 0], optcond[:, 2])[0]]
    ax.legend(lines, [r"$\mathcal{R}[\sigma_{xy}]$", r"$\mathcal{I}[\sigma_{xy}]$"])
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.savefig(file_out)
    plt.close(fig)


def make_figure_cond_dc(file_data="condDC.dat", title="DC conductivity", xlabel="E (eV)",
                         ylabel=r"$\sigma (2e^2/h)$", file_out="optcond.pdf"):
    """Make a figure for the DOS"""
    optcond = np.loadtxt(file_data)
    fig = plt.figure()
    ax = fig.subplots()
    lines = [ax.plot(optcond[:, 0], optcond[:, 1])[0], ax.plot(optcond[:, 0], optcond[:, 2])[0]]
    ax.legend(lines, [r"$\sigma_{xx}$", r"$\sigma_{xy}$"])
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.savefig(file_out)
    plt.close(fig)


def make_figure_cond_dc_ss(file_data="condDC.dat", title="DC conductivity", xlabel="E (eV)",
                         ylabel=r"$\sigma (2e^2/h)$", file_out="optcond.pdf"):
    """Make a figure for the DOS"""
    optcond = np.loadtxt(file_data)
    fig = plt.figure()
    ax = fig.subplots()
    lines = [ax.plot(optcond[:, 0], optcond[:, 2])[0]]
    ax.legend(lines, [r"$\sigma_{xx}$"])
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.savefig(file_out)
    plt.close(fig)


def print_command(text=""):
    print("\033[94m{0}\033[0m".format(text))


def print_title(text=""):
    print_command("-" * 74 + "\n" + " " * 10 + text + "\n" + "-" * 74 + "\n")


def main(selection=None):
    if selection is None:
        selection = np.arange(17, dtype=np.int) + 1
    import matplotlib as mpl
    import seaborn as sns

    mpl.rcParams['figure.dpi'] = 100
    mpl.rcParams['savefig.dpi'] = 100
    sns.set_style("white")

    colors = ["dusty purple", "faded green","windows blue", "amber", "greyish"]
    current_palette = sns.xkcd_palette(colors)
    sns.set_palette(current_palette)
    sns.set_style("ticks")
    sns.set_context("talk", font_scale=1.1)

    # Header
    print_command("##########################################################################")
    print_command("#                         Copyright 2022, KITE                           #")
    print_command("#                         Home page: quantum-kite.com                    #")
    print_command("##########################################################################")
    print_command()
    print_title("Run all the examples in the kite/examples-folder")

    if 1 in selection:
        # Example 1: dos_square_lattice.py
        print_command("======= Example 1: DOS for a square lattice                      =========")
        import dos_square_lattice as example
        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main()
        pre_file_name = hdf5_file.replace("-output.h5", "")
        dos_data = "{0}-dos.dat".format(pre_file_name)
        dos_figure = "{0}-dos.pdf".format(pre_file_name)
        run_calculation(hdf5_file)
        run_tools(hdf5_file)
        terminal("mv dos.dat {0}".format(dos_data))
        make_figure_dos(file_data=dos_data, title="DOS Square lattice", xlabel="Energy (|t|)",
                        file_out=dos_figure)

    if 2 in selection:
        # Example 2: dos_cube.py
        print_command("======= Example 2: DOS for a cubic lattice                       =========")
        import dos_cube as example
        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main()
        pre_file_name = hdf5_file.replace("-output.h5", "")
        dos_data = "{0}-dos.dat".format(pre_file_name)
        dos_figure = "{0}-dos.pdf".format(pre_file_name)
        run_calculation(hdf5_file)
        run_tools(hdf5_file)
        terminal("mv dos.dat {0}".format(dos_data))
        make_figure_dos(file_data=dos_data, title="DOS cube lattice", xlabel="Energy (|t|)",
                        file_out=dos_figure)

    if 3 in selection:
        # Example 3: dos_square_lattice_twisted_bc.py
        print_command("    ======= Example 3: DOS for a square lattice with twisted BC  =========")
        import dos_square_lattice_twisted_bc as example
        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main()
        pre_file_name = hdf5_file.replace("-output.h5", "")
        dos_data = "{0}-dos.dat".format(pre_file_name)
        dos_figure = "{0}-dos.pdf".format(pre_file_name)
        run_calculation(hdf5_file)
        run_tools(hdf5_file)
        terminal("mv dos.dat {0}".format(dos_data))
        make_figure_dos(file_data=dos_data, title="DOS Square lattice & twisted BC", xlabel="Energy (|t|)",
                        file_out=dos_figure)

    if 4 in selection:
        # Example 4: dos_checkerboard_lattice.py
        print_command("======= Example 4: DOS for a checkerboard lattice                =========")
        import dos_checkerboard_lattice as example
        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main()
        pre_file_name = hdf5_file.replace("-output.h5", "")
        dos_data = "{0}-dos.dat".format(pre_file_name)
        dos_figure = "{0}-dos.pdf".format(pre_file_name)
        run_calculation(hdf5_file)
        run_tools(hdf5_file)
        terminal("mv dos.dat {0}".format(dos_data))
        make_figure_dos(file_data=dos_data, title="DOS Checkerboard lattice", xlabel="Energy (|t|)",
                        file_out=dos_figure)

    if 5 in selection:
        # Example 5: dos_optcond_graphene.py
        print_command("======= Example 5: DOS for graphene                              =========")
        import dos_optcond_graphene as example
        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main()
        pre_file_name = hdf5_file.replace("-output.h5", "")
        dos_data = "{0}-dos.dat".format(pre_file_name)
        dos_figure = "{0}-dos.pdf".format(pre_file_name)
        run_calculation(hdf5_file)
        run_tools(hdf5_file)
        terminal("mv dos.dat {0}".format(dos_data))
        make_figure_dos(file_data=dos_data, title="DOS Graphene",
                        file_out=dos_figure)

    if 6 in selection:
        # Example 6: ldos_graphene.py
        print_command("======= Example 6: LDOS for graphene                             =========")
        import ldos_graphene as example
        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main()
        pre_file_name = hdf5_file.replace("-output.h5", "")
        run_calculation(hdf5_file)
        run_tools(hdf5_file)
        from pybinding.repository.graphene import monolayer as lattice
        for energy in np.linspace(-1, 1, 101):
            print_command("Making a figure for LDOS at energy {energy:.6f}".format(energy=energy))
            ldos_data = "{0}.dat".format(pre_file_name)
            terminal("mv ldos{energy:.6f}.dat ldos-{energy:.6f}{data}".format(energy=energy, data=ldos_data))
            example.analyze_results("ldos-{energy:.6f}{data}".format(energy=energy, data=ldos_data), lattice())

    if 7 in selection:
        # Example 7: dos_cubic_lattice_twisted_bc.py
        print_command("======= Example 7: DOS for a cubic lattice                       =========")
        import dos_cubic_lattice_twisted_bc as example
        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main()
        pre_file_name = hdf5_file.replace("-output.h5", "")
        dos_data = "{0}-dos.dat".format(pre_file_name)
        dos_figure = "{0}-dos.pdf".format(pre_file_name)
        run_calculation(hdf5_file)
        run_tools(hdf5_file)
        terminal("mv dos.dat {0}".format(dos_data))
        make_figure_dos(file_data=dos_data, title="DOS Cubic lattice", xlabel="Energy (|t|)",
                        file_out=dos_figure)

    if 8 in selection:
        # Example 8: dos_square_lattice_disorder.py
        print_command("======= Example 8: DOS for a square lattice with on-site disorder ========")
        import dos_square_lattice_disorder as example
        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main()
        pre_file_name = hdf5_file.replace("-output.h5", "")
        dos_data = "{0}-dos.dat".format(pre_file_name)
        dos_figure = "{0}-dos.pdf".format(pre_file_name)
        run_calculation(hdf5_file)
        run_tools(hdf5_file)
        terminal("mv dos.dat {0}".format(dos_data))
        make_figure_dos(file_data=dos_data, title="DOS square lattice with uniform on-site disorder",
                        file_out=dos_figure)

    if 9 in selection:
        # Example 9: dos_on_site_disorder.py
        print_command("======= Example 9: DOS for a graphene lattice with on-site disorder ======")
        import dos_on_site_disorder as example
        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main()
        pre_file_name = hdf5_file.replace("-output.h5", "")
        dos_data = "{0}-dos.dat".format(pre_file_name)
        dos_figure = "{0}-dos.pdf".format(pre_file_name)
        run_calculation(hdf5_file)
        run_tools(hdf5_file)
        terminal("mv dos.dat {0}".format(dos_data))
        make_figure_dos(file_data=dos_data, title="DOS Graphene with mixed on-site disorder",
                        file_out=dos_figure)

    if 10 in selection:
        # Example 10: dos_vacancies.py
        print_command("======= Example 10: DOS for lattice with vacancies               =========")
        import dos_vacancies as example
        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main()
        pre_file_name = hdf5_file.replace("-output.h5", "")
        dos_data = "{0}-dos.dat".format(pre_file_name)
        dos_figure = "{0}-dos.pdf".format(pre_file_name)
        run_calculation(hdf5_file)
        run_tools(hdf5_file)
        terminal("mv dos.dat {0}".format(dos_data))
        make_figure_dos(file_data=dos_data, title="DOS Honeycomb lattice with vacancies", xlabel="Energy (|t|)",
                        file_out=dos_figure)

    if 11 in selection:
        # Example 11: dos_mixed_disorder.py
        print_command("======= Example 11: DOS for lattice with vacancies and onsite disorder ===")
        import dos_mixed_disorder as example
        print_command("- - - -            Making the configuration file                -  - - - -")
        hdf5_file = example.main()
        pre_file_name = hdf5_file.replace("-output.h5", "")
        dos_data = "{0}-dos.dat".format(pre_file_name)
        dos_figure = "{0}-dos.pdf".format(pre_file_name)
        run_calculation(hdf5_file)
        run_tools(hdf5_file)
        terminal("mv dos.dat {0}".format(dos_data))
        make_figure_dos(file_data=dos_data, title="DOS Graphene with vacancies & on-site disorder",
                        file_out=dos_figure)

    if 12 in selection:
        # Example 12: dos_optcond_gaussian_disorder.py
        print_command("======= Example 12: Optical conductivity with onsite disorder    =========")
        import dos_optcond_gaussian_disorder as example
        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main()
        pre_file_name = hdf5_file.replace("-output.h5", "")
        dos_data = "{0}-dos.dat".format(pre_file_name)
        opt_cond_data = "{0}-optcond.dat".format(pre_file_name)
        dos_figure = "{0}-dos.pdf".format(pre_file_name)
        opt_cond_figure = "{0}-optcond.pdf".format(pre_file_name)
        run_calculation(hdf5_file)
        run_tools(hdf5_file)
        terminal("mv dos.dat {0}".format(dos_data))
        terminal("mv optcond.dat {0}".format(opt_cond_data))
        make_figure_dos(file_data=dos_data, title="DOS Graphene Gaussian On-site disorder",
                        file_out=dos_figure)
        make_figure_opt_cond(file_data=opt_cond_data, title="Optical Conductivity Graphene",
                             file_out=opt_cond_figure)

    if 13 in selection:
        # Example 13: dos_dccond_square_lattice.py
        print_command("======= Example 13: DOS & DC conductivity for the a square lattice =======")
        import dos_dccond_square_lattice as example
        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main()
        pre_file_name = hdf5_file.replace("-output.h5", "")
        dos_data = "{0}-dos.dat".format(pre_file_name)
        cond_dc_data = "{0}-condDC.dat".format(pre_file_name)
        dos_figure = "{0}-dos.pdf".format(pre_file_name)
        cond_dc_figure = "{0}-condDC.pdf".format(pre_file_name)
        run_calculation(hdf5_file)
        run_tools(hdf5_file, "--CondDC -F -4.5 -2 512 -E 4096")
        terminal("mv dos.dat {0}".format(dos_data))
        terminal("mv condDC.dat {0}".format(cond_dc_data))
        make_figure_dos(file_data=dos_data, title="DOS Square Lattice",
                        file_out=dos_figure)
        make_figure_cond_dc(file_data=cond_dc_data, title="DC Conductivity Square Lattice",
                            file_out=cond_dc_figure)

    if 14 in selection:
        # Example 14: dos_dccond_haldane.py
        print_command("======= Example 14: DOS & DC conductivity for the Haldene model  =========")
        import dos_dccond_haldane as example
        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main()
        pre_file_name = hdf5_file.replace("-output.h5", "")
        dos_data = "{0}-dos.dat".format(pre_file_name)
        cond_dc_data = "{0}-condDC.dat".format(pre_file_name)
        dos_figure = "{0}-dos.pdf".format(pre_file_name)
        cond_dc_figure = "{0}-condDC.pdf".format(pre_file_name)
        run_calculation(hdf5_file)
        run_tools(hdf5_file)
        terminal("mv dos.dat {0}".format(dos_data))
        terminal("mv condDC.dat {0}".format(cond_dc_data))
        make_figure_dos(file_data=dos_data, title="DOS Haldane",
                        file_out=dos_figure)
        make_figure_cond_dc(file_data=cond_dc_data, title="DC Conductivity Haldane",
                            file_out=cond_dc_figure)

    if 15 in selection:
        # Example 15: dccond_phosphorene.py
        print_command("======= Example 15: DC coductivity for phosphorene in XX         =========")
        import dccond_phosphorene as example
        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main()
        pre_file_name = hdf5_file.replace("-output.h5", "")
        cond_dc_data = "{0}-condDC.dat".format(pre_file_name)
        cond_dc_figure = "{0}-condDC.pdf".format(pre_file_name)
        run_calculation(hdf5_file)
        from process_single_shot import post_process_singleshot_conductivity_dc as post_process_dccond
        post_process_dccond(hdf5_file)
        terminal("mv {hdf5} {data}".format(hdf5=hdf5_file[:-3] + ".dat", data=cond_dc_data))
        make_figure_cond_dc_ss(file_data=cond_dc_data, title="DC Conductivity Phosphorene XX",
                               file_out=cond_dc_figure)

    if 16 in selection:
        # Example 16: dccond_phosphorene.py
        print_command("======= Example 16: DC coductivity for phosphorene in YY         =========")
        import dccond_phosphorene as example
        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main(direction='yy')
        pre_file_name = hdf5_file.replace("-output.h5", "")
        cond_dc_data = "{0}-condDC.dat".format(pre_file_name)
        cond_dc_figure = "{0}-condDC.pdf".format(pre_file_name)
        run_calculation(hdf5_file)
        from process_single_shot import post_process_singleshot_conductivity_dc as post_process_dccond
        post_process_dccond(hdf5_file)
        terminal("mv {hdf5} {data}".format(hdf5=hdf5_file[:-3] + ".dat", data=cond_dc_data))
        make_figure_cond_dc_ss(file_data=cond_dc_data, title="DC Conductivity Phosphorene YY",
                               file_out=cond_dc_figure)

    if 17 in selection:
        # Example 17: dos_twisted_bilayer.py
        print_command("======= Example 17: DOS for twisted bilayer graphene at 21.787 degrees ===")
        import dos_twisted_bilayer as example
        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main()
        pre_file_name = hdf5_file.replace("-output.h5", "")
        dos_data = "{0}-dos.dat".format(pre_file_name)
        dos_figure = "{0}-dos.pdf".format(pre_file_name)
        run_calculation(hdf5_file)
        run_tools(hdf5_file)
        terminal("mv dos.dat {0}".format(dos_data))
        make_figure_dos(file_data=dos_data, title="DOS Twisted Bilayer graphene at 21.787 degrees",
                        file_out=dos_figure)
    if 18 in selection:
        # Example 18: dos_t_symmetric_cubic_weyl_sm.py
        print_command("======= Example 18: DOS for T symmetric weyl sm                    =======")
        import dos_t_symmetric_cubic_weyl_sm as example
        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main()
        pre_file_name = hdf5_file.replace("-output.h5", "")
        dos_data = "{0}-dos.dat".format(pre_file_name)
        dos_figure = "{0}-dos.pdf".format(pre_file_name)
        run_calculation(hdf5_file)
        run_tools(hdf5_file)
        terminal("mv dos.dat {0}".format(dos_data))
        make_figure_dos(file_data=dos_data, title="DOS for T symmetric weyl semi-metal",
                        file_out=dos_figure)

    if 19 in selection:
        # Example 19: optcond_t_symmetric_cubic_weyl_sm.py
        print_command("======= Example 19: Optical confuctivity  for T symmetric weyl sm    =====")
        import optcond_t_symmetric_cubic_weyl_sm as example
        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main()
        pre_file_name = hdf5_file.replace("-output.h5", "")
        opt_cond_data = "{0}-optcond.dat".format(pre_file_name)
        opt_cond_figure = "{0}-optcond.pdf".format(pre_file_name)
        run_calculation(hdf5_file)
        run_tools(hdf5_file)
        terminal("mv optcond.dat {0}".format(opt_cond_data))
        make_figure_opt_cond(file_data=opt_cond_data, title="Optical Conductivity for T symmetric weyl semi-metal",
                             file_out=opt_cond_figure)

    if 20 in selection:
        # Example 20: dos_fu_kane_mele_model.py
        print_command("======= Example 20: DOS Fu-Kane-Mele model                         =======")
        import dos_fu_kane_mele_model as example
        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main()
        pre_file_name = hdf5_file.replace("-output.h5", "")
        dos_data = "{0}-dos.dat".format(pre_file_name)
        dos_figure = "{0}-dos.pdf".format(pre_file_name)
        run_calculation(hdf5_file)
        run_tools(hdf5_file)
        terminal("mv dos.dat {0}".format(dos_data))
        make_figure_dos(file_data=dos_data, title="DOS Fu-Kane-Mele model",
                        file_out=dos_figure)

    if 21 in selection:
        # Example 21: ARPES in bilayer graphene
        print_command("======= Example 21: ARPES in bilayer graphene model                         =======")
        import arpes_bilayer as example

        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main()

        pre_file_name = hdf5_file.replace("-output.h5", "")
        arpes_data = "{0}-image.png".format(pre_file_name)

        run_calculation(hdf5_file)
        run_tools(hdf5_file, options="--ARPES -K green 0.1 -E -10 10 2048 -F 100")
        pa.process_arpes("arpes.dat")
        terminal("mv arpes.png {0}".format(arpes_data))

    if 22 in selection:
        # Example 22: ARPES in cubic lattice
        print_command("======= Example 22: ARPES in cubic lattice model                         =======")
        import arpes_cubic as example

        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main()

        pre_file_name = hdf5_file.replace("-output.h5", "")
        arpes_data = "{0}-image.png".format(pre_file_name)

        run_calculation(hdf5_file)
        run_tools(hdf5_file, options="--ARPES -K green 0.1 -E -10 10 2048 -F 100")
        pa.process_arpes("arpes.dat")
        terminal("mv arpes.png {0}".format(arpes_data))

    if 23 in selection:
        # Example 23: second-order optical conductivity
        print_command("======= Example 23: second-order photoconductivity in hexagonal Boron Nitride =======")
        import hbn_optcond2_vacancies as example

        print_command("- - - -            Making the configuration file                 - - - - -")
        hdf5_file = example.main()

        pre_file_name = hdf5_file.replace("-output.h5", "")
        opt_cond_data = "{0}-optcond2.dat".format(pre_file_name)
        opt_cond_figure = "{0}-optcond2.pdf".format(pre_file_name)

        run_calculation(hdf5_file)
        run_tools(hdf5_file, options="--CondOpt2 -E 512 -R -1.0 -T 0.0001 -S 0.03 -O 0 4.0 256 -F 0.2")

        terminal("mv nonlinear_cond.dat {0}".format(opt_cond_data))
        make_figure_opt_cond(
                file_data=opt_cond_data, 
                title="Nonlinear photoconductivity for hBN", 
                ylabel=r"$\sigma_{xxy}(\omega, -\omega)$",
                file_out=opt_cond_figure)

        # make_figure_opt_cond(file_data="optcond.dat", title="Optical conductivity", xlabel=r"$\hbar\omega$ (ev)",
# ylabel=r"$\sigma_{xx}(\omega)$", file_out="optcond.pdf"):


    print_title("Ran all the KITE-examples")


def clean():
    # Header
    print_title("Remove the previous calculations from the examples folder")
    terminal("rm *.h5 *.dat *.pdf")


if __name__ == "__main__":
    clean()
    main()
