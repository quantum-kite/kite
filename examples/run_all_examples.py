""" Run all the examples in the kite/examples-folder

    ##########################################################################
    #                         Copyright 2022, KITE                           #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Last updated: 18/07/2022
"""

import numpy as np
import matplotlib.pyplot as plt
from os import system as terminal

_kitex_dir = "../build"
_kite_tools_dir = "../tools/build"


def run_calculation(input_file="output.h5"):
    """Run KITEx"""
    print_command("- - - -             Doing the KITEx-calculation                - - - - -")
    terminal("{0}/KITEx {1}".format(_kitex_dir, input_file))


def run_tools(input_file="output.h5"):
    """Run KITE-tools to obtain the DOS"""
    print_command("- - - -             Doing the KITE-tools postprocessing        - - - - -")
    terminal("{0}/KITE-tools {1}".format(_kite_tools_dir, input_file))


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


def make_figure_cond_cd(file_data="condDC.dat", title="DC conductivity", xlabel="E (eV)",
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


def print_command(text=""):
    print("\033[94m {0} \033[0m".format(text))


def print_title(text=""):
    print_command("-" * 74 + "\n" + " " * 10 + text + "\n" + "-" * 74 + "\n")


def main():
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

    # Example 1: dos_square_lattice.py
    print_command("======= Example 1: DOS for a square lattice                      =========")
    import dos_square_lattice
    print_command("- - - -            Making the configuration file                 - - - - -")
    hdf5_file = dos_square_lattice.main()
    pre_file_name = hdf5_file.replace("-output.h5", "")
    dos_data = "{0}-dos.dat".format(pre_file_name)
    dos_figure = "{0}-dos.pdf".format(pre_file_name)
    run_calculation(hdf5_file)
    run_tools(hdf5_file)
    terminal("mv dos.dat {0}".format(dos_data))
    make_figure_dos(file_data=dos_data, title="DOS Square lattice", xlabel="Energy (|t|)",
                    file_out=dos_figure)

    # Example 2: dos_square_lattice_twisted_bc.py
    print_command("    ======= Example 2: DOS for a square lattice with twisted BC  =========")
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

    # Example 3: dos_checkerboard_lattice.py
    print_command("======= Example 3: DOS for a checkerboard lattice                =========")
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

    # Example 4: dos_graphene.py
    print_command("======= Example 4: DOS for graphene                              =========")
    import dos_graphene as example
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

    # Example 5: dos_cubic_lattice_twisted_bc.py
    print_command("======= Example 5: DOS for a cubic lattice                       =========")
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

    # Example 6: basic_on_site_disorder.py
    print_command("======= Example 6: DOS for lattice with onsite disorder          =========")
    import basic_on_site_disorder as example
    print_command("- - - -            Making the configuration file                 - - - - -")
    hdf5_file = example.main()
    pre_file_name = hdf5_file.replace("-output.h5", "")
    dos_data = "{0}-dos.dat".format(pre_file_name)
    dos_figure = "{0}-dos.pdf".format(pre_file_name)
    run_calculation(hdf5_file)
    run_tools(hdf5_file)
    terminal("mv dos.dat {0}".format(dos_data))
    make_figure_dos(file_data=dos_data, title="DOS Graphene with on-site disorder",
                    file_out=dos_figure)

    # Example 7: basic_vacancies.py
    print_command("======= Example 7: DOS for lattice with vacancies                =========")
    import basic_vacancies as example
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

    # Example 8: basic_mixed_disorder.py
    print_command("======= Example 8: DOS for lattice with vacancies and onsite disorder ====")
    import basic_mixed_disorder as example
    print_command("- - - -            Making the configuration file                 - - - - -")
    hdf5_file = example.main()
    pre_file_name = hdf5_file.replace("-output.h5", "")
    dos_data = "{0}-dos.dat".format(pre_file_name)
    dos_figure = "{0}-dos.pdf".format(pre_file_name)
    run_calculation(hdf5_file)
    run_tools(hdf5_file)
    terminal("mv dos.dat {0}".format(dos_data))
    make_figure_dos(file_data=dos_data, title="DOS Graphene with vacancies & on-site disorder",
                    file_out=dos_figure)

    # Example 9: optcond_gaussian_disorder.py
    print_command("======= Example 9: Optical conductivity with onsite disorder     =========")
    import optcond_gaussian_disorder as example
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

    # Example 10: haldane.py
    print_command("======= Example 10: DOS & optical conductivity for the Haldene model =====")
    import haldane as example
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
    make_figure_dos(file_data=dos_data, title="DOS Graphene Gaussian On-site disorder",
                    file_out=dos_figure)
    make_figure_cond_cd(file_data=cond_dc_data, title="DC Conductivity Haldene",
                        file_out=cond_dc_figure)

    print_command("-" * 74 + "\n" + " " * 10 + "Ran all the examples in the kite/examples-folder\n" + "-" * 74 + "\n")

    # Example 11: phosphorene.py
    print_command("======= Example 11: DC coductivity for phosphorene in XX         =========")
    import phosphorene as example
    print_command("- - - -            Making the configuration file                 - - - - -")
    hdf5_file = example.main()
    pre_file_name = hdf5_file.replace("-output.h5", "")
    cond_dc_data = "{0}-condDC.dat".format(pre_file_name)
    cond_dc_figure = "{0}-condDC.pdf".format(pre_file_name)
    run_calculation(hdf5_file)
    example.post_process(hdf5_file)
    terminal("mv condDC.dat {0}".format(cond_dc_data))
    make_figure_cond_cd(file_data=cond_dc_data, title="DC Conductivity Phosphorene XX",
                        file_out=cond_dc_figure)

    # Example 12: phosphorene.py
    print_command("======= Example 12: DC coductivity for phosphorene in YY         =========")
    import phosphorene as example
    print_command("- - - -            Making the configuration file                 - - - - -")
    hdf5_file = example.main(direction='yy')
    pre_file_name = hdf5_file.replace("-output.h5", "")
    cond_dc_data = "{0}-condDC.dat".format(pre_file_name)
    cond_dc_figure = "{0}-condDC.pdf".format(pre_file_name)
    run_calculation(hdf5_file)
    example.post_process(hdf5_file)
    terminal("mv condDC.dat {0}".format(cond_dc_data))
    make_figure_cond_cd(file_data=cond_dc_data, title="DC Conductivity Phosphorene YY",
                        file_out=cond_dc_figure)

    # Example 13: twisted_bilayer.py
    print_command("======= Example 12: DOS for twisted bilayer graphene at 21.787 degrees ===")
    import twisted_bilayer as example
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

    print_title("Ran all the KITE-examples")


def clean():
    # Header
    print_title("Remove the previous calculations from the examples folder")
    terminal("rm *.h5 *.dat *.pdf")


if __name__ == "__main__":
    clean()
    main()
