""" Script to Extract Data for the longitudinal DC-conductitivy obtained by the single-shot
    method (note that this method does not require post-processing by KITE-tools)

    ##########################################################################
    #                         Copyright 2022, KITE                           #
    #                         Home page: quantum-kite.com                    #
    ##########################################################################

    Last updated: 28/07/2022
"""

__all__ = ["post_process_singleshot_conductivity_dc"]

import numpy as np
import h5py as hdf
import sys


def post_process_singleshot_conductivity_dc(filename="output.h5"):
    """Read the HDF5-file and process the singlesho_conductivity_dc"""

    # open the HDF5-file and read the data
    file = hdf.File(filename, 'r+')
    data = np.array(file['Calculation/singleshot_conductivity_dc/SingleShot'], dtype=np.float64)
    file.close()

    # export the file
    output = (filename[:-3] if filename[-3:] == ".h5" else filename) + ".dat"
    data2 = np.array([data[:, 0], data[:, 2], data[:, -1]]).transpose()
    np.savetxt(output, data2)
    print("Done")


if __name__ == "__main__":
    filename = str(sys.argv[1])
    post_process_singleshot_conductivity_dc(filename=filename)
