""" Script to Extract Data for the longitudinal DC-conductitivy obtained by the single-shot
    method (note that this method does not require post-processing by KITE-tools)

    ##########################################################################
    #                     Copyright 2022, KITE                               #
    #                 Home page: quantum-kite.com                            #
    ##########################################################################
    Last updated: 27/07/2022
"""


import h5py as hdf
import numpy as np
import sys

filename = str(sys.argv[1])

file = hdf.File(filename, 'r+')
data = np.array(file['Calculation/singleshot_conductivity_dc/SingleShot']).astype(float)
file.close()

filename = filename[:-2] + "dat"
data2 = np.array([data[:,0],data[:,2],data[:,-1]]).transpose()
np.savetxt(filename,data2)
print("Done")
