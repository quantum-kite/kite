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
