import h5py
import sys
import numpy as np

f = h5py.File(sys.argv[1], 'r')
dset = f[sys.argv[2]]
npset = np.array(dset[:])
sumall = np.sqrt((npset**2).sum())
maxim  = np.sqrt(np.amax(npset**2))

print(sumall, maxim)
