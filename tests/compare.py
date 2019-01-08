import h5py
import sys
import numpy as np

f1 = h5py.File(sys.argv[1], 'r')
f2 = h5py.File(sys.argv[3], 'r')
dset1 = f1[sys.argv[2]]
dset2 = f2[sys.argv[4]]
npset1 = np.array(dset1[:])
npset2 = np.array(dset2[:])
set_dif = np.absolute(npset2 - npset1)

# sum of the differences
sumall = np.sqrt((set_dif**2).sum())

# maximum difference between the same elements of both arrays
maxim  = np.amax(set_dif)

norm1 = np.linalg.norm(npset1)
norm2 = np.linalg.norm(npset2)
pct = sumall/np.sqrt(norm1*norm2)

print("{:<11f} {:<11f} {:<11f} {:<11f} {:<11f}".format(sumall, maxim, norm1, norm2, pct))
