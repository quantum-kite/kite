import h5py
import sys
import numpy as np

def compare(argv):

    f1 = h5py.File(argv[0], 'r')
    f2 = h5py.File(argv[2], 'r')
    dset1 = f1[argv[1]]
    dset2 = f2[argv[3]]

    npset1 = np.array(dset1[:])
    npset2 = np.array(dset2[:])
    set_dif = np.absolute(npset2 - npset1)

    # sum of the differences
    sumall = np.sqrt((set_dif**2).sum())

    # maximum difference between the same elements of both arrays
    maxim  = np.amax(set_dif)

    norm1 = np.linalg.norm(npset1)
    norm2 = np.linalg.norm(npset2)

    # relative error
    # pct = sumall/np.sqrt(norm1*norm2)*0
    pct = 0

    return sumall, maxim, norm1, norm2, pct
    # print("{:<11f} {:<11f} {:<11f} {:<11f} {:<11f}".format(sumall, maxim, norm1, norm2, pct))

def compare_txt(argv):
    # Compare text files

    f1 = argv[0]
    f2 = argv[1]

    # Load datasets
    d1 = np.loadtxt(f1)
    d2 = np.loadtxt(f2)


    # sum of the differences
    set_dif = np.absolute(d1 - d2)
    sumall = np.sqrt((set_dif**2).sum())

    # maximum difference between the same elements of both arrays
    maxim  = np.amax(set_dif)

    norm1 = np.linalg.norm(d1)
    norm2 = np.linalg.norm(d2)

    # relative error
    pct = 0

    return sumall, maxim, norm1, norm2, pct

# if __name__ == "__main__":
