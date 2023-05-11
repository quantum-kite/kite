import numpy as np
import subprocess
import sys
sys.path.insert(1, '..')
import compare

# Parameters
file1 = "config.h5"
file2 = "configREF.h5"
file3 = "dos.dat"
file4 = "dosREF.dat"
dset1 = dset2 = "/Calculation/dos/MU"
tol = 1e-8

commands = ["SEED=3 ../KITEx config.h5", "../KITE-tools config.h5"]

for command in commands:
    result = subprocess.run([command], capture_output=True, shell=True)
    # print("stdout:", result.stdout)
    # print("stderr:", result.stderr)
    error_code = result.returncode
    if error_code != 0:
        print("ERROR")
        exit(0)


# Perform the first comparison
argv = [file1,dset1,file2,dset2]
res = compare.compare(argv)

if res[0] > tol:
    print("Problem", end=" ")
else:
    print("OK", end=" ")


# Perform the second comparison
argv = [file3,file4]
res = compare.compare_txt(argv)

if res[0] > tol:
    print("Problem")
else:
    print("OK")
