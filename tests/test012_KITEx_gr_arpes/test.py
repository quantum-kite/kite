import numpy as np
import subprocess
import sys
sys.path.insert(1, '..')
import compare

# Parameters
file1 = "config.h5"
file2 = "configREF.h5"
dset1 = dset2 = "/Calculation/arpes/kMU"
tol = 1e-8

result = subprocess.run(["SEED=3 ../KITEx config.h5"], capture_output=True, shell=True)
# print("stdout:", result.stdout)
# print("stderr:", result.stderr)
error_code = result.returncode
if error_code != 0:
    print("ERROR")
    exit(0)

# Perform the comparison
argv = [file1,dset1,file2,dset2]
res = compare.compare(argv)

if res[0] > tol:
    print("Problem")
else:
    print("OK")
