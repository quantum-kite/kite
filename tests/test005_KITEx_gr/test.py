import numpy as np
import subprocess
import sys
sys.path.insert(1, '..')
import compare

# Parameters
file1 = "config.h5"
file2 = "configREF.h5"
tol = 1e-8

dsets = ["/Calculation/conductivity_optical/Gammayy", "/Calculation/conductivity_optical/Lambdayy"]
for dset in dsets:

    dset1 = dset2 = dset
    result = subprocess.run(["SEED=ones ../KITEx config.h5"], capture_output=True, shell=True)
    error_code = result.returncode
    if error_code != 0:
        print("ERROR")
        exit(0)

    # Perform the comparison
    argv = [file1,dset1,file2,dset2]
    res = compare.compare(argv)

    if res[0] > tol:
        print("Problem", *res)
        exit(0)
    else:
        print("OK ", end="")
print("")




