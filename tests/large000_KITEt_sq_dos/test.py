import numpy as np
import subprocess
import sys
sys.path.insert(1, '..')
import compare

# Parameters
file1 = "dos.dat"
file2 = "dosREF.dat"
tol = 1e-8

h5file = "../large000_KITEx_sq_dos/config.h5"
result = subprocess.run(["../KITE-tools " + h5file], capture_output=True, shell=True)
# print("stdout:", result.stdout)
# print("stderr:", result.stderr)
error_code = result.returncode
if error_code != 0:
    print("ERROR")
    exit(0)

# Perform the comparison
argv = [file1, file2]
res = compare.compare_txt(argv)

if res[0] > tol:
    print("Problem")
else:
    print("OK")
