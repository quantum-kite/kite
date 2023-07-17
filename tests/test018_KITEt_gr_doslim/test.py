import numpy as np
import subprocess
import sys
sys.path.insert(1, '..')
import compare

# Parameters
tol = 1e-8
N1 = 51
N2 = N1*2-1
d = (N1-1)//2

result1 = subprocess.run([f"../KITE-tools config.h5 --DOS -E -0.1 0.1 {N1} -N dos1.dat -X"], capture_output=True, shell=True)
result2 = subprocess.run([f"../KITE-tools config.h5 --DOS -E -0.2 0.2 {N2} -N dos2.dat -X"], capture_output=True, shell=True)
# print("stdout:", result.stdout)
# print("stderr:", result.stderr)
error_code = result1.returncode*result2.returncode
if error_code != 0:
    print("ERROR")
    exit(0)

# Perform the comparison
# argv = [file1,dset1,file2,dset2]
# res = compare.compare(argv)
dos1 = np.loadtxt("dos1.dat")[:,1]
dos2 = np.loadtxt("dos2.dat")[d:d+N1,1]
res = np.linalg.norm(dos2 - dos1)

if res > tol:
    print("Problem")
else:
    print("OK")
