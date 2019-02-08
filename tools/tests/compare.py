import sys
import numpy as np

f1 = open(sys.argv[1], "r")
f2 = open(sys.argv[2], "r")

a1 = f1.readlines()
a2 = f2.readlines()

b1 = np.zeros([len(a1),3])
b2 = np.zeros([len(a2),3])

for i in range(len(a1)):
  line = a1[i]
  c1, c2, c3 = line[:-1].split(" ")
  d1 = float(c1)
  d2 = float(c2)
  d3 = float(c3)

  b1[i][0] = d1
  b1[i][1] = d2
  b1[i][2] = d3



for i in range(len(a2)):
  line = a2[i]
  c1, c2, c3 = line[:-1].split(" ")
  d1 = float(c1)
  d2 = float(c2)
  d3 = float(c3)

  b2[i][0] = d1
  b2[i][1] = d2
  b2[i][2] = d3


dif=np.absolute(b1-b2)
M = np.amax(dif)
N = np.linalg.norm(dif)
print(M, N)

