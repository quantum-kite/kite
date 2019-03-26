import numpy as np
# np.set_printoptions(precision=2)

a1 = np.array([1,0])*np.pi*2
a2 = np.array([0,1])*np.pi*2
scale = 4.1

def eps(k):
  return -2.0*(np.cos(np.dot(k,a1)) + np.cos(np.dot(k,a2)))/scale

def cheb(n,x):
  return np.cos(n*np.arccos(x))

for p in [(0,0),(0.125,0.125),(0.25,0.25),(0.375,0.375),(0.5,0.5)]:
  momentum = np.array(p)
  energy = eps(np.array(p))
  # print(str(p) + " " + str(energy),end="")
  # print("momentum: " + str(momentum))
  # print("energy: " + str(energy))
  for n in range(5):
    moment = cheb(n, energy)
    print("{0:0.4f} ".format(moment), end="")
  print("")
