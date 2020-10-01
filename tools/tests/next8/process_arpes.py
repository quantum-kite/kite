import numpy as np
import sys
from matplotlib import rc


rc('font',size=50)
rc('axes',labelsize=50,linewidth=4)

# open file and read from it
name = sys.argv[1]
f = open(name,"r")
whole = f.readlines()


# inside that file, find where is each field. There are three fields: 
# list of energies, list of k-vectors and the arpes matrix
print("Processing file: finding location of each dataset. ", end="")
cutarpes, cutenergies, cutvectors = -1, -1, -1
count = 0
found = 0
for i in whole:
  if(i == "ARPES:\n"):
    cutarpes = count + 1
    found += 1
  elif(i == "Energies:\n"):
    cutenergies = count + 1
    found += 1
  elif(i == "k-vectors:\n"):
    cutvectors = count + 1
    found += 1

  count +=1
  if(found == 3):
    break
print("Done.")
print(cutarpes, cutenergies, cutvectors)

# process the energies and arpes matrix into numpy arrays
print("Fetching energies. ", end="")
energies = np.array([float(i[:-1]) for i in whole[cutenergies:cutarpes - 1]])
print("Number of energies: ", len(energies))
print("Done.")
# A = [[float(y) for y in i[:-1].split(" ") if y!=""] for i in whole[cutarpes:]]
# for a in A:
  # print(len(a))
print("Fetching ARPES matrix. ", end="")
arpes = np.array([[float(y) for y in i[:-1].split(" ") if y!=""] for i in whole[cutarpes:]])
print("Arpes matrix dimensions:", len(arpes), len(arpes[0]), "\n")
print("Done.")
print(energies)
print("arpes ")
print(np.max(arpes))
print("done ")

# find the number of k-vectors
print("Fetching k vectors. ", end="")
kvectors = whole[cutvectors:cutenergies - 1]
num_vectors = len(kvectors)
nvectors = np.linspace(1, num_vectors, num_vectors)
print("Number of k vectors: ", num_vectors,"\n")
print(kvectors)
kx = np.zeros(num_vectors) 
ky = np.zeros(num_vectors) 
for i,k in enumerate(kvectors):
    kk = k[:-1].split(" ")
    kx[i] = float(kk[0])
    ky[i] = float(kk[1])
    np.savetxt("dos_k" + str(kx[i]) + ".dat", np.c_[energies,arpes[:,i]])

print(kx)
print(ky)

print("Done.")
