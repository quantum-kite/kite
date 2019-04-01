import numpy as np
import matplotlib.pyplot as plt
import sys

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
  elif(i == "Energies: \n"):
    cutenergies = count + 1
    found += 1
  elif(i == "k-vectors: \n"):
    cutvectors = count + 1
    found += 1

  count +=1
  if(found == 3):
    break
print("Done.")
# print(cutarpes, cutenergies, cutvectors)

# process the energies and arpes matrix into numpy arrays
print("Fetching energies. ", end="")
energies = np.array([float(i[:-1]) for i in whole[cutenergies:cutarpes - 1]])
print("Done.")
# A = [[float(y) for y in i[:-1].split(" ") if y!=""] for i in whole[cutarpes:]]
# for a in A:
  # print(len(a))
print("Fetching ARPES matrix. ", end="")
arpes = np.array([[float(y) for y in i[:-1].split(" ") if y!=""] for i in whole[cutarpes:]])
print("Done.")
# print(energies)
# print(arpes)

# find the number of k-vectors
print("Fetching k vectors. ", end="")
kvectors = whole[cutvectors:cutenergies - 1]
num_vectors = len(kvectors)
nvectors = np.linspace(1, num_vectors, num_vectors)
print("Done.")
# print(nvectors)

# plot the matrix
print("Making meshgrid. ", end="")
X,Y = np.meshgrid(nvectors, energies)
print("Done.")
print("Plotting 2D color map. ", end="")
plt.pcolor(X,Y,arpes)
print("Done.")
print("Saving figure.", end="")
plt.savefig("arpes.png")
print("Done.")
