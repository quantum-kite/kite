import numpy as np
import matplotlib.pyplot as plt
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
  if(i == "ARPES:\n" or i == "ARPES: \n"):
    cutarpes = count + 1
    found += 1
  elif(i == "Energies:\n" or i == "Energies: \n"):
    cutenergies = count + 1
    found += 1
  elif(i == "k-vectors:\n" or i == "k-vectors: \n"):
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
print("Done.")
print(nvectors)

# plot the matrix
print("Making meshgrid. ", end="")
X,Y = np.meshgrid(nvectors, energies)
print("Done.")
print("Plotting 2D color map. ", end="")

plt.figure(figsize=(30,20))
plt.pcolormesh(X,Y,arpes,cmap='hot',vmin=0,vmax=300000)
n1 = 568
n2 = 284
xcoords = [ (n1-1)*6, (n1-1+n2-1)*6]

# for xc in xcoords:
    # plt.axvline(x=xc,ymin=-10, ymax=10,lw=5,color='w')
#plt.xlabel(r'K', labelpad=1, size=60)
plt.ylabel(r'E(eV)', labelpad=1, size=60)
# plt.annotate(r'K',xy=(0.565,-0.1),xycoords='axes fraction',size=70)
# plt.annotate(r'$\Gamma$',xy=(0.0,-0.1),xycoords='axes fraction',size=70)
# plt.annotate(r'$\Gamma$',xy=(0.995,-0.1),xycoords='axes fraction',size=70)
# plt.annotate(r'M',xy=(0.35,-0.1),xycoords='axes fraction',size=70)
yticks = np.linspace(-9,9,7)
plt.yticks(yticks)
plt.xticks([])
print("Done.")
print("Saving figure.", end="")
plt.savefig("arpes.png")
print("Done.")
