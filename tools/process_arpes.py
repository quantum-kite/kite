"""       
        ##############################################################################      
        #                        KITE | Release  1.1                                 #      
        #                                                                            #      
        #                        Kite home: quantum-kite.com                         #           
        #                                                                            #      
        #  Developed by: Simao M. Joao, Joao V. Lopes, Tatiana G. Rappoport,         #       
        #  Misa Andelkovic, Lucian Covaci, Aires Ferreira, 2018-2022                 #      
        #                                                                            #      
        ##############################################################################      
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib import rc

rc('font',size=30)
rc('axes',labelsize=30,linewidth=4)

def process_arpes(name):
    # open file and read from it
    f = open(name,"r")
    whole = f.readlines()


    # inside that file, find where is each field. There are three fields: 
    # list of energies, list of k-vectors and the arpes matrix
    print("Processing file: finding location of each dataset. ", end="")
    cutarpes, cutenergies, cutvectors = -1, -1, -1
    count = 0
    found = 0
    for i in whole:
      if("ARPES:" in i):
        cutarpes = count + 1
        found += 1
      elif("Energies:" in i):
        cutenergies = count + 1
        found += 1
      elif("k-vectors:" in i):
        cutvectors = count + 1
        found += 1

      count +=1
      if(found == 3):
        break
    print(f"Done. Locations: {cutarpes} {cutenergies} {cutvectors}")

    # process the energies and arpes matrix into numpy arrays
    print("Fetching energies. ", end="")
    energies = np.array([float(i[:-1]) for i in whole[cutenergies:cutarpes - 1]])
    print(f"Done. Number of energies: {len(energies)}")

    # process the arpes matrix
    print("Fetching ARPES matrix. ", end="")
    arpes = np.array([[float(y) for y in i[:-1].split(" ") if y!=""] for i in whole[cutarpes:]])
    print(f"Done. Arpes matrix dimensions: {len(arpes)} {len(arpes[0])}")

    # find the number of k-vectors
    print("Fetching k vectors. ", end="")
    kvectors = whole[cutvectors:cutenergies - 1]
    num_vectors = len(kvectors)
    nvectors = np.linspace(1, num_vectors, num_vectors)
    print(f"Done. Number of k vectors: {num_vectors}")
    # print(nvectors)

    # plot the matrix
    print("Making meshgrid. ", end="")
    X,Y = np.meshgrid(nvectors, energies)
    print("Done.")

    print("Plotting 2D color map. ", end="")
    fig, axs = plt.subplots(1,1,figsize=(20,15))
    axs.pcolor(X,Y,arpes,cmap='hot', rasterized=True)
    axs.set_xlabel("k-points", fontsize=30)
    axs.set_ylabel("Energy", fontsize=30)
    print("Done.")

    print("Saving figure.", end="")
    plt.savefig("arpes.pdf", format="pdf", dpi=50)
    print("Done.")

if __name__ == "__main__":
    name = sys.argv[1]
    process_arpes(name)
