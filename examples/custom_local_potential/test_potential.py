#!/usr/bin/env python
# coding: utf-8

# # Test potential
# 
# This is a short script to quickly test the local custom potential features in KITE

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import sys


# In[ ]:


# Check if running inside a jupyter notebook
if 'ipykernel' in sys.modules:
    get_ipython().system('jupyter-nbconvert --to script test_potential.ipynb ')


# In[ ]:


# Info about the files
DIM = 2 # dimension

# Find all the potential files
files = os.listdir()
pots = []
for filename in files:
    if "local_potential" in filename:
        pots.append(filename)
        
# Extract the data from the files
data = {}
for pot in pots:
    print("Extracting from:", pot)
    
    with open(pot, 'r') as f:
        dat = f.readlines()        
        
        for line in dat:
            # orb = line[DIM]
            x,y,orb,V = line[:-1].split(" ")
                
            x = float(x)
            y = float(y)
            orb = int(orb)
            V = float(V)      
            
            if orb not in data.keys():
                data[orb] = []
            else:
                data[orb].append([x,y,V])
    
# Number of orbitals
Norbs = len(data.keys())
print("Number of orbitals:", Norbs)

for key in data.keys():
    data[key] = np.array(data[key])


# In[ ]:


# Represent the potential for each orbital
for o in range(Norbs):
    x = data[o][:,0]
    y = data[o][:,1]
    V = data[o][:,2]

    fig, axs = plt.subplots(1,1,figsize=(10,8))
    axs.set_title("Local potential at orbital " + str(o), fontsize=20)
    scat = axs.scatter(x,y,c = V, cmap="Spectral", s=0.5)
    plt.colorbar(scat)
    plt.show()


# In[ ]:





# In[ ]:





# In[ ]:




