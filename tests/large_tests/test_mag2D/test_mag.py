# Compare the DoS calculated with KITE (2D square lattice with 
# applied magneteic field along the z direction) with the DoS 
# of the 2DEG with Landau levels

import matplotlib.pyplot as plt
import numpy as np

h = 6.62607015e-34 # J/s
e = 1.60217662e-19 # Coulomb
a = 1e-9  # meters (typical length scale)

N = 512
mult = 1

t = 1
T = h/e/a/a 
print("Tesla", T/N)

# spacings
de1 = 4*np.pi/N*mult
print(de1)

# band edge inand 2D
E0 = -4*t

# calculate the Landau levels
NI = 10
landaus = []
for i in range(NI):    
    E = E0 + de1*(0.5 + i)
    landaus.append(E) 
landaus.sort()

# fetch DoS calculated with KITE
basedir = ""
file = basedir + "dos.dat"
dos = np.loadtxt(file)


# check if the Landau levels are in the predicted spot
fig, axs = plt.subplots(1,1, figsize=(10,5))
axs.plot(dos[:,0], dos[:,1], label="KITEx")
axs.set_xlim([-4,-3.8])
for landau in landaus:
    axs.axvline(landau,c='k', linestyle='--')
axs.legend(fontsize=15)
axs.set_xlabel("Energy (t)", fontsize=15)
axs.set_xlabel("DoS (1/t)", fontsize=15)
plt.show()
