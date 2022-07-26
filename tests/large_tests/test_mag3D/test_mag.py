# Compare the density of states calculated with KITE (cubic lattice
# with applied magnetic field along y) with the 3DEG with Landau 
# levels

import matplotlib.pyplot as plt
import numpy as np

h = 6.62607015e-34 # J/s
e = 1.60217662e-19 # Coulomb
a = 1e-9  # meters (typical length scale)

N = 64
Ny = 32
mult = 1

t = 1
T = h/e/a/a 
print("Magnetic field in Tesla", T/N)

# spacings
de1 = 4*np.pi/N*mult
de2 = 2*np.pi/Ny
print(de1, de2)

# calculate the Landau levels
NI = 10
NJ = N
landaus = []
colors = []
c = ['black', 'red', 'blue', 'green', 'cyan','orange']*50
for i in range(NI):
    for j in range(NJ):
        HO = -4*t + de1*(0.5 + i) # harmonic oscillator part
        E = HO - 2*np.cos(de2*j)
        landaus.append(E) 
        colors.append(c[j])

# fetch DoS calculated with KITE
dos = np.loadtxt("dos.dat")

# check if the Landau levels are in the predicted spot
fig, axs = plt.subplots(1,1, figsize=(10,5))
axs.plot(dos[:,0], dos[:,1], label="KITEx")
axs.set_xlim([-6,-5.6])
axs.set_ylim([-0.001, 0.04])
for i,c in zip( landaus, colors):
    axs.axvline(i,c=c, linestyle='-')
axs.legend(fontsize=15)
axs.set_xlabel("Energy (t)", fontsize=15)
axs.set_ylabel("DoS (1/t)", fontsize=15)
plt.show()
