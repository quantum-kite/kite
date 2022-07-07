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
from numpy import linalg as LA
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys

name = sys.argv[1]
f = open(name,"r")
a = f.readlines()
f.close()
N = len(a)
theta = np.pi/3.0
a1 = np.array([1 + np.cos(theta),np.sin(theta)])
a2 = np.array([0, 2*np.sin(theta)])

xmin, xmax, Nx = 40, 55, 300
ymin, ymax, Ny = 75, 90, 300
gridxs = np.linspace(xmin,xmax,Nx)
gridys = np.linspace(ymin,ymax,Ny)
meshx, meshy = np.meshgrid(gridxs, gridys)

def triangle(x):
    return np.exp(x)

def process_orb(ORB):
    count = 0
    xs = []
    ys = []
    non_normalized_values = []

    for i in range(N):
        line = a[i].split(" ") 
        i0 = int(line[0])
        i1 = int(line[1])
        orb = int(line[2])
        value = float(line[3])
        if(orb == ORB):
            r = a1*i0 + a2*i1
            xs.append(r[0])
            ys.append(r[1])
            non_normalized_values.append(value)
            count += 1

    # Now we have to find the points that are inside this window
    total  = np.zeros([Ny, Nx])
    totalu = np.zeros([Ny, Nx])
    MAX = max(non_normalized_values)*2.0
    for i in range(count):
        print(i,count)
        x, y, variance, amplitudeu = xs[i], ys[i], non_normalized_values[i]/MAX, non_normalized_values[i]
        x, y, variance, amplitudeu = xs[i], ys[i], 0.2, non_normalized_values[i]
        if xmin <= x <= xmax and ymin <= y <= ymax:
            print(x,y)
            exponent = -((meshx - x)**2 + (meshy - y)**2)/variance
            totalu += triangle(exponent)*amplitudeu

    return totalu

unnormalized_orb0 = process_orb(0)
unnormalized_orb1 = process_orb(1)
unnormalized_orb2 = process_orb(2)
unnormalized = unnormalized_orb0 + unnormalized_orb1 + unnormalized_orb2

cdict1 = {'blue':   ((0.0, 0.0, 0.0),
                   (0.3, 0.5, 0.5),
                   (1.0, 1.0, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'red':  ((0.0, 0.0, 0.0),
                   # (0.5, 0.0, 0.0),
                   # (0.6, 0.8, 0.8),
                   (1.0, 0.0, 0.0))
         }

blue_red1 = mpl.colors.LinearSegmentedColormap('BlueRed1', cdict1)

ma = blue_red1
plt.figure(figsize=(6,6))
m = 0
M = 0.1
# plt.contourf(gridxs, gridys, unnormalized_orb0, 100, cmap = ma, vmin = m, vmax = M)
plt.contourf(gridxs, gridys, unnormalized_orb0, 100)
# plt.contourf(gridxs, gridys, np.log(unnormalized_orb0+0.001), 100)
# plt.contourf(gridxs, gridys, unnormalized_orb0, 100, cmap = ma)
plt.colorbar()
plt.savefig("tmd_ldos_u0.png")

plt.figure(figsize=(6,6))
# plt.contourf(gridxs, gridys, unnormalized_orb1, 100, cmap = ma, vmin = m, vmax = M)
# plt.contourf(gridxs, gridys, unnormalized_orb1, 100, cmap = ma)
plt.contourf(gridxs, gridys, unnormalized_orb1, 100)
# plt.contourf(gridxs, gridys, np.log(unnormalized_orb1+0.001), 100)
plt.colorbar()
plt.savefig("tmd_ldos_u1.png")

plt.figure(figsize=(6,6))
# plt.contourf(gridxs, gridys, unnormalized_orb2, 100, cmap = ma, vmin = m, vmax = M)
# plt.contourf(gridxs, gridys, np.log(unnormalized_orb2+0.001), 100)
plt.contourf(gridxs, gridys, unnormalized_orb2, 100)
# plt.contourf(gridxs, gridys, unnormalized_orb2, 100, cmap = ma)
plt.colorbar()
plt.savefig("tmd_ldos_u2.png")

plt.figure(figsize=(6,6))
# plt.contourf(gridxs, gridys, unnormalized,  100, cmap = ma, vmin = m, vmax = M)
# plt.contourf(gridxs, gridys, unnormalized,  100, cmap = ma)
plt.contourf(gridxs, gridys, unnormalized,  100)
# plt.contourf(gridxs, gridys, np.log(unnormalized+0.001),  100)
plt.colorbar()
plt.savefig("tmd_ldos_u.png")
