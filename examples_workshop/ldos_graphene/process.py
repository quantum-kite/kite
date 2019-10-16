"""       
        ##############################################################################      
        #                        KITE | Release  1.0                                 #      
        #                                                                            #      
        #                        Kite home: quantum-kite.com                         #           
        #                                                                            #      
        #  Developed by: Simao M. Joao, Joao V. Lopes, Tatiana G. Rappoport,         #       
        #  Misa Andelkovic, Lucian Covaci, Aires Ferreira, 2018-2019                 #      
        #                                                                            #      
        ##############################################################################      
"""
import numpy as np
import matplotlib.pyplot
import sys

name = sys.argv[1]
f = open(name,"r")
a = f.readlines()
f.close()
N = len(a)
theta = np.pi/3.0
a1 = np.array([1 + np.cos(theta),np.sin(theta)])
a2 = np.array([0, 2*np.sin(theta)])
delta = np.array([1,0])

xs = []
ys = []
values = []

for i in range(N):
    line = a[i].split(" ") 
    i0 = int(line[0])
    i1 = int(line[1])
    orb = int(line[2])
    value = float(line[3])
    r = a1*i0 + a2*i1 + orb*delta
    xs.append(r[0])
    ys.append(r[1])
    values.append(value)

M = max(values)
m = min(values)
colors=[]
for i in values:
    v = (i-m)/(M-m)
    colors.append(v)


F = 6
S = 70
f1 = F*(np.cos(theta*0.5))
f2 = F*(np.sin(theta*0.5) + 1)
matplotlib.pyplot.figure(figsize=(f1,f2))
matplotlib.pyplot.scatter(xs, ys, c=colors, s=S, cmap="hot",vmin=0.0, vmax=1.0)
matplotlib.pyplot.savefig(name[:-4]+".png")
