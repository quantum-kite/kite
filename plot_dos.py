import matplotlib.pyplot as plt
import h5py
import matplotlib as mpl
import numpy as np
#************************************
#keep these definitions for kite website
import seaborn as sns
mpl.rcParams['figure.dpi'] = 100
mpl.rcParams['savefig.dpi'] = 100
sns.set_style("white")
#Kite color scheme
colors = ["dusty purple", "faded green","windows blue", "amber", "greyish"]
current_palette=sns.xkcd_palette(colors)
sns.set_palette(current_palette)
sns.set_style("ticks")
sns.set_context("talk",font_scale=1.1)
#*************************************

data=np.loadtxt('dos.dat')
#data2=np.loadtxt('condphxx.txt')
#data3=np.loadtxt('condphyy.txt')
#plt.title('Phosphorene')
#plt.plot(data2[:,0]+(7.0-6.16)/2,data2[:,1],label='$\sigma_{xx}$')
#plt.plot(data3[:,0]+(7.0-6.16)/2,data3[:,1],label='$\sigma_{yy}$')
plt.plot(data[:,0],data[:,1])
plt.xlabel('E (eV)')
plt.ylabel('DOS (a.u)')
plt.ylim(ymin=0)
#plt.xlim(-3.5,3.5)  
plt.legend()
plt.tight_layout()
plt.show()

