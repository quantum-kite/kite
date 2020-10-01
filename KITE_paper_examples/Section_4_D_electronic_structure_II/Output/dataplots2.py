#!/usr/bin/env python3

# This is a comment
"""	
This is a Python Script created to plot 2D density plot given a x,y,z file. 
It was created by Jose Hugo Garcia  Aguilar during his post-doctoral fellowship,
 at the ICTP-SAIFR institue of IFT UNESP (2015-2016) and modified completely by Tatiana Rappoport in 2019 :)
	
"""
#IMPORTING THE NECESSARY MODULES
import argparse							#Used to parse parameters to the Python Script
import numpy as np						#Nummerical library
import matplotlib.pyplot as plt			#Plot library
from scipy.interpolate import griddata	#Used for griddata interpolation method

from matplotlib import rc
rc('font',size=20)
rc('axes',labelsize=20,linewidth=3)

parser = argparse.ArgumentParser(prog='color_plot_xyz.py', description='Make a (x,y,z) plot where (x,y) represents the coordinates in a plane z is the color.')
parser.add_argument('--filename' ,'-f',help='path to the xyz file', required=True)
parser.add_argument('--output','-o', help='Indicates the name of the outputfile', default='xyz.pdf')
parser.add_argument('--normalize','-N', help='Indicates if the data should be normalized', action='store_true' , default=False)
parser.add_argument('--normalize_inv','-NI', help='Indicates if the data should be normalized and inverted', action='store_true' , default=False)
parser.add_argument('--mean','-M', help='Indicates if the data should rescaled by the mean', action='store_true' , default=True)
parser.add_argument('--scale','-sc', help='Indicates the scale of the data, options are: lin, quad, log, exp ', default='lin')
parser.add_argument('--xplot','-X', help='Shows the data into de screen', action='store_true' , default=True )
parser.add_argument('--sample','-s',help=' Number of point in each direction of the grid ', default=500  , type=int)
parser.add_argument('--dinterpolation','-dinterpolation',help='Indicates the method for the interpolation of the data: linear, nearest, cubic', default='cubic')
parser.add_argument('--pinterpolation','-pinterpolation',help='Indicates the method for the interpolation of the plot points: none, nearest,bilinear, bicubic, spline16, spline36, hanning, hamming, hermite, kaiser, quadric, catrom, gaussian, bessel, mitchell, sinc, lanczos', default='bicubic')

#Get the parameters into the var namespace
var = parser.parse_args()

#Load the data
n1,n2,orb,ldos =np.loadtxt( var.filename , unpack=True )
ac=1.0
theta = np.pi/3.0
a1x=ac*(1 + np.cos(theta))
a2x=0
a1y=ac*np.sin(theta)
a2y=2*ac*np.sin(theta)
xmin, xmax= 40*ac, 55*ac
ymin, ymax= 75*ac, 90*ac


boolArr = (orb==0)
newArr = orb[boolArr]
x=a1x*n1[boolArr]+a2x*n2[boolArr]
y=a1y*n1[boolArr]+a2y*n2[boolArr]


x=x-xmin
y=y-ymin
Z0=ldos[(orb==0)]
Z1=ldos[(orb==1)]
Z2=ldos[(orb==2)]
Z=Z0+Z1+Z2
print(Z.max(), Z0.max(), Z1.max(), Z2.max())
#boolArr =((x>=xmin)&(x <= xmax)&(y <= ymax) & (y>=ymin))
#x=x[boolArr]
#y=y[boolArr]
#Z=Z[boolArr]
print(x,y)

#plt.scatter(nkx[::4], nky[::4], marker='.')
#print('Contents of the Bool Numpy Array : ', newArr)

if (var.mean):
	Z =  Z/Z.mean()
	Z0 =  Z0/Z.mean()
	Z1 =  Z1/Z.mean()
	Z2 =  Z2/Z.mean()


else:
	print ("Option in --scale is not defined, linear used")
	
#Multiply the var.sample parameter by 1j so that the integer part of its magnitude is interpreted as 
#specifying the number of points to create between the start and stop values, where the stop value is inclusive.

Z=np.log(Z)
Z0=np.log(Z0)
Z1=np.log(Z1)
Z2=np.log(Z2)
print(Z.max(), Z0.max(), Z1.max(), Z2.max())

var.sample=var.sample*1j
#Creates the grid
grid_x, grid_y = np.mgrid[0:(xmax-xmin):var.sample, 0:(ymax-ymin):var.sample]
#Creates the interpolation of the z value:
#						Point (X ,Y ,     Z ;  points to interp.	Used Method 
grid_z = griddata( np.transpose([x,y]) ,  Z, (grid_x, grid_y)	, method=var.dinterpolation )
grid_z0 = griddata( np.transpose([x,y]) , Z0, (grid_x, grid_y)	, method=var.dinterpolation )
grid_z1 = griddata( np.transpose([x,y]) , Z1, (grid_x, grid_y)	, method=var.dinterpolation )
grid_z2 = griddata( np.transpose([x,y]) , Z2, (grid_x, grid_y)	, method=var.dinterpolation )


#Plot the structure
plt.figure(figsize=(20,5))
ax1 = plt.subplot2grid((1,20),(0,0), rowspan=1,colspan=4)
ax2 = plt.subplot2grid((1,20),(0,5), rowspan=1,colspan=4)
ax3 = plt.subplot2grid((1,20),(0,10), rowspan=1,colspan=4)
ax4 = plt.subplot2grid((1,20),(0,15), rowspan=1,colspan=4)
zmax=Z2.max()
zmin=Z0.min()

ax1.imshow(grid_z.T,extent=(0,(xmax-xmin),0,(ymax-ymin)), origin='lower',interpolation=var.pinterpolation,cmap='jet' )
#ax1.xlabel(r'x(nm)', labelpad=5, size=50)
#ax1.ylabel(r'y(nm)',size=50)
ax2.imshow(grid_z0.T,extent=(0,(xmax-xmin),0,(ymax-ymin)), origin='lower',  interpolation=var.pinterpolation,cmap='jet' )
ax3.imshow(grid_z1.T,extent=(0,(xmax-xmin),0,(ymax-ymin)), origin='lower',  interpolation=var.pinterpolation ,cmap='jet')
ax4.imshow(grid_z2.T,extent=(0,(xmax-xmin),0,(ymax-ymin)), origin='lower', interpolation=var.pinterpolation ,cmap='jet')

#ax2.xlabel(r'x(nm)', labelpad=5, size=50)
#ax2.ylabel(r'y(nm)',size=50)
ax1.xticks = np.linspace(0,15,4)
#ax1.yticks = np.linspace(-0.8,0.1,4)
ax1.set_xticks(ax1.xticks)
#ax1.set_yticks(ax1.yticks)

#ax1.set_xlabel(r'x(nm)', labelpad=5, size=30)
plt.annotate(r'x(nm)',xy=(0.45,0.02),xycoords='figure fraction',size=30)

#yticks = np.linspace(0,ymax-ymin,10)
#plt.contour(grid_z0.T,extent=(kx.min(),kx.max(),ky.min(),ky.max()),levels = [-1.74],
#                 colors=('k',),linestyles=('-',),linewidths=(2,))

plt.savefig('ldos.png',dpi=100,bbox_inches='tight')
if (var.xplot):
	plt.show()
else:
	plt.savefig(var.output,bbox_inches='tight')



#	cb.set_label('Log(LDOS(E)/DOS(E))')
#	cb.set_label('LDOS(E)')

	#plt.gcf().set_size_inches(1000, 1000)
	#plt.show()
#	plt.tight_layout()
#	plt.savefig(str(foutname)+'En'+str(En[idx])+'.pdf',bbox_inches='tight')

