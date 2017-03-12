#!/usr/bin/python2.7                           
                                               
import matplotlib.pyplot as plt			
import numpy as np				
import sys					
						
						                                               
import pyshtools as sht                        	
from pyshtools.expand import SHExpandDH;       	
from pyshtools.expand import MakeGridDH;       	
from scipy.special import sph_harm             	
					       	
from matplotlib import cm, colors		
from mpl_toolkits.mplot3d import Axes3D		
					       	
import dft_su_sfera as dfts		       	



def get_distro(filename, lmax):
	"reads file formatted as: # cos/sin ## l ## m ## coeff # ... coefficients for the starting distribution"
	INPUT = open(filename,'r')
	
	Clm = np.zeros((2,lmax+1,lmax+1), order='F')
	
	for line in INPUT:
		fields = line.split()
		u = int(fields[0])
		l = int(fields[1])
		m = int(fields[2])
		C = float(fields[3])
		if l < lmax+1:
			Clm[u,l,m] = C
	
	INPUT.close()
	
	return Clm


def show(F):
	"show function on 2D plane"
	
	
	N = F.shape[0]
	print F.shape

	theta = np.linspace(0.,np.pi,N)
	phi   = np.linspace(0.,2*np.pi,2*N)

	PHI,THETA = np.meshgrid(phi, theta)
	
	fig, ax = plt.subplots(1,1,figsize=(10,5))
	ax.imshow(F, extent=(0,2*np.pi,0,np.pi), cmap='magma')
	ax.set(xlabel='longitude', ylabel='latitude');
	
	plt.show()
	
def show_sphere(F):
	"show function on the sphere"
	
	N = F.shape[0]
	print N
	
	theta = np.linspace(0.,np.pi,N)
	phi   = np.linspace(0.,2*np.pi,2*N)
	
	phi, theta = np.meshgrid(phi,theta)
	
	x = np.sin(theta) * np.cos(phi)
	y = np.sin(theta) * np.sin(phi)
	z = np.cos(theta)
	
	fig = plt.figure(figsize=plt.figaspect(1.))
	ax = fig.add_subplot(111, projection='3d')
	
	ax.plot_surface(x, y, z,  rstride=1, cstride=1, facecolors=cm.seismic(F))
	
	ax.set_axis_off()
	plt.show()
	
	
###########################################################################

L_MAX = int(sys.argv[1])
DISTRO = sys.argv[2]

rho_ = get_distro(DISTRO, L_MAX)

rho = MakeGridDH(rho_, sampling=2, norm=4)

show(rho)

show_sphere(rho)


















					       	
