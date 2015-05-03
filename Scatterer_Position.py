#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

###############################################################################
def cart2pol(x,y):
###############################################################################
    rad = np.sqrt(x**2 + y**2)
    phi = np.degrees(np.arctan2(y,x))
    if phi < 0:
  	    phi = abs(phi)+180
    return rad, phi
###############################################################################
def pol2cart(rad, phi):
###############################################################################
    x = rad * np.cos(np.radians(phi))
    y = rad * np.sin(np.radians(phi))
    return x, y

###############################################################################
class Scatterer_Position(object):
###############################################################################

# Initialize with the name of the .dat file containing the tracer information.
# Creates arrays in cartesian and polar.
###############################################################################
    def __init__(self,MORB_array,HARZ_array,TEMP_array):
###############################################################################
        self.MORB_array = np.loadtxt(MORB_array)
        self.HARZ_array = np.loadtxt(HARZ_array)
        self.TEMP_array = np.loadtxt(TEMP_array)
        
        self.MORB_polar = np.empty((self.MORB_array.shape[0],
            self.MORB_array.shape[1]))
        self.HARZ_polar = np.empty((self.HARZ_array.shape[0],
            self.HARZ_array.shape[1]))
        self.TEMP_polar = np.empty((self.TEMP_array.shape[0],
            self.TEMP_array.shape[1]))

        for ii in range(0,len(self.MORB_array)):
            self.MORB_polar[ii,0], self.MORB_polar[ii,1] = \
            cart2pol(self.MORB_array[ii,0],self.MORB_array[ii,1])
        for ii in range(0,len(self.HARZ_array)):
            self.HARZ_polar[ii,0], self.HARZ_polar[ii,1] = \
            cart2pol(self.HARZ_array[ii,0],self.HARZ_array[ii,1])
        for ii in range(0,len(self.TEMP_array)):
            self.TEMP_polar[ii,0], self.TEMP_polar[ii,1] = \
            cart2pol(self.TEMP_array[ii,0],self.TEMP_array[ii,1])

        self.TEMP_polar[:,2] = self.TEMP_array[:,2]
        self.TEMP_polar[:,3] = self.TEMP_array[:,3]

# Determine spacing intervals for radius and theta.
        self.min_rad = np.amin(np.floor(self.TEMP_polar[:,0])) 
        self.imin_rad = np.where(np.floor(self.TEMP_polar[:,0]) == self.min_rad)[0].size
        self.rad_interval = self.TEMP_polar[:,0].size/self.imin_rad

        self.min_th = np.amin(np.floor(self.TEMP_polar[:,1]))
        self.imin_th = np.where(np.floor(self.TEMP_polar[:,1]) == self.min_th)[0].size
        self.th_interval = self.TEMP_polar[:,1].size/self.imin_th

        self.polar_temp_array = -1.*np.ones((self.rad_interval,self.th_interval+1))
        self.radius = list()
        self.theta = list()

        for ii in range(0,len(self.TEMP_polar)):
			radius = np.round((self.TEMP_polar[ii,0]-self.min_rad)/29.)
			theta = np.round(self.TEMP_polar[ii,1])
			
			self.radius.append(int(radius))
			self.theta.append(int(theta))

			self.polar_temp_array[radius,theta] = self.TEMP_polar[ii,2]


###############################################################################
    def preview(self):
###############################################################################
		fig = plt.figure(figsize=(12,12))
		ax1 = fig.add_subplot(221)
		c = plt.scatter(self.MORB_array[:,0],self.MORB_array[:,1],s=0.1,marker='.')
		ax1.set_title('MORB tracers')
		ax1.set_xlim(-6371,6371)
	   	ax1.set_ylim(-6371,6371)
	   	ax1.plot([0,0],[7000,0],'y-',lw=2)
	   	ax1.plot([0,0],[-7000,0],'y-',lw=2)
        

		ax2 = fig.add_subplot(222)
		c = plt.scatter(self.HARZ_array[:,0],self.HARZ_array[:,1],s=0.1,marker='.')
		ax2.set_title('HARZ tracers')
		ax2.set_xlim(-6371,6371)
		ax2.set_ylim(-6371,6371)
	   	ax2.plot([0,0],[7000,0],'y-',lw=2)
	   	ax2.plot([0,0],[-7000,0],'y-',lw=2)
		
		ax3 = fig.add_subplot(223)
		c = plt.scatter(self.MORB_array[::3,0],self.MORB_array[::3,1],s=0.1,marker='.')
		ax3.set_title('Every 3rd MORB tracer')
		ax3.set_xlim(-6371,6371)
	   	ax3.set_ylim(-6371,6371)
	   	ax3.plot([0,0],[7000,0],'y-',lw=2)
	   	ax3.plot([0,0],[-7000,0],'y-',lw=2)

		ax4 = fig.add_subplot(224)

		azimuths = np.radians(np.linspace(min(self.TEMP_polar[:,1]),max(self.TEMP_polar[:,1]),
				num=self.TEMP_polar.shape[0]))
		zeniths = np.linspace(min(self.TEMP_polar[:,0]),max(self.TEMP_polar[:,0]),
				num=self.TEMP_polar.shape[0])
		r, theta = np.meshgrid(zeniths,azimuths)

		c = plt.pcolor(self.TEMP_array[:,0],self.TEMP_array[:,1],self.TEMP_array[:,2])
		ax4.set_title('Temperature')
		plt.show()
###############################################################################
    def rotate(self,degrees):
###############################################################################

# Rotate scatterer positions by degrees counterclockwise.
        self.MORB_polar[:,1] = self.MORB_polar[:,1]+degrees 
        self.HARZ_polar[:,1] = self.HARZ_polar[:,1]+degrees 
        self.TEMP_polar[:,1] = self.TEMP_polar[:,1]+degrees 
        for ii in range(0,len(self.MORB_polar)):
            if self.MORB_polar[ii,1] > 360:
                self.MORB_polar[ii,1] = self.MORB_polar[ii,1]-360.
            if self.HARZ_polar[ii,1] > 360:
                self.HARZ_polar[ii,1] = self.HARZ_polar[ii,1]-360.
    	    
        for ii in range(0,len(self.MORB_polar)):
            self.MORB_array[ii,0], self.MORB_array[ii,1] = \
            pol2cart(self.MORB_polar[ii,0],self.MORB_polar[ii,1])
            self.HARZ_array[ii,0], self.HARZ_array[ii,1] = \
            pol2cart(self.HARZ_polar[ii,0],self.HARZ_polar[ii,1])
