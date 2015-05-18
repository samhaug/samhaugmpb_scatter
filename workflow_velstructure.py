import numpy as np
import seis_vel as sv
import os
import subprocess as subp
import matplotlib.pyplot as plt
import matplotlib
import Scatterer_Position as sp

#matplotlib.use('PDF')

TandC_array = np.loadtxt('TEMP_output.dat')
vp_file = ('cpl.pyrolite_anelastic_vp.dat')
vs_file = ('cpl.pyrolite_anelastic_vs.dat')
rho_file = ('cpl.pyrolite_density.dat')

T0, P0, dT, dP, nT, nP = sv.get_base_values('cpl.pyrolite_anelastic_vp.dat')

vp_lookup, vs_lookup, density_lookup = sv.array_generator(vp_file, vs_file, rho_file)

# Make XY into polar coordinates with theta clockwise from northpole
Theta_array = np.abs(np.degrees(np.arctan2(TandC_array[:,0],TandC_array[:,1])))

# Determine Radius
Radius_array = np.sqrt(TandC_array[:,0]**2+TandC_array[:,1]**2)

Coord_array = np.vstack((Radius_array,Theta_array)).transpose()
Polar_Temp_array = np.round(np.hstack((Coord_array,TandC_array[:,2][:,None])),decimals=1)

#Place unique radius and theta values into array.
Radius_array = np.unique(Polar_Temp_array[:,0])
Theta_array = np.unique(Polar_Temp_array[:,1])

Temp_array = np.zeros((len(Radius_array),len(Theta_array)))

for ii in range(0,len(Polar_Temp_array)):
    rad = np.argmin(np.abs(Radius_array-Polar_Temp_array[ii,0]))
    th = np.argmin(np.abs(Theta_array-Polar_Temp_array[ii,1]))
    Temp_array[rad,th] = Polar_Temp_array[ii,2]

# Iteratively find pressure conditions throughout array. Use PT conditions to find
# Density which is used for lithostatic pressure calculations. P=0 at surface.

#Allocate arrays for pressure, velocity and density. Make same shape as Temp_array
Density_array = np.zeros(Temp_array.shape)
Pressure_array = np.zeros(Temp_array.shape)
Vp_array = np.zeros(Temp_array.shape)
Vs_array = np.zeros(Temp_array.shape)

#Flip upside down to make for loop indexing more intuitive
Temp_array = np.flipud(Temp_array)
Radius_array = np.flipud(Radius_array)

for ii in range(0,Temp_array.shape[1]):
    Pressure_array[0,ii] = 0
    Density_array[0,ii] = sv.find_density(Temp_array[0,ii], Pressure_array[0,ii],
            density_lookup, P0, T0, dP, dT)
    Vp_array[0,ii] = sv.find_vp(Temp_array[0,ii], Pressure_array[0,ii],
            vp_lookup, P0, T0, dP, dT)
    Vs_array[0,ii] = sv.find_vs(Temp_array[0,ii], Pressure_array[0,ii],
            vs_lookup, P0, T0, dP, dT)

#Now compute values iteratively starting on the second row of the array.
for ii in range(1,Temp_array.shape[0]):
    for jj in range(0,Temp_array.shape[1]):
        Pressure_array[ii,jj] = Pressure_array[ii-1,jj]+Density_array[ii-1,jj]*\
                np.round((np.abs(Radius_array[ii]-Radius_array[ii-1])),decimals=1)*1000.*\
                9.81/pow(10.,9.)
        Density_array[ii,jj] = sv.find_density(Temp_array[ii-1,jj], Pressure_array[ii-1,jj],
                density_lookup, P0, T0, dP, dT)
        Vp_array[ii,jj] = sv.find_vp(Temp_array[ii,jj], Pressure_array[ii,jj],
                vp_lookup, P0, T0, dP, dT)
        Vs_array[ii,jj] = sv.find_vs(Temp_array[ii,jj], Pressure_array[ii,jj],
                vs_lookup, P0, T0, dP, dT)

theta, radius = np.meshgrid(Theta_array,Radius_array)

Pressure_column = Pressure_array.reshape((Pressure_array.size,1))
Density_column = Density_array.reshape((Density_array.size,1))
Vp_column = Vp_array.reshape((Vp_array.size,1))
Vs_column = Vs_array.reshape((Vs_array.size,1))
theta_column = theta.reshape((theta.size,1))
radius_column = radius.reshape((radius.size,1))

file = open('Synthetic_vel_structure.csv','w+')
#file.write('Radius (km)   theta (deg from NP)      Vp (km/s)     Vs (km/s)     rho (kg/m^3) \n')
for ii in range(0,len(Pressure_column)):
    file.write(str(radius_column[ii][0])+' '+str(theta_column[ii][0])+' '+str(Vp_column[ii][0])+' '+
        str(Vs_column[ii][0])+' '+str(Density_column[ii][0])+'\n')

file.close()
