import numpy as np
import seis_vel as sv
import os
import subprocess as sp
import matplotlib.pyplot as plt
import matplotlib

#matplotlib.use('PDF')
TandC_array = np.loadtxt('TEMP_output.dat')
vp_file = ('cpl.pyrolite_anelastic_vp.dat')
vs_file = ('cpl.pyrolite_anelastic_vs.dat')
rho_file = ('cpl.pyrolite_density.dat')

T0, P0, dT, dP, nT, nP = sv.get_base_values('cpl.pyrolite_anelastic_vp.dat')

vp_array, vs_array, density_array = sv.array_generator(vp_file, vs_file, rho_file)

height = TandC_array.shape[0]
width = TandC_array.shape[1]
x_coord = TandC_array[:,0]
y_coord = TandC_array[:,1]
temperature = TandC_array[:,2]

pressure = list()
radius = list()
theta = list()

for ii in range(0,len(x_coord)):
      p, r, t = sv.xy_2_pressure(x_coord[ii],y_coord[ii])
      pressure.append(p)
      radius.append(r)
      theta.append(t)

pressure = np.array(pressure)
radius = np.array(radius)
theta = np.array(theta)

#stack the pressure, radius, theta columns and sort by largest radius to smallest
reshape = np.column_stack((radius,theta,temperature))
#reshape = np.flipud(reshape[reshape[:,0].argsort()])
reshape = np.flipud(reshape)
reshape[:,1]=np.round(reshape[:,1]*180/np.pi)
reshape[:,0]=np.round(reshape[:,0]/10.)*10
#Determine how many "incriments" the radial unit is stepped in. divide by 10
#to only account for large stepwise discontinutities. see plt.plot(reshape[:,0]). k
#is the number of radial layers.

k=0
for ii in range(0,reshape.shape[0]-1):
   if reshape[ii,1] == 0.:
     k=k+1
radial_step = (reshape[:,0].shape[0])/k

#Convert reshape into Temp_array: an x-y array on which to determine lithostatic pressure.

Temp_array = np.zeros([k+1,182])
for ii in range(0,k+1):
   for jj in range(0,182):
      Temp_array[ii,jj] = reshape[int(180*ii+jj)+ii-1,2]

Temp_array = Temp_array

#Pressure at top layer is 0 Pa. Use temperature and pressure to find density of top layer:

Rho_grid = np.zeros([k+1,181])
vs_grid = np.zeros([k+1,181])
vp_grid = np.zeros([k+1,181])

theta_line = np.linspace(180,0,181)
theta_grid = np.array([theta_line]*101)

radius_line = np.linspace(6371,3481,101)
radius_grid = np.array([radius_line]*181).transpose()

#meters per radial step:

mpr = 2890000./k
P_layer = np.zeros(181)*P0
for jj in range(0,181):
   Rho_grid[0,jj] = sv.find_density(Temp_array[0,jj],0,density_array,P0,T0,dP,dT)
   vs_grid[0,jj] = sv.find_vs(Temp_array[0,jj],0,vs_array,P0,T0,dP,dT)
   vp_grid[0,jj] = sv.find_vp(Temp_array[0,jj],0,vp_array,P0,T0,dP,dT)

P_layer = P_layer+(Rho_grid[0,:]*9.81*mpr*1.)/(1e9)

#Now for each sucsessive layer, find new density from P_layer pressure and Temp_array
#temperature. Refine P_layer with new lithostatic contribution from the last layer.

for ii in range(1,k+1):
   for jj in range(0,181):
      Rho_grid[ii,jj] = sv.find_density(Temp_array[ii,jj],P_layer[jj],density_array,P0,T0,dP,dT)
      vs_grid[ii,jj] = sv.find_vs(Temp_array[ii,jj],P_layer[jj],vs_array,P0,T0,dP,dT)
      vp_grid[ii,jj] = sv.find_vp(Temp_array[ii,jj],P_layer[jj],vp_array,P0,T0,dP,dT)

   
   P_layer = P_layer+(Rho_grid[ii,:]*9.81*mpr)/(1e9)
   
#Now to organize into the output file
#Convert radius_grid, theta_grid, Rho_grid, vs_grid, vp_grid to columns

Rho_column = np.reshape(Rho_grid,Rho_grid.size)
vs_column = np.reshape(vs_grid,vs_grid.size)
vp_column = np.reshape(vp_grid,vp_grid.size)
radius_column = np.reshape(radius_grid,radius_grid.size)
theta_column = np.reshape(theta_grid,theta_grid.size)

#Remove lithosphere (top 10 rows) and replace with the average of the 11th row
Temp_array_nolith = Temp_array
Temp_array_nolith[0:9,:] = Temp_array[10,:].mean()
Avg_Temp_nolith = Temp_array_nolith.mean(1)

vs_grid_nolith = vs_grid      
vs_grid_nolith[0:9,:] = vs_grid[10,:].mean()
vs_column_nolith = np.reshape(vs_grid_nolith,vs_grid_nolith.size)
avg_vs_nolith = vs_grid_nolith.mean(1)

vp_grid_nolith = vp_grid      
vp_grid_nolith[0:9,:] = vp_grid[10,:].mean()
vp_column_nolith = np.reshape(vp_grid_nolith,vp_grid_nolith.size)
avg_vp_nolith = vp_grid_nolith.mean(1)

Rho_grid_nolith = Rho_grid
Rho_grid_nolith[0:9,:] = Rho_grid[10,:].mean()
Rho_column_nolith = np.reshape(Rho_grid_nolith,Rho_grid_nolith.size)
avg_rho_nolith = Rho_grid_nolith.mean(1)


Rho_column_nolith = np.reshape(Rho_grid_nolith,Rho_grid_nolith.size)
vs_column_nolith = np.reshape(vs_grid_nolith,vs_grid_nolith.size)
vp_column_nolith = np.reshape(vp_grid_nolith,vp_grid_nolith.size)
#Make 1D array for reference. Use average values from 2D array.

Avg_Temp = Temp_array.mean(1)
avg_rho = Rho_grid.mean(1)
avg_vs = vs_grid.mean(1)
avg_vp = vp_grid.mean(1)

'''
#Make 1D reference file NOLITH
file = open('Sam_ref_nolith.bm','w+')

file.write('ANELASTIC     F \n')
file.write('ANISOTROPIC   F \n')
file.write('UNITS        m \n')
file.write('COLUMNS    radius    rho    vpv    vsv\n')

for ii in range(0,len(radius_grid)):
   radius = int(radius_grid[ii,0]*1000)
   rho = int(avg_rho_nolith[ii])
   vpv = int(avg_vp_nolith[ii]*1000)
   vsv = int(avg_vs_nolith[ii]*1000)
   file.write(str(radius)+' '+str(rho)+' '+str(vpv)+' '+str(vsv)+'\n')
file.close()

#Make 1D reference file WITH_LITH
file = open('Sam_ref.bm','w+')

file.write('ANELASTIC     F \n')
file.write('ANISOTROPIC   F \n')
file.write('UNITS        m \n')
file.write('COLUMNS    radius    rho    vpv    vsv\n')

for ii in range(0,len(radius_grid)):
   radius = int(radius_grid[ii,0]*1000)
   rho = int(avg_rho[ii])
   vpv = int(avg_vp[ii]*1000)
   vsv = int(avg_vs[ii]*1000)
   file.write(str(radius)+' '+str(rho)+' '+str(vpv)+' '+str(vsv)+'\n')
file.close()

#Make velocity structure file NOLITH
file_length=Rho_column.size

if os.path.exists('./vel_structure_nolith.sph'):
   sp.call('rm ./vel_structure_nolith.sph',shell=True)
   
file = open('vel_structure_nolith.sph','w+')
file.write(str(file_length)+'\n')

for ii in range(0,Rho_grid_nolith.size):
   file.write(str(radius_column[ii])+' '+str(theta_column[ii])+' '+ \
   str(vp_column_nolith[ii]*1000.)+' '+str(vs_column_nolith[ii]*1000.)+' '+str(Rho_column_nolith[ii]))
   file.write('\n')

file.close()
'''
#Make velocity structure file WITH_LITH
file_length=Rho_column.size

if os.path.exists('./vel_structure.sph'):
   sp.call('rm ./vel_structure.sph',shell=True)
   
file = open('vel_structure.sph','w+')
file.write(str(file_length)+'\n')

for ii in range(0,Rho_grid.size):
   file.write(str(radius_column[ii])+' '+str(theta_column[ii])+' '+ \
   str(vp_column[ii]*1000.)+' '+str(vs_column[ii]*1000.)+' '+str(Rho_column[ii]))
   file.write('\n')

file.close()

