import numpy as np
import matplotlib.pylab as plt
from scipy.signal import resample
from find_density import find_density
from find_velocity import find_vs
from find_velocity import find_vp
from find_g import find_g

#Read lookup tables
vp_pyrolite = np.loadtxt("Pyrolite_Stx08/cpl.pyrolite_anelastic_vp.dat", skiprows=4)
vs_pyrolite = np.loadtxt("Pyrolite_Stx08/cpl.pyrolite_anelastic_vs.dat", skiprows=4)
density_pyrolite = np.loadtxt("Pyrolite_Stx08/cpl.pyrolite_density.dat", skiprows=4)
g_array = np.loadtxt("PREM_g.dat")
nP = 1347        #points in P space
nT = 420         #points in T space
P0 = 1.37        #starting P (bar)
T0 = 300.004     #starting T (K)
dP = 1002.710    #step size in P (bar)
dT = 10.0238     #step size in T (K)
#g = 10.0

# Put data from lookup tables into an evenly spaced array
vp_array = np.zeros((nP,nT))
vs_array = np.zeros((nP,nT))
density_array = np.zeros((nP,nT))
for i in range(0,nP):
    for j in range(0,nT):
        global_index = i*nT + j
        vp_array[i,j] = vp_pyrolite[global_index]
        vs_array[i,j] = vs_pyrolite[global_index]
        density_array[i,j] = density_pyrolite[global_index]

#Read in temperature snapshot after postprocessing
T_Plume_data = np.loadtxt("Tpoints0055.dat")
T_Plume = T_Plume_data[:,2]
  
NPTH = 201                              #number of points in theta
NPR = 288                               #number of points in R
DPTH = 8.72664677444e-4                 #theta spacing
DPR = 3.48432055749e-3
DPR_dim = 10052.26478675
T_Plume_array = np.zeros((NPR,NPTH))      #Array containing T snapshot
x_axis = np.linspace(0.0,0.40,NPTH)       #From GMT.log
x_axis = x_axis*2885.0
r_axis = np.linspace(2.21369,1.21369,NPR)
r_axis_dim = 2885.0*np.linspace(2.21369,1.21369,NPR)
z_axis = (2.21369-r_axis) * (2885.0)    #dimensional depth axis

#Tpoints0052.dat contains points T values at points that are evenly spaced in r and theta.

#Fill T_Plume_array from T snapshot
for i in range(0,NPR):
    for j in range(0,NPTH):
        global_index = i*NPTH + j
        T_Plume_array[(NPR-1)-i,j] = T_Plume[global_index]

#Build the Pressure structure of the plume by finding the
#density of each grid point, and adding the weight of each
#vertical column
Rho_Plume_array = np.zeros((NPR,NPTH))
P_Plume_array = np.zeros((NPR,NPTH))
P_Plume_array[0,:] = 1.37   #Pressure is atmospheric for the first layer
Rho_Plume_array[0,:] = 3326.0

#From the lookup tables density depends on P and T.  Start
#by finding the density of the first layer (P = T = 0)
for i in range(1,NPR):
    for j in range(0,NPTH):
        if (r_axis_dim[i] > 6371.0):
           g = 9.8156
        else:
           g = find_g(r_axis_dim[i], g_array)
        #print "interpolated gravity is", g
        #g = 10.0
        if (i < (NPR-1)):
           Rho_Plume_array[i,j] = find_density(T_Plume_array[i,j], P_Plume_array[i,j], density_array)
           P_Plume_array[i+1,j] = (sum(Rho_Plume_array[:,j])*g*DPR_dim)/1e5   #P in bar
           if (P_Plume_array[i+1,j] > 1350000.0):  #for deep mantle, make sure P doesn't exceed lookup table range
              P_Plume_array[i+1,j] = 1350000.0
           #print "gravity is ", g, "pressure is ",P_Plume_array[i+1,j]
        if (i == (NPR-1)):
           Rho_Plume_array[i,j] = find_density(T_Plume_array[i,j], P_Plume_array[i,j], density_array)

#Use P and T conditions in the plume to find seismic velocity
Vs_Plume_array = np.zeros((NPR,NPTH))
Vp_Plume_array = np.zeros((NPR,NPTH))

for i in range(0,len(z_axis)):
    for j in range(0,len(x_axis)):
        #if(T_Plume_array[i,j] >  0):
        Vs_Plume_array[i,j] = find_vs(T_Plume_array[i,j], P_Plume_array[i,j], vs_array)
        Vp_Plume_array[i,j] = find_vp(T_Plume_array[i,j], P_Plume_array[i,j], vp_array)
        print Vs_Plume_array[i,j]

#Check reference density profile
#for i in range(0,NPR):
    #print z_axis[i], Rho_Plume_array[ i, (NPTH-1) ]


#Plots
plt.figure(1)
plt.pcolor(x_axis, r_axis_dim, T_Plume_array)
plt.xlabel('x (km)')
plt.ylabel('radius (km)')
plt.colorbar()

plt.figure(2)
plt.pcolor(x_axis, r_axis_dim, Vs_Plume_array)
plt.xlabel('x (km)')
plt.ylabel('radius (km)')
plt.colorbar()
plt.show()


