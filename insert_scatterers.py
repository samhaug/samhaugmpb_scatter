#!/usr/bin/env python

import numpy as np
import Scatterer_Position as sp
import hetero_structure as hs

#Make velocity structure

vel_structure = hs.Mantle_Structure('try','Synthetic_vel_structure.csv')

interp_radius = np.linspace(3500.4,6356.6,num=4000)
interp_theta = np.linspace(0,180,num=1800)

vel_structure.array_interp2D(interp_radius,interp_theta)

MORB = np.loadtxt('MORB_output.dat')

MORB_pol = np.zeros(MORB.shape)

MORB_pol[:,0] = np.sqrt(MORB[:,0]**2+MORB[:,1]**2)
MORB_pol[:,1] = np.degrees(np.arctan2(MORB[:,0],MORB[:,1]))

vp_reference = vel_structure.vp_2D
vs_reference = vel_structure.vs_2D
rho_reference = vel_structure.rho_2D

for ii in range(0,len(MORB_pol)):
    r = np.argmin(np.abs(vel_structure.radius-MORB_pol[ii,0]))
    th = np.argmin(np.abs(vel_structure.theta-MORB_pol[ii,1]))
    new_vp = vp_reference[r,th]*(1-0.04)
    new_vs = vs_reference[r,th]*(1-0.04)
    new_rho = rho_reference[r,th]*(1-0.04)

    vel_structure.add_hetero_point(MORB_pol[ii,0],MORB_pol[ii,1],new_vp,
            new_vs,new_rho,10.)
