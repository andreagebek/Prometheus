"""
Author: Andrea Gebek
Created on 16.8.2021
Plot the line-of-sight velocity due to
stellar rotation over the stellar disk.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import json
import sys

matplotlib.rcParams['axes.linewidth'] = 2.5
matplotlib.rcParams['xtick.major.size'] = 10
matplotlib.rcParams['xtick.minor.size'] = 6
matplotlib.rcParams['xtick.major.width'] = 2.5
matplotlib.rcParams['xtick.minor.width'] = 1.5
matplotlib.rcParams['ytick.major.size'] = 10
matplotlib.rcParams['ytick.minor.size'] = 6
matplotlib.rcParams['ytick.major.width'] = 2.5
matplotlib.rcParams['ytick.minor.width'] = 1.5
matplotlib.rcParams['image.origin'] = 'lower'
matplotlib.rcParams.update({'font.size': 26, 'font.weight': 'bold'})

"""
Read in settings file and stored optical depth values
"""

paramsFilename = sys.argv[1]

with open('../../' + paramsFilename + '.txt') as file:
    param = json.load(file)

sphericalSymmetry = param['sphericalSymmetry']

R_0 = param['Architecture']['R_0']
R_s = param['Architecture']['R_star']
period_starrot = param['Architecture']['period_starrot']
inclination_starrot = param['Architecture']['inclination_starrot']
azimuth_starrot = param['Architecture']['azimuth_starrot']

phi_steps = int(param['Grids']['phi_steps'])
z_steps = int(param['Grids']['z_steps'])


phi = np.linspace(0, 2 * np.pi, int(phi_steps) + 1, dtype = np.dtype('f4'))[:-1] + np.pi / float(phi_steps)
rho = np.linspace(0, R_s, int(z_steps) + 1, dtype = np.dtype('f4'))[:-1] + 0.5 * R_s / float(z_steps)

phiphi, rhorho = np.meshgrid(phi, rho)

dir_omega = np.array([-np.sin(inclination_starrot) * np.cos(azimuth_starrot), -np.sin(inclination_starrot) * np.sin(azimuth_starrot), np.cos(inclination_starrot)])
omega = 2. * np.pi / period_starrot * dir_omega # Angular velocity vector of the stellar rotation
r_surface = np.array([np.tensordot(np.sqrt(R_s**2 - rho**2), np.ones(len(phi)), axes = 0), np.tensordot(rho, np.sin(phi), axes = 0), np.tensordot(rho, np.cos(phi), axes = 0)]) 
# Vector to the surface of the star
v_los = np.cross(omega, r_surface, axisb = 0)[:, :, 0] * 1e-5 # The line-of-sight velocity is the one along the x-axis (in km/s)


yy = rhorho * np.sin(phiphi) / R_0
zz = rhorho * np.cos(phiphi) / R_0

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

plt.scatter(yy, zz, c = v_los, s = 2, vmin = np.min(v_los), vmax = np.max(v_los), cmap = 'Spectral_r')

cbar = plt.colorbar()
cbar.set_label(r'$v_{\rm{los}}\,[\rm{km/s}]$')


ax.set_xlabel(r'$y\,[R_0]$')
ax.set_ylabel(r'$z\,[R_0]$')

ax.set_aspect('equal', adjustable='box')
ax.minorticks_on()
ax.tick_params(which = 'both', direction = 'in', right = True, top = True)

plt.tight_layout()

plt.savefig('../../' + paramsFilename + '_stellarRotationPlot.pdf', dpi = 50)