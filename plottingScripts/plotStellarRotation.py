# coding=utf-8
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
import os
import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable
SCRIPTPATH = os.path.realpath(__file__)
GITPATH = os.path.dirname(os.path.dirname(SCRIPTPATH))
PARENTPATH = os.path.dirname(GITPATH)
sys.path.append(GITPATH) 
import prometheusScripts.fluxDecrease as flux
import prometheusScripts.geometryHandler as geom

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
Plotting settings - change ad lib
"""

plotStar = True

"""
Read in settings file and stored optical depth values
"""

paramsFilename = sys.argv[1]

with open(PARENTPATH + '/setupFiles/' + paramsFilename + '.txt') as file:
    param = json.load(file)

architectureDict = param['Architecture']
gridsDict = param['Grids']

R_0 = architectureDict['R_0']
R_star = architectureDict['R_star']

phi_axis = flux.constructAxis(gridsDict, 'phi')[0]
rho_axis = flux.constructAxis(gridsDict, 'rho')[0]

phiGrid, rhoGrid = np.meshgrid(phi_axis, rho_axis, indexing = 'ij')

v_los = []

for phi in phi_axis:
    for rho in rho_axis:

        v_los.append(geom.calculateStellarLOSvelocity(architectureDict, phi, rho) * 1e-5) # The line-of-sight velocity is the one along the x-axis (in km/s)

v_los = np.array(v_los).reshape(len(phi_axis), len(rho_axis))

y, z = geom.getCartesianFromCylinder(phiGrid, rhoGrid)


fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

im = plt.scatter(y / R_0, z / R_0, c = v_los, s = 2, vmin = np.min(v_los), vmax = np.max(v_los), cmap = 'Spectral_r')

if plotStar:
    starCircle = plt.Circle((0, 0), R_star / R_0, edgecolor = 'black', fill = False, linewidth = 2)
    ax.add_patch(starCircle)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.25)

cbar = plt.colorbar(im, cax=cax)
cbar.set_label(r'$v_{\rm{los}}\,[\rm{km/s}]$')
cbar.ax.minorticks_on()


ax.set_xlabel(r'$y\,[R_0]$')
ax.set_ylabel(r'$z\,[R_0]$')

ax.set_aspect('equal', adjustable='box')
ax.minorticks_on()
ax.tick_params(which = 'both', direction = 'in', right = True, top = True)

plt.tight_layout()

plt.savefig(PARENTPATH + '/figures/' + paramsFilename + '_stellarRotationPlot.pdf', dpi = 50)