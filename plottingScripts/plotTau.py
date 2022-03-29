# coding=utf-8
"""
Author: Andrea Gebek
Created on 9.7.2021
Plot the optical depth at a certain
wavelength over the stellar disk.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import json
import sys
import os
SCRIPTPATH = os.path.realpath(__file__)
GITPATH = os.path.dirname(os.path.dirname(SCRIPTPATH))
PARENTPATH = os.path.dirname(GITPATH)
sys.path.append(GITPATH) 
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
Read in settings file and stored optical depth values
"""

plotPlanet = False
plotStar = True

paramsFilename = sys.argv[1]

with open(PARENTPATH + '/setupFiles/' + paramsFilename + '.txt') as file:
    param = json.load(file)

architectureDict = param['Architecture']
gridsDict = param['Grids']

phi_steps = int(gridsDict['phi_steps'])
rho_steps = int(gridsDict['rho_steps'])

R_0 = architectureDict['R_0']
R_star = architectureDict['R_star']

tauFile = np.loadtxt(PARENTPATH + '/output/' + paramsFilename + '_tau.txt')


phi = tauFile[:, 0].reshape(phi_steps, rho_steps)
rho = tauFile[:, 1].reshape(phi_steps, rho_steps)
tau = tauFile[:, 2].reshape(phi_steps, rho_steps)

x, y, z = geom.getCartesianFromCylinder(0, phi, rho) # x is not used here


fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)



with np.errstate(divide = 'ignore'):
    plt.scatter(y / R_0, z / R_0, c = np.log10(tau), s = 2, vmin = -3, vmax = 3, cmap = 'Spectral_r')

if plotPlanet:
    planetCircle = plt.Circle((0, 0), 1, color = 'black', linewidth = 0)
    ax.add_patch(planetCircle)

if plotStar:
    starCircle = plt.Circle((0, 0), R_star / R_0, color = 'black', fill = False, linewidth = 1)
    ax.add_patch(starCircle)

cbar = plt.colorbar()
cbar.set_label(r'$\log_{10}(\tau)$')
cbar.ax.minorticks_on()

ax.set_xlabel(r'$y\,[R_0]$')
ax.set_ylabel(r'$z\,[R_0]$')

ax.set_xlim(-1.05 * R_star / R_0, 1.05 * R_star / R_0)
ax.set_ylim(-1.05 * R_star / R_0, 1.05 * R_star / R_0)

ax.set_aspect('equal', adjustable='box')

ax.minorticks_on()
ax.tick_params(which = 'both', direction = 'in', right = True, top = True)

plt.tight_layout()

plt.savefig(PARENTPATH + '/figures/' + paramsFilename + '_tauPlot.pdf', dpi = 50)