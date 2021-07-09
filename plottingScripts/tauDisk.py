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


with open('../../settings.txt') as file:
    param = json.load(file)

sphericalSymmetry = param['sphericalSymmetry']

R_0 = param['Architecture']['R_0']

if not sphericalSymmetry:
    phi_steps = int(param['Grids']['phi_steps'])

z_steps = int(param['Grids']['z_steps'])

outputFilename = param['Output']['outputFilename']

tauFile = np.loadtxt('../../' + outputFilename + '_tau.txt')

if sphericalSymmetry:
    #SOMETHINGSOMETHING
    0

else:

    phiphi = tauFile[:, 0].reshape(phi_steps, z_steps)
    rhorho = tauFile[:, 1].reshape(phi_steps, z_steps)
    tau = tauFile[:, 2].reshape(phi_steps, z_steps)
    
    yy = rhorho * np.sin(phiphi) / R_0
    zz = rhorho * np.cos(phiphi) / R_0

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)

    with np.errstate(divide = 'ignore'):
        plt.scatter(yy, zz, c = np.log10(tau), s = 2, vmin = -3, vmax = 3, cmap = 'Spectral')

cbar = plt.colorbar()
cbar.set_label(r'$\log_{10}(\tau)$')


ax.set_xlabel(r'$y\,[R_0]$')
ax.set_ylabel(r'$z\,[R_0]$')

ax.minorticks_on()
ax.tick_params(which = 'both', direction = 'in', right = True, top = True)

ax.set_aspect('equal', adjustable='box')

plt.tight_layout()

plt.savefig('../../' + outputFilename + '_tau.pdf', dpi = 150)