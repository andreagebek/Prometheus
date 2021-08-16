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

if not sphericalSymmetry:
    phi_steps = int(param['Grids']['phi_steps'])

z_steps = int(param['Grids']['z_steps'])

tauFile = np.loadtxt('../../' + paramsFilename + '_tau.txt')

if sphericalSymmetry:
    z = tauFile[:, 0]
    tau = tauFile[:, 1]

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)

    with np.errstate(divide = 'ignore'):
        ax.plot(z / R_0, np.log10(tau), linewidth = 4, color = 'blue')

    ax.set_xlabel(r'$z\,[R_0]$')
    ax.set_ylabel(r'$\log_{10}(\tau)$')

    ax.set_xlim(1, R_s / R_0)
    ax.set_ylim(-6, np.log10(np.max(tau)))

else:

    phiphi = tauFile[:, 0].reshape(phi_steps, z_steps)
    rhorho = tauFile[:, 1].reshape(phi_steps, z_steps)
    tau = tauFile[:, 2].reshape(phi_steps, z_steps)

    yy = rhorho * np.sin(phiphi) / R_0
    zz = rhorho * np.cos(phiphi) / R_0

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)

    with np.errstate(divide = 'ignore'):
        plt.scatter(yy, zz, c = np.log10(tau), s = 2, vmin = -3, vmax = 3, cmap = 'Spectral_r')

    cbar = plt.colorbar()
    cbar.set_label(r'$\log_{10}(\tau)$')


    ax.set_xlabel(r'$y\,[R_0]$')
    ax.set_ylabel(r'$z\,[R_0]$')

    ax.set_aspect('equal', adjustable='box')
ax.minorticks_on()
ax.tick_params(which = 'both', direction = 'in', right = True, top = True)

plt.tight_layout()

plt.savefig('../../' + paramsFilename + '_tauPlot.pdf', dpi = 50)