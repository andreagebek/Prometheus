# coding=utf-8
"""
Author: Andrea Gebek
Created on 17.8.2021
Plot a field of transit depth
as a function of wavelength and time (orbital phase).
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
import prometheusScripts.fluxDecrease as flux

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
Read in settings file and stored light curve
"""

paramsFilename = sys.argv[1]

with open(PARENTPATH + '/setupFiles/' + paramsFilename + '.txt') as file:
    param = json.load(file)

architectureDict = param['Architecture']
gridsDict = param['Grids']

wavelength_axis = flux.constructAxis(gridsDict, architectureDict, 'wavelength')
orbphase_axis = flux.constructAxis(gridsDict, architectureDict, 'orbphase')

wavelength, orbphase = np.meshgrid(wavelength_axis, orbphase_axis, indexing = 'ij')

LightcurveFile = np.loadtxt(PARENTPATH + '/output/' + paramsFilename + '_lightcurve.txt')


R = LightcurveFile[:, 2].reshape(len(wavelength_axis), len(orbphase_axis))

"""
Plot the light curve and store the figure
"""

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

cmap = matplotlib.cm.get_cmap('Accent')
im = ax.pcolormesh(wavelength * 1e8, orbphase / (2. * np.pi), R, cmap = 'Spectral')

cbar = fig.colorbar(im, ax = ax)
cbar.set_label(r'$\Re$')
cbar.ax.minorticks_on()

ax.set_xlabel(r'$\lambda\,[\AA]$')
ax.set_ylabel(r'$\rm{Orbital \ Phase}$')


ax.minorticks_on()
ax.tick_params(which = 'both', direction = 'in', right = True, top = True)

plt.tight_layout()

plt.savefig(PARENTPATH + '/figures/' + paramsFilename + '_spectralFieldPlot.pdf', dpi = 150)