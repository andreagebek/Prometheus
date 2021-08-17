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

with open('../../' + paramsFilename + '.txt') as file:
    param = json.load(file)

grids_dict = param['Grids']

wavelength = (np.arange(grids_dict['lower_w'], grids_dict['upper_w'] + grids_dict['resolution'], grids_dict['resolution']) - grids_dict['resolution']  / 2.) * 1e8 # In Angstrom

LightcurveFile = np.loadtxt('../../' + paramsFilename + '_lightcurve.txt')

orbphase_border = grids_dict['orbphase_border']
orbphase_steps = grids_dict['orbphase_steps']
orbphase = np.linspace(-orbphase_border, orbphase_border + 2 * orbphase_border / (float(orbphase_steps) - 1.), 1 + int(orbphase_steps)) - orbphase_border / (float(orbphase_steps) - 1.)

wavelength_grid, orbphase_grid = np.meshgrid(wavelength, orbphase)

lightcurve_list = []

for idx in range(len(wavelength) - 1):
    lightcurve_list.append(LightcurveFile[:, idx + 1])

lightcurves = np.array(lightcurve_list)

"""
Plot the light curve and store the figure
"""

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

cmap = matplotlib.cm.get_cmap('Accent')

im = ax.pcolormesh(wavelength_grid, orbphase_grid / 2 * np.pi, lightcurves.T, cmap = 'Spectral')

cbar = fig.colorbar(im, ax = ax)
cbar.set_label(r'$\Re$')


ax.set_xlabel(r'$\lambda\,[\AA]$')
ax.set_ylabel(r'$\rm{Orbital \ Phase}$')


ax.minorticks_on()
ax.tick_params(which = 'both', direction = 'in', right = True, top = True)

plt.tight_layout()

plt.savefig('../../' + paramsFilename + '_spectralFieldPlot.pdf', dpi = 150)