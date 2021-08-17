"""
Author: Andrea Gebek
Created on 16.7.2021
Plot a light curve
from a txt file.
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

wavelength = np.arange(grids_dict['lower_w'], grids_dict['upper_w'], grids_dict['resolution']) * 1e8 # In Angstrom

LightcurveFile = np.loadtxt('../../' + paramsFilename + '_lightcurve.txt')

orbphase = LightcurveFile[:, 0]

lightcurve_list = []

for idx in range(len(wavelength)):
    lightcurve_list.append(LightcurveFile[:, idx + 1])

"""
Plot the light curve and store the figure
"""

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

cmap = matplotlib.cm.get_cmap('Accent')

for idx, lightcurve in enumerate(lightcurve_list):
    ax.plot(orbphase, lightcurve, color = cmap(float(idx) / (len(wavelength) - 1.)), linewidth = 2, label = str(np.round(wavelength[idx], 2)) + r'$\,\AA$')

lg = ax.legend(loc = 'center right')
lg.get_frame().set_alpha(0)
lg.get_frame().set_linewidth(0)

ax.set_xlim(np.min(orbphase), np.max(orbphase))
ax.set_ylim(np.min(lightcurve_list) - 0.05 * (1 - np.min(lightcurve_list)), 1 + 0.05 * (1 - np.min(lightcurve_list)))

ax.set_xlabel(r'$\rm{Orbital \ Phase}$')
ax.set_ylabel(r'$\Re$')

ax.minorticks_on()
ax.tick_params(which = 'both', direction = 'in', right = True, top = True)

plt.tight_layout()

plt.savefig('../../' + paramsFilename + '_lightcurvePlot.pdf', dpi = 150)