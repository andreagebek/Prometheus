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

wavelength = flux.constructAxis(gridsDict, architectureDict, 'wavelength')
orbphase = flux.constructAxis(gridsDict, architectureDict, 'orbphase')


LightcurveFile = np.loadtxt(PARENTPATH + '/output/' + paramsFilename + '_lightcurve.txt')

R = LightcurveFile[:, 2].reshape(len(wavelength), len(orbphase))

lightcurveList = []

for idx in range(len(wavelength)):
    lightcurveList.append(R[idx, :])


"""
Plot the light curve and store the figure
"""

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

cmap = matplotlib.cm.get_cmap('Accent')

for idx, lightcurve in enumerate(lightcurveList):

    ax.plot(orbphase, lightcurve, color = cmap(float(idx) / (len(wavelength) - 1.)), linewidth = 2, label = str(np.round(wavelength[idx] * 1e8, 2)) + r'$\,\AA$')

lg = ax.legend(loc = 'center right')
lg.get_frame().set_alpha(0)
lg.get_frame().set_linewidth(0)

ax.set_xlim(np.min(orbphase), np.max(orbphase))
ax.set_ylim(np.min(R) - 0.05 * (1 - np.min(R)), 1 + 0.05 * (1 - np.min(R)))

ax.set_xlabel(r'$\rm{Orbital \ Phase}$')
ax.set_ylabel(r'$\Re$')

ax.minorticks_on()
ax.tick_params(which = 'both', direction = 'in', right = True, top = True)

plt.tight_layout()

plt.savefig(PARENTPATH + '/figures/' + paramsFilename + '_lightcurvePlot.pdf', dpi = 150)