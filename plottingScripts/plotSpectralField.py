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
Read in settings file and stored lightcurve.
"""

paramsFilename = sys.argv[1]

with open(PARENTPATH + '/setupFiles/' + paramsFilename + '.txt') as file:
    param = json.load(file)

architectureDict = param['Architecture']
gridsDict = param['Grids']



R_star = architectureDict['R_star']
R_0 = architectureDict['R_0']
a_p = architectureDict['a_p']


LightcurveFile = np.loadtxt(PARENTPATH + '/output/' + paramsFilename + '_lightcurve.txt')

N_orbphase = int(param['Grids']['orbphase_steps']) # Get the number of orbital phase steps from the setup file
N_wavelength = int(np.size(LightcurveFile, axis = 0) / float(N_orbphase))

wavelength = LightcurveFile[::N_orbphase, 0] # In Angstrom
orbphase = LightcurveFile[0:N_orbphase, 1] / (2. * np.pi) # Convert to unity
R = LightcurveFile[:, 2].reshape(N_wavelength, N_orbphase) # Transit depth R(orbphase, wavelength)

"""
Plot the light curve and store the figure
"""

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

cmap = matplotlib.cm.get_cmap('Accent')
im = ax.pcolormesh(wavelength, orbphase, R.T, cmap = 'Spectral', shading = 'nearest')

cbar = fig.colorbar(im, ax = ax)
cbar.set_label(r'$\Re$')
cbar.ax.minorticks_on()

ax.set_xlabel(r'$\lambda\,[\AA]$')
ax.set_ylabel(r'$\rm{Orbital \ Phase}$')


ax.minorticks_on()
ax.tick_params(which = 'both', direction = 'in', right = True, top = True)

plt.tight_layout()

plt.savefig(PARENTPATH + '/figures/' + paramsFilename + '_spectralFieldPlot.pdf', dpi = 150)