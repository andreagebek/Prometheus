# coding=utf-8
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
Plotting settings - change ad lib
bandwidth: Bandwidth in Angstrom used to calculate the average lightcurve within that band
centerWavelengths: List of wavelengths (in Angstrom) which denote the centers of the bands 
"""

bandwidth = 0.4 # In Angstrom
centerWavelengths = [5880., 5890., 5891., 5893.] # In Angstrom


"""
Read in settings file and stored light curve
"""

paramsFilename = sys.argv[1]

with open(PARENTPATH + '/setupFiles/' + paramsFilename + '.txt') as file:
    param = json.load(file)

N_orbphase = int(param['Grids']['orbphase_steps']) # Get the number of orbital phase steps from the setup file


LightcurveFile = np.loadtxt(PARENTPATH + '/output/' + paramsFilename + '_lightcurve.txt')

wavelength = LightcurveFile[::N_orbphase, 0] # In Angstrom
N_wavelength = len(wavelength)
orbphase = LightcurveFile[0:N_orbphase, 1] / (2. * np.pi) # Convert to unity

R = LightcurveFile[:, 2].reshape(N_wavelength, N_orbphase) # Transit depth R(orbphase, wavelength)


lightcurveList = []

for w in centerWavelengths:

    SEL = (wavelength > w - bandwidth / 2.) * (wavelength < w + bandwidth / 2.)

    if len(SEL[SEL]) == 0:
        print(w)

    lightcurveList.append(np.mean(R[SEL, :], axis = 0))


"""
Plot the light curves and store the figure
"""

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

cmap = matplotlib.cm.get_cmap('Accent')

for idx, lightcurve in enumerate(lightcurveList):

    ax.plot(orbphase, lightcurve, color = cmap(float(idx) / (len(centerWavelengths) - 1.)), linewidth = 2, label = str(np.round(centerWavelengths[idx], 2)) + r'$\,\AA$')

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