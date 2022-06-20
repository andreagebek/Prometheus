# coding=utf-8
"""
Author: Andrea Gebek
Created on 3.9.2021
Plot a transmission spectrum
from a lightcurve txt file.
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
"""

Normalize = False # Do a quick normalization of the spectra

"""
Read in settings file and stored lightcurve.
Calculate the transmission spectrum by averaging
the depth in the lightcurve during full occultation
of the planet's R_0.
"""

paramsFilename = sys.argv[1]

with open(PARENTPATH + '/setupFiles/' + paramsFilename + '.txt') as file:
    param = json.load(file)

architectureDict = param['Architecture']
gridsDict = param['Grids']



R_star = architectureDict['R_star']
R_0 = architectureDict['R_0']
a_p = architectureDict['a_p']

orbphaseFullIngress = np.arcsin((R_star - R_0) / a_p) / (2. * np.pi) # Orbital phase at which the planet's R_0 is fully within the stellar disk


LightcurveFile = np.loadtxt(PARENTPATH + '/output/' + paramsFilename + '_lightcurve.txt')

N_orbphase = int(param['Grids']['orbphase_steps']) # Get the number of orbital phase steps from the setup file
N_wavelength = int(np.size(LightcurveFile, axis = 0) / float(N_orbphase))

wavelength = LightcurveFile[::N_orbphase, 0] # In Angstrom
orbphase = LightcurveFile[0:N_orbphase, 1] / (2. * np.pi) # Convert to unity
R = LightcurveFile[:, 2].reshape(N_wavelength, N_orbphase) # Transit depth R(orbphase, wavelength)

SEL = (orbphase > -orbphaseFullIngress) * (orbphase < orbphaseFullIngress)

if not any(SEL):

    print('No simulation output during the main planetary transit. The plotting script exits now.')
    sys.exit()


R_mean =  np.mean(R[:, SEL], axis = 1)

if Normalize:

    R_mean += (1. - np.max(R_mean)) # Normalize the transit depth such that the smallest flux decrease is shifted to 1

"""
Plot the spectrum and store the figure
"""


fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

ax.plot(wavelength, R_mean, color = 'black', linewidth = 2)

ax.set_xlabel(r'$\lambda\,[\AA]$')
ax.set_ylabel(r'$\Re$')

ax.minorticks_on()
ax.tick_params(which = 'both', direction = 'in', right = True, top = True)

ax.set_xlim(np.min(wavelength), np.max(wavelength))
ax.set_ylim(np.min(R_mean) - 0.05 * (1 - np.min(R_mean)), 1 + 0.05 * (1 - np.min(R_mean)))

plt.tight_layout()

plt.savefig(PARENTPATH + '/figures/' + paramsFilename + '_spectrumPlot.pdf', dpi = 150)
