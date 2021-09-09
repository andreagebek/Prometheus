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
Calculate the transmission spectrum by averaging
the depth in the lightcurve during full occultation
of the planet's R_0.
"""

paramsFilename = sys.argv[1]

with open('../../' + paramsFilename + '.txt') as file:
    param = json.load(file)

architecture_dict = param['Architecture']

R_s = architecture_dict['R_star']
R_0 = architecture_dict['R_0']
a_p = architecture_dict['a_p']

grids_dict = param['Grids']

wavelength = np.arange(grids_dict['lower_w'], grids_dict['upper_w'], grids_dict['resolution']) * 1e8 # In Angstrom

LightcurveFile = np.loadtxt('../../' + paramsFilename + '_lightcurve.txt')

orbphase = LightcurveFile[:, 0]  # Not in rad, ranges from -0.5 to 0.5

orbphaseFullIngress = np.arcsin((R_s - R_0) / a_p) / (2. * np.pi) # Orbital phase at which the planet's R_0 is fully within the stellar disk

SEL = np.abs(orbphase) <= orbphaseFullIngress

transit_depth = []
for idx in range(len(wavelength)):

    lightcurve = LightcurveFile[:, idx + 1][SEL]
    transit_depth.append(np.mean(lightcurve))

transit_depth = np.array(transit_depth) + (1 - np.max(transit_depth)) # Normalize the transit depth such that the smallest flux decrease is shifted to 1

"""
Plot the spectrum and store the figure
"""


fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

ax.plot(wavelength, transit_depth, color = 'blue', linewidth = 1)


ax.set_xlabel(r'$\lambda\,[\AA]$')
ax.set_ylabel(r'$\Re$')

ax.minorticks_on()
ax.tick_params(which = 'both', direction = 'in', right = True, top = True)

ax.set_xlim(np.min(wavelength), np.max(wavelength))
ax.set_ylim(np.min(transit_depth) - 0.05 * (1 - np.min(transit_depth)), 1 + 0.05 * (1 - np.min(transit_depth)))

plt.tight_layout()

plt.savefig('../../' + paramsFilename + '_spectrumFromLightcurvePlot.pdf', dpi = 150)
