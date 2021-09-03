"""
Author: Andrea Gebek
Created on 9.7.2021
Plot a transmission spectrum
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
Read in settings file and stored spectrum
"""

paramsFilename = sys.argv[1]

with open('../../' + paramsFilename + '.txt') as file:
    param = json.load(file)


benchmark = param['Output']['benchmark']

SpectrumFile = np.loadtxt('../../' + paramsFilename + '_spectrum.txt')

wavelength = SpectrumFile[:, 0] * 1e8 # In Angstrom
transit_depth = SpectrumFile[:, 1] + (1 - np.max(SpectrumFile[:, 1]))
print((1 - np.max(SpectrumFile[:, 1]))*100)
"""
Plot the spectrum and store the figure
"""


fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

ax.plot(wavelength, transit_depth, color = 'blue', linewidth = 1, label = 'Spectrum')

if benchmark:

    transit_depth_benchmark = SpectrumFile[:, 2] + (1 - np.max(SpectrumFile[:, 2]))
    ax.plot(wavelength, transit_depth_benchmark, color = 'red', linewidth = 1, label = 'Benchmark')


    lg = ax.legend(loc = 'lower center')
    lg.get_frame().set_linewidth(0)

ax.set_xlabel(r'$\lambda\,[\AA]$')
ax.set_ylabel(r'$\Re$')

ax.minorticks_on()
ax.tick_params(which = 'both', direction = 'in', right = True, top = True)

ax.set_xlim(np.min(wavelength), np.max(wavelength))
ax.set_ylim(np.min(transit_depth) - 0.05 * (1 - np.min(transit_depth)), 1 + 0.05 * (1 - np.min(transit_depth)))

plt.tight_layout()

plt.savefig('../../' + paramsFilename + '_spectrumPlot.pdf', dpi = 150)