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
Read in settings file and stored lightcurve.
Calculate the transmission spectrum by averaging
the depth in the lightcurve during full occultation
of the planet's R_0.
"""

plotMeanSpectrum = True # Plot the spectrum averaged over all orbital phases between planetary ingress and egress
plotSpectra = True  # Plot the spectra at all orbital phases
plotBenchmarkSpectra = True # Plot barometric benchmark spectra at all orbital phases (only if the barometric benchmark option is True in the setup file)
Normalize = False # Do a quick normalization of the spectra

paramsFilename = sys.argv[1]

with open(PARENTPATH + '/setupFiles/' + paramsFilename + '.txt') as file:
    param = json.load(file)

architectureDict = param['Architecture']
gridsDict = param['Grids']

orbphase = flux.constructAxis(gridsDict, architectureDict, 'orbphase') # In radians
wavelength = flux.constructAxis(gridsDict, architectureDict, 'wavelength') # In cm

R_star = architectureDict['R_star']
R_0 = architectureDict['R_0']
a_p = architectureDict['a_p']

benchmark = param['Output']['benchmark']

LightcurveFile = np.loadtxt(PARENTPATH + '/output/' + paramsFilename + '_lightcurve.txt')

R = LightcurveFile[:, 2].reshape(len(wavelength), len(orbphase))

R_plot = []
for idx in range(len(orbphase)):

    if Normalize:
        R_plot.append(R[:, idx] + (1 - np.max(R[:, idx]))) # Normalize the transit depth such that the smallest flux decrease is shifted to 1
    else:
        R_plot.append(R[:, idx])

if benchmark:

    R_benchmark = np.loadtxt(PARENTPATH + '/output/' + paramsFilename + '_barometricBenchmark.txt')[:, 2].reshape(len(wavelength), len(orbphase))

    R_benchmark_plot = []
    for idx in range(len(orbphase)):

        if Normalize:
            R_benchmark_plot.append(R_benchmark[:, idx] + (1 - np.max(R_benchmark[:, idx]))) # Normalize the transit depth such that the smallest flux decrease is shifted to 1
        else:
            R_benchmark_plot.append(R_benchmark[:, idx])
"""
Plot the spectrum and store the figure
"""


fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

cmap = matplotlib.cm.get_cmap('viridis')

for idx, r in enumerate(R_plot):

    if plotSpectra:
        ax.plot(wavelength * 1e8, r, color = cmap(float(idx) / float(len(orbphase))), linewidth = 1, label = r'$\phi=$' + str(np.round(orbphase[idx] / (2. * np.pi), 3)))

    if benchmark and plotBenchmarkSpectra:

        ax.plot(wavelength * 1e8, R_benchmark_plot[idx], color = cmap(float(idx) / float(len(orbphase))), linewidth = 1, linestyle = '--')

orbphaseFullIngress = np.arcsin((R_star - R_0) / a_p) # Orbital phase at which the planet's R_0 is fully within the stellar disk
SEL = np.abs(orbphase) <= orbphaseFullIngress

if len(orbphase[SEL]) > 1 and plotMeanSpectrum:

    R_avg = np.mean(R[:, SEL], axis = 1)

    if Normalize:
        R_avg += 1 - np.max(R_avg)

    ax.plot(wavelength * 1e8, R_avg, color = 'black', linewidth = 2, label = 'Mean')




lg = ax.legend(loc = 'center right')
lg.get_frame().set_alpha(0)
lg.get_frame().set_linewidth(0)

ax.set_xlabel(r'$\lambda\,[\AA]$')
ax.set_ylabel(r'$\Re$')

ax.minorticks_on()
ax.tick_params(which = 'both', direction = 'in', right = True, top = True)

ax.set_xlim(np.min(wavelength) * 1e8, np.max(wavelength) * 1e8)
ax.set_ylim(np.min(R_plot) - 0.05 * (1 - np.min(R_plot)), 1 + 0.05 * (1 - np.min(R_plot)))

plt.tight_layout()

plt.savefig(PARENTPATH + '/figures/' + paramsFilename + '_spectrumPlot.pdf', dpi = 150)
