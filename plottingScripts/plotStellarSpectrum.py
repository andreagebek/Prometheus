# coding=utf-8
"""
Author: Andrea Gebek
Created on 18.8.2021
Plot the Phoenix input 
spectrum from the host star.
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
import prometheusScripts.gasProperties as gasprop
import prometheusScripts.stellarSpectrum as stellar
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

plotRMDopplerShift = False

"""
Read in settings file and calculate stellar spectrum
"""

paramsFilename = sys.argv[1]

with open(PARENTPATH + '/setupFiles/' + paramsFilename + '.txt') as file:
    param = json.load(file)

architectureDict = param['Architecture']
gridsDict = param['Grids']
scenarioDict = param['Scenarios']
speciesDict = param['Species']

R_star = architectureDict['R_star']
T_starrot = architectureDict['period_starrot']


wavelength = flux.constructWavelengthGrid(gridsDict, scenarioDict, speciesDict)

T_eff = architectureDict['T_eff']
log_g = architectureDict['log_g']
Fe_H = architectureDict['Fe_H']
alpha_Fe = architectureDict['alpha_Fe']

PHOENIX_output = stellar.readSpectrum(T_eff, log_g, Fe_H, alpha_Fe)

w_star = PHOENIX_output[0]
v_max = 2. * np.pi * R_star / T_starrot
w_max = np.max(wavelength * gasprop.calculateDopplerShift(-v_max))
w_min = np.min(wavelength * gasprop.calculateDopplerShift(v_max))
SEL_w = np.argwhere((w_star > w_min) * (w_star < w_max))[:, 0]

SEL = np.concatenate((np.array([np.min(SEL_w) - 1]), SEL_w, np.array([np.max(SEL_w) + 1])))

w_star = w_star[SEL]
F_0 = PHOENIX_output[1][SEL]



"""
Plot the spectrum and store the figure
"""


fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

ax.plot(w_star * 1e8, F_0, color = 'blue', linewidth = 2)


if plotRMDopplerShift:
    delta_w_max = (gasprop.calculateDopplerShift(v_max) - 1.) * np.mean(wavelength) * 1e8
    plt.arrow(np.mean(wavelength) * 1e8, np.min(F_0), delta_w_max, 0, linewidth = 4, facecolor = 'black')
    ax.annotate(r'$\Delta \lambda$', xy = (np.mean(wavelength) * 1e8, 1.5 * np.min(F_0)), ha = 'center')  

ax.set_xlabel(r'$\lambda\,[\AA]$')
ax.set_ylabel(r'$F_{\lambda}\,[\rm{erg\,s^{-1}\,cm^{-2}\,cm^{-1}}]$')

ax.minorticks_on()
ax.tick_params(which = 'both', direction = 'in', right = True, top = True)

ax.set_xlim(np.min(w_star * 1e8), np.max(w_star * 1e8))


plt.tight_layout()

plt.savefig(PARENTPATH + '/figures/' + paramsFilename + '_stellarSpectrumPlot.pdf', dpi = 150)