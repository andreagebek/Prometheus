"""
Author: Andrea Gebek
Created on 18.8.2021
Plot the absorption cross section
for the selected absorbers.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import json
import sys
import os
from scipy.interpolate import RegularGridInterpolator
SCRIPTPATH = os.path.realpath(__file__)
GITPATH = os.path.dirname(os.path.dirname(SCRIPTPATH))
PARENTPATH = os.path.dirname(GITPATH)
sys.path.append(GITPATH) 
import prometheusScripts.gasProperties as gasprop
import prometheusScripts.fluxDecrease as flux
import prometheusScripts.constants as const


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
Read in settings file and calculate the absorption cross section
"""

plotRayleighScatt = True

T = 1000 # Temperature in K to get velocity dispersion
P = 1e0 # Pressure in cgs units (for molecular absorption), 1 bar is 1e6 barye (cgs-unit) 

paramsFilename = sys.argv[1]

with open(PARENTPATH + '/setupFiles/' + paramsFilename + '.txt') as file:
    param = json.load(file)

architectureDict = param['Architecture']
gridsDict = param['Grids']
speciesDict = param['Species']

wavelength = flux.constructAxis(gridsDict, architectureDict, 'wavelength')

sigmaAbsSpecies = []
evaporativeSpeciesList = [] # Atoms/Ions
molecularSpeciesList = [] # Molecules

for key_scenario in speciesDict.keys():

    for key_species in speciesDict[key_scenario].keys():

        if not key_species in evaporativeSpeciesList and 'sigma_v' in speciesDict[key_scenario][key_species].keys():

            evaporativeSpeciesList.append(key_species)

        if not key_species in molecularSpeciesList and not 'sigma_v' in speciesDict[key_scenario][key_species].keys():

            molecularSpeciesList.append(key_species)

labelList = evaporativeSpeciesList + molecularSpeciesList

for key_species in evaporativeSpeciesList:

    sigma_v = np.sqrt(const.k_B * T / const.speciesInfoDict[key_species][2])

    line_wavelength, line_gamma, line_f = gasprop.readLineList(key_species, wavelength)
    sigmaAbsSpecies.append(gasprop.calculateLineAbsorption(wavelength, line_wavelength, line_gamma, line_f, {'chi': 1., 'sigma_v': sigma_v}))

    labelList.append(key_species)

for key_species in molecularSpeciesList:

    P_mol, T_mol, wavelength_mol, sigma_mol = gasprop.readMolecularAbsorption(key_species)

    sigma_mol_function = RegularGridInterpolator((P_mol, T_mol, wavelength_mol[::-1]), sigma_mol[:, :, ::-1], bounds_error = False, fill_value = 0.)

    sigma_mol = sigma_mol_function((P, T, wavelength))

    sigmaAbsSpecies.append(sigma_mol)

"""
Plot the absorption cross section and store the figure
"""

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

cmap = matplotlib.cm.get_cmap('viridis')

for idx, sigmaAbs in enumerate(sigmaAbsSpecies):

    ax.plot(wavelength * 1e8, sigmaAbs, color = cmap(float(idx) / float(len(labelList))), linewidth = 2, label = labelList[idx])


if plotRayleighScatt:

    sigmaRayleigh = 8.49e-45 / wavelength**4

    ax.plot(wavelength * 1e8, sigmaRayleigh, color = 'black', linewidth = 2, linestyle = '--', label = 'Rayleigh')

lg = ax.legend(loc = 'upper right')
lg.get_frame().set_alpha(0)
lg.get_frame().set_linewidth(0)

ax.set_xlabel(r'$\lambda\,[\AA]$')
ax.set_ylabel(r'$\sigma_{\mathrm{abs}}\,[\mathrm{cm^2/particle}]$')

ax.minorticks_on()
ax.tick_params(which = 'both', direction = 'in', right = True, top = True)

ax.set_xlim(min(wavelength) * 1e8, max(wavelength) * 1e8)
ax.set_yscale('log')

plt.tight_layout()

plt.savefig(PARENTPATH + '/figures/' + paramsFilename + '_SigmaAbsPlot.pdf', dpi = 150)