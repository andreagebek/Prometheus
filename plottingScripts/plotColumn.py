"""
Author: Andrea Gebek
Created on 4.12.2021
Plot the total line-of-sight column density
(between star and observer) of each absorbing species.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import json
import sys
import os
from mpl_toolkits import axes_grid1
SCRIPTPATH = os.path.realpath(__file__)
GITPATH = os.path.dirname(os.path.dirname(SCRIPTPATH))
PARENTPATH = os.path.dirname(GITPATH)
sys.path.append(GITPATH) 
import prometheusScripts.geometryHandler as geom
import prometheusScripts.gasProperties as gasprop
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
Read in settings file
"""

plotPlanet = False
plotStar = True

paramsFilename = sys.argv[1]

with open(PARENTPATH + '/setupFiles/' + paramsFilename + '.txt') as file:
    param = json.load(file)

architectureDict = param['Architecture']
scenarioDict = param['Scenarios']
speciesDict = param['Species']
gridsDict = param['Grids']
fundamentalsDict = param['Fundamentals']

"""
Calculate the number density of each species by summing it over all scenarios
"""

x, phi, rho, orbphase = flux.constructSpatialGrid(gridsDict, architectureDict)

n_speciesDict = {}

for key_scenario in scenarioDict.keys():

    n_total = gasprop.getNumberDensity(x, phi, rho, orbphase, key_scenario, scenarioDict[key_scenario], architectureDict, fundamentalsDict)

    for key_species in speciesDict[key_scenario].keys():

        n_species = n_total * speciesDict[key_scenario][key_species]['chi']

        if key_species in n_speciesDict:
            
            n_speciesDict[key_species] += n_species

        else:

            n_speciesDict[key_species] = n_species

"""
Calculate the line-of-sight column by integrating along the x-axis
"""

delta_x = 2 * gridsDict['x_border'] / float(gridsDict['x_steps'])

N_species = len(n_speciesDict.keys())
N_speciesDict = {}
speciesList = []

for key_species in n_speciesDict.keys():

    N_speciesDict[key_species] = delta_x * np.sum(n_speciesDict[key_species], axis = 0)
    speciesList.append(key_species)


"""
Plot the column for each species in a different panel
"""

def add_colorbar(im, aspect=20, pad_fraction=0.5, **kwargs): # Add colorbar to a plot which matches the graph size
    """Add a vertical color bar to an image plot."""
    divider = axes_grid1.make_axes_locatable(im.axes)
    width = axes_grid1.axes_size.AxesY(im.axes, aspect=1./aspect)
    pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
    current_ax = plt.gca()
    cax = divider.append_axes("right", size=width, pad=pad)
    plt.sca(current_ax)
    return im.axes.figure.colorbar(im, cax=cax, **kwargs)

R_0 = architectureDict['R_0']
R_star = architectureDict['R_star']

y, z = geom.getCartesianFromCylinder(0, phi, rho)[1:] # x is not used here

y = y[0, :, :, 0].flatten() # Create a 1D-array with all the different y-coordinates of the 2D-grid (ignore x-coordinate and orbital phase)
z = z[0, :, :, 0].flatten() # Same for z-coordinate

fig, axes = plt.subplots(figsize=(2 + 10 * N_species, 8), nrows = 1, ncols = N_species)

for idx, ax in enumerate(axes):

    column = N_speciesDict[speciesList[idx]].flatten()


    with np.errstate(divide = 'ignore'):
        im = ax.scatter(y / R_0, z / R_0, c = np.log10(column), vmin = 5, s = 2, cmap = 'Spectral_r')
    
    if plotPlanet:
        planetCircle = plt.Circle((0, 0), 1, color = 'black', linewidth = 0)
        ax.add_patch(planetCircle)

    if plotStar:
        starCircle = plt.Circle((0, 0), R_star / R_0, color = 'black', fill = False, linewidth = 1)
        ax.add_patch(starCircle)

    cbar = add_colorbar(im, ax = ax)
    cbar.set_label(r'$\log_{10}(N_\mathrm{los})\,[$' + speciesList[idx] + r'$\,\mathrm{cm}^{-2}]$')
    cbar.ax.minorticks_on()

    ax.set_xlabel(r'$y\,[R_0]$')
    ax.set_ylabel(r'$z\,[R_0]$')

    ax.set_xlim(-1.05 * R_star / R_0, 1.05 * R_star / R_0)
    ax.set_ylim(-1.05 * R_star / R_0, 1.05 * R_star / R_0)

    ax.set_aspect('equal', adjustable='box')

    ax.minorticks_on()
    ax.tick_params(which = 'both', direction = 'in', right = True, top = True)

plt.tight_layout()

plt.savefig(PARENTPATH + '/figures/' + paramsFilename + '_columnPlot.pdf', dpi = 50)