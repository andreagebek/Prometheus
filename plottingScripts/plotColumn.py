"""
Author: Andrea Gebek
Created on 4.12.2021
Plot the total line-of-sight column density
(between star and observer) of each absorbing species.
If there is a time-dependency (i.e. multiple orbital
phases), create an mp4 movie. Note that for mp4
movies you need the ffmpeg program.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import json
import sys
import os
from mpl_toolkits import axes_grid1
import matplotlib.animation as animation
SCRIPTPATH = os.path.realpath(__file__)
GITPATH = os.path.dirname(os.path.dirname(SCRIPTPATH))
PARENTPATH = os.path.dirname(GITPATH)
sys.path.append(GITPATH) 
import prometheusScripts.geometryHandler as geom
import prometheusScripts.gasProperties as gasprop
import prometheusScripts.fluxDecrease as flux
import datetime

startTime = datetime.datetime.now()

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
Prepare quantities for the plots
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

"""
Plot single figure if there is no time dependency,
or save an mp4 movie for each species if there is a time dependency
"""

if gridsDict['orbphase_steps'] == 1: 

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



else:
    
    for idx1 in range(N_species):

        fig = plt.figure(figsize = (10, 8))
        ax = fig.add_subplot(111)

        ims = []

        max_column = np.max(np.ma.masked_invalid(N_speciesDict[speciesList[idx1]])) # Ignore inf values

        for idx2 in range(gridsDict['orbphase_steps']):

            column = N_speciesDict[speciesList[idx1]][:, :, idx2].flatten()

            with np.errstate(divide = 'ignore'):
                im = ax.scatter(y / R_0, z / R_0, c = np.log10(column), vmin = 5, vmax = np.log10(max_column), s = 2, cmap = 'Spectral_r', animated = True)

            if idx2 == 0: # Else the colorbar is added to every single frame and the program takes FOREVER to execute
                cbar = add_colorbar(im, ax = ax)
                cbar.set_label(r'$\log_{10}(N_\mathrm{los})\,[$' + speciesList[idx1] + r'$\,\mathrm{cm}^{-2}]$')
                cbar.ax.minorticks_on()

            ax.set_xlabel(r'$y\,[R_0]$')
            ax.set_ylabel(r'$z\,[R_0]$')

            ax.set_xlim(-1.05 * R_star / R_0, 1.05 * R_star / R_0)
            ax.set_ylim(-1.05 * R_star / R_0, 1.05 * R_star / R_0)

            ax.set_aspect('equal', adjustable='box')

            ax.minorticks_on()
            ax.tick_params(which = 'both', direction = 'in', right = True, top = True)

            ims.append([im])

        ani = animation.ArtistAnimation(fig, ims, interval = 50, blit = True, repeat_delay = 1000)
        writervideo = animation.FFMpegWriter(fps=60)
        ani.save(PARENTPATH + '/figures/' + paramsFilename + '_' + speciesList[idx1] + 'columnMovie.mp4', writer = writervideo)
