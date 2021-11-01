"""
Author: Andrea Gebek
Created on 29.10.2021
Plot the spatial grid.
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
import eliteScripts.geometryHandler as geom
import eliteScripts.fluxDecrease as flux

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

plotStar = True
plotPlanets = True
plotRhoArrow = True
plotPhi = True
plotXArrow = True
writeRadii = True
plotXwobble = True

paramsFilename = sys.argv[1]

with open(PARENTPATH + '/setupFiles/' + paramsFilename + '.txt') as file:
    param = json.load(file)

architectureDict = param['Architecture']
gridsDict = param['Grids']

R_0 = architectureDict['R_0']
R_star = architectureDict['R_star']
a_p = architectureDict['a_p']

x_border = gridsDict['x_border']
x_steps = gridsDict['x_steps']
phi_steps = gridsDict['phi_steps']
rho_steps = gridsDict['rho_steps']
orbphase_border = gridsDict['orbphase_border']
orbphase_steps = gridsDict['orbphase_steps']

x_axis = flux.constructAxis(gridsDict, architectureDict, 'x')
rho_axis = flux.constructAxis(gridsDict, architectureDict, 'rho')
phi_axis = flux.constructAxis(gridsDict, architectureDict, 'phi')
orbphase_axis = flux.constructAxis(gridsDict, architectureDict, 'orbphase') # In radians



fig, axes = plt.subplots(figsize=(18, 8), nrows = 1, ncols = 2)

"""
Plot birds-eye view
"""

x, y = np.meshgrid(x_axis, np.concatenate((-rho_axis, rho_axis)), indexing = 'ij')

rectangle = matplotlib.patches.Rectangle(((a_p - x_border) / R_0, -R_star / R_0), 2 * x_border / R_0, 2 * R_star / R_0, linewidth = 1, edgecolor = 'lightblue', facecolor = 'None')
axes[0].add_patch(rectangle)

if plotStar:
    starCircle = plt.Circle((0, 0), R_star / R_0, facecolor = 'orange', edgecolor = 'black', fill = True, alpha = 0.25, linewidth = 2)
    axes[0].add_patch(starCircle)
 
if plotPlanets:

    x_p, y_p = geom.getPlanetPosition(architectureDict, orbphase_axis)

    for idx in range(int(orbphase_steps)):

        planetCircle = plt.Circle((x_p[idx] / R_0, y_p[idx] / R_0), 1, color = 'black', fill = True, linewidth = 1, alpha = 0.3)
        axes[0].add_patch(planetCircle)

if plotXArrow:
    axes[0].arrow(0, 0, 1.25 * a_p / R_0, 0, length_includes_head = True, head_width = 2.5, linewidth = 1, facecolor = 'black')
    axes[0].annotate(r'$x$', xy = (0.4 * a_p / R_0, 0.03 * a_p / R_0))  

if writeRadii:
    axes[0].annotate(r'$R_{\ast}$', xy = (-2.5 * R_star / R_0, R_star / R_0))
    axes[0].annotate(r'$a_p$', xy = (-0.82 * a_p / R_0, -0.82 * a_p / R_0))   

if plotXwobble:
    axes[0].arrow((a_p - x_border)/ R_0, 0, -0.15 * a_p / R_0, 0, length_includes_head = True, head_width = 2.5, width = 1, linewidth = 1, facecolor = 'lightblue')

orbitCircle = plt.Circle((0, 0), a_p / R_0, fill = False, linewidth = 2, alpha = 0.5)
axes[0].add_patch(orbitCircle)

axes[0].scatter((a_p + x) / R_0, y / R_0, color = 'blue', s = 0.5)

axes[0].set_xlabel(r'$x\,[R_0]$')
axes[0].set_ylabel(r'$y\,[R_0]$')

axes[0].set_xlim(-1.3 * a_p / R_0, 1.3 * a_p / R_0)
axes[0].set_ylim(-1.3 * a_p / R_0, 1.3 * a_p / R_0)

"""
Plot face-on view
"""

y = np.tensordot(rho_axis, np.sin(phi_axis), axes = 0)
z = np.tensordot(rho_axis, np.cos(phi_axis), axes = 0)

axes[1].scatter(y / R_0, z / R_0, color = 'blue', s = 1.5)



if plotPlanets:

    y_p = geom.getPlanetPosition(architectureDict, orbphase_axis)[1]

    for y in y_p:

        planetCircle = plt.Circle((y / R_0, 0), 1, color = 'black', fill = True, alpha = 0.5, linewidth = 1)
        axes[1].add_patch(planetCircle)

if plotStar:
    starCircle = plt.Circle((0, 0), R_star / R_0, facecolor = 'orange', edgecolor = 'black', fill = True, alpha = 0.25, linewidth = 2)
    axes[1].add_patch(starCircle)

if plotRhoArrow:
    plt.arrow(0, 0, 0, R_star / R_0, length_includes_head = True, head_width = 0.2, linewidth = 0.3, facecolor = 'black')
    axes[1].annotate(r'$\varrho$', xy = (0.04 * R_star / R_0, 0.7 * R_star / R_0))

if writeRadii:
    axes[1].annotate(r'$R_{\ast}$', xy = (-0.8 * R_star / R_0, -0.8 * R_star / R_0))
    axes[1].annotate(r'$R_0$', xy = (0.9, -1.6))

if plotPhi:
    arc = matplotlib.patches.Arc((0, 0), 0.6 * R_star / R_0, 0.6 * R_star / R_0, theta1 = 63.43, theta2 = 90, linewidth = 2)
    axes[1].add_patch(arc)
    plt.arrow(0, 0, 0.15 * R_star / R_0, 0.3 * R_star / R_0, linewidth = 0.3)
    axes[1].annotate(r'$\varphi$', xy = (0.04 * R_star / R_0, 0.35 * R_star / R_0))
 

axes[1].set_xlabel(r'$y\,[R_0]$')
axes[1].set_ylabel(r'$z\,[R_0]$')

axes[1].set_xlim(-1.05 * R_star / R_0, 1.05 * R_star / R_0)
axes[1].set_ylim(-1.05 * R_star / R_0, 1.05 * R_star / R_0)



for ax in axes:

    ax.set_aspect('equal', adjustable='box')

    ax.minorticks_on()
    ax.tick_params(which = 'both', direction = 'in', right = True, top = True)

plt.tight_layout()

plt.savefig(PARENTPATH + '/figures/' + paramsFilename + '_SpatialGrids.pdf', dpi = 50)

