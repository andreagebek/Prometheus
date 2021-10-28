"""
Main code to run the radiative transfer calculation.
Created on 15. July 2021 by Andrea Gebek.
"""

import sys
import json
import numpy as np
from datetime import datetime
import eliteScripts.constants as const

N_arguments = len(sys.argv)

if N_arguments == 1:

    import eliteScripts.setup

    sys.exit(0)

import eliteScripts.fluxDecrease as flux

startTime = datetime.now()

"""
Read in the json parameters and perform the radiative transfer calculation
"""

paramsFilename = sys.argv[1]

with open('../SetupFiles/' + paramsFilename + '.txt') as file:
    param = json.load(file)

fundamentalsDict = param['Fundamentals']
scenarioDict = param['Scenarios']
architectureDict = param['Architecture']
speciesDict = param['Species']
gridsDict = param['Grids']
outputDict = param['Output']

resultsDict = flux.calculateTransitDepth(fundamentalsDict, architectureDict, scenarioDict, speciesDict, gridsDict, outputDict)
R = resultsDict['R'].flatten()

"""
Store the output in .txt files
"""

orbphase_axis = np.linspace(-gridsDict['orbphase_border'], gridsDict['orbphase_border'], int(gridsDict['orbphase_steps']))
wavelength_axis = np.arange(gridsDict['lower_w'], gridsDict['upper_w'], gridsDict['resolution']) * 1e8 # Conversion from cm to Angstrom

wavelength, orbphase = np.meshgrid(wavelength_axis, orbphase_axis, indexing = 'ij')
wavelength = wavelength.flatten()
orbphase = orbphase.flatten()


header = 'Wavelength grid (Å), Orbital phase grid, R'

np.savetxt('../output/' + paramsFilename + '_lightcurve.txt', np.array([wavelength, orbphase, R]).T, header = header)


if outputDict['benchmark']:

    R_benchmark = resultsDict['R_benchmark'].flatten()

    np.savetxt('../output/' + paramsFilename + '_barometricBenchmark.txt', np.array([wavelength, orbphase, R_benchmark]).T, header = header)


if outputDict['recordTau']:

    rho_axis = np.linspace(0, architectureDict['R_star'], int(gridsDict['rho_steps']) + 1)[:-1] + 0.5 * architectureDict['R_star'] / float(gridsDict['rho_steps'])

    phi_axis = np.linspace(0, 2 * np.pi, int(gridsDict['phi_steps']) + 1)[:-1] + np.pi / float(gridsDict['phi_steps'])

    phi, rho = np.meshgrid(phi_axis, rho_axis, indexing = 'ij')
    phi = phi.flatten()
    rho = rho.flatten()

    tauDisk = resultsDict['tauDisk'].flatten()

    np.savetxt('../output/' + paramsFilename + '_tau.txt', np.array([phi, rho, tauDisk]).T, header = 'phi grid [rad], rho grid [cm], tau')  

elapsedTime = datetime.now() - startTime

print("DISHOOM-ELITE finished, yay! Elapsed time is:", elapsedTime)

print("The maximal flux decrease due to atmospheric/exospheric absorption in percent is:", np.abs(np.round(100 * (1 - np.min(R)), 5)))

print("The minimal flux decrease due to atmospheric/exospheric absorption in percent is:", np.abs(np.round(100 * (1 - np.max(R)), 5)))