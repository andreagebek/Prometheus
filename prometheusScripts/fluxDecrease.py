"""
Calculate the light curve of a transiting exoplanet.
Created on 15. July 2021 by Andrea Gebek.
"""

import numpy as np
from scipy.interpolate import interp1d
import sys
import os
from datetime import datetime
SCRIPTPATH = os.path.realpath(__file__)
GITPATH = os.path.dirname(os.path.dirname(SCRIPTPATH))
sys.path.append(GITPATH) 
import prometheusScripts.constants as const
import prometheusScripts.geometryHandler as geom
import prometheusScripts.gasProperties as gasprop
import prometheusScripts.stellarSpectrum as stellar

def constructAxis(gridsDict, architectureDict, axisName):

    if axisName == 'x':

        x_midpoint = gridsDict['x_midpoint']
        x_border = gridsDict['x_border']
        x_steps = gridsDict['x_steps']
        x_axis = np.linspace(x_midpoint - x_border, x_midpoint + x_border, int(x_steps) + 1)[:-1] + x_border / float(x_steps)  

        return x_axis

    elif axisName == 'rho':   

        rho_steps = gridsDict['rho_steps']
        rho_axis = np.linspace(0, architectureDict['R_star'], int(rho_steps) + 1)[:-1] + 0.5 * architectureDict['R_star'] / float(rho_steps)

        return rho_axis

    elif axisName == 'phi':

        phi_steps = gridsDict['phi_steps']
        phi_axis = np.linspace(0, 2 * np.pi, int(phi_steps) + 1)[:-1] + np.pi / float(phi_steps)

        return phi_axis

    elif axisName == 'orbphase':

        orbphase_border = gridsDict['orbphase_border']
        orbphase_axis = np.linspace(-orbphase_border, orbphase_border, int(gridsDict['orbphase_steps']))

        return orbphase_axis  

    elif axisName == 'wavelength':

        wavelength = np.arange(gridsDict['lower_w'], gridsDict['upper_w'], gridsDict['resolution'])

        return wavelength

def calculateCLV(rho, R_star, u1, u2):

    arg = 1. - np.sqrt(1. - rho**2 / R_star**2)

    return 1. - u1 * arg - u2 * arg**2

def calculateRM(phi, rho, wavelengthArray, architectureDict, Fstar_function):

    v_los = geom.calculateStellarLOSvelocity(architectureDict, phi, rho)
    w_shift = gasprop.calculateDopplerShift(v_los)

    F_shifted = Fstar_function(wavelengthArray / w_shift) # This contains the shifted flux incident on the exoplanet, instead of shifting the spectra
    # at all positions on the stellar surface just read in the spectrum at different wavelengths depending on the position in the spatial grid (phi and rho)

    return F_shifted

def getPHOENIX_output(wavelengthArray, fundamentalsDict, architectureDict):

    if fundamentalsDict['RM_effect']:

        T_starrot = architectureDict['period_starrot']
        R_star = architectureDict['R_star']

        PHOENIX_output = stellar.readSpectrum(architectureDict['T_eff'], architectureDict['log_g'], architectureDict['Fe_H'], architectureDict['alpha_Fe'])
    
        w_star = PHOENIX_output[0]
        v_max = 2. * np.pi * R_star / T_starrot
        w_max = np.max(wavelengthArray * gasprop.calculateDopplerShift(-v_max))
        w_min = np.min(wavelengthArray * gasprop.calculateDopplerShift(v_max))
        SEL_w = np.argwhere((w_star > w_min) * (w_star < w_max))[:, 0]

        SEL = np.concatenate((np.array([np.min(SEL_w) - 1]), SEL_w, np.array([np.max(SEL_w) + 1])))

        w_starSEL = w_star[SEL]
        F_0 = PHOENIX_output[1][SEL]

        Fstar_function = interp1d(w_starSEL, F_0, kind = 'cubic')

    else:

        Fstar_function = None

    return Fstar_function

def getFstarIntegrated(wavelengthArray, fundamentalsDict, architectureDict, gridsDict, PHOENIX_output):
    # Return F_star(wavelength, phi, rho) integrated over the stellar disk

    R_star = architectureDict['R_star']

    if not fundamentalsDict['CLV_variations'] and not fundamentalsDict['RM_effect']:

        FstarIntegrated = np.pi * R_star**2 * np.ones_like(wavelengthArray)


    elif fundamentalsDict['CLV_variations'] and not fundamentalsDict['RM_effect']:

        FstarIntegrated = np.pi * R_star**2 * (1. - architectureDict['u1'] / 3. - architectureDict['u2'] / 6.) * np.ones_like(wavelengthArray)


    else: # RM_effect = True

        phiArray = constructAxis(gridsDict, architectureDict, 'phi')
        rhoArray = constructAxis(gridsDict, architectureDict, 'rho')

        phiGrid, rhoGrid = np.meshgrid(phiArray, rhoArray, indexing = 'ij')

        GRID = np.stack((phiGrid.flatten(), rhoGrid.flatten()), axis = -1)

        delta_rho = R_star / float(gridsDict['rho_steps'])
        delta_phi = 2 * np.pi / float(gridsDict['phi_steps'])

        FstarIntegrated = np.zeros_like(wavelengthArray)

        for point in GRID:

            phi = point[0]
            rho = point[1]
        
            Fstar = calculateRM(phi, rho, wavelengthArray, architectureDict, PHOENIX_output)

            if fundamentalsDict['CLV_variations']:

                Fstar *= calculateCLV(rho, R_star, architectureDict['u1'], architectureDict['u2'])
            
            FstarIntegrated += Fstar * delta_phi * delta_rho * rho

    return FstarIntegrated

def getFstar(phi, rho, wavelengthArray, architectureDict, fundamentalsDict, PHOENIX_output):

    # Return F_star(wavelength, phi, rho)

    R_star = architectureDict['R_star']

    if not fundamentalsDict['CLV_variations'] and not fundamentalsDict['RM_effect']:

        Fstar = np.ones_like(wavelengthArray)

    elif fundamentalsDict['CLV_variations'] and not fundamentalsDict['RM_effect']:

        Fstar = calculateCLV(rho, R_star, architectureDict['u1'], architectureDict['u2']) * np.ones_like(wavelengthArray)

    else: # RM_effect = True

        Fstar = calculateRM(phi, rho, wavelengthArray, architectureDict, PHOENIX_output)

        if fundamentalsDict['CLV_variations']:

            Fstar *= calculateCLV(rho, R_star, architectureDict['u1'], architectureDict['u2'])

    return Fstar


def checkBlock(phi, rho, orbphase, architectureDict, fundamentalsDict):
    # Check if the chord is blocked by the planet or the moon

    y, z = geom.getCartesianFromCylinder(phi, rho)
    y_p = geom.getPlanetPosition(architectureDict, orbphase)[1]

    blockingPlanet = (np.sqrt((y - y_p)**2 + z**2) < architectureDict['R_0'])

    if blockingPlanet:

        return True
    
    if fundamentalsDict['ExomoonSource']:

        y_moon = geom.getMoonPosition(architectureDict, orbphase)[1]

        blockingMoon = ((y - y_moon)**2 + z**2 < architectureDict['R_moon']**2)

        if blockingMoon:

            return True
        
    return False


def calculateOpticalDepth(phi, rho, orbphase, xArray, wavelengthArray, fundamentalsDict, architectureDict, scenarioDict, speciesDict, gridsDict, sigmaLookupDict):

    delta_x = 2. * gridsDict['x_border'] / float(gridsDict['x_steps'])

    tau = 0.

    for key_scenario in scenarioDict.keys(): # Loop over scenarios

        specificScenarioDict = scenarioDict[key_scenario]

        n = gasprop.getNumberDensity(phi, rho, orbphase, xArray, key_scenario, specificScenarioDict, architectureDict, fundamentalsDict)     

        sigma_abs = gasprop.getAbsorptionCrossSection(phi, rho, orbphase, xArray, wavelengthArray, key_scenario, fundamentalsDict, specificScenarioDict, architectureDict, speciesDict, sigmaLookupDict)

        tau += delta_x * np.sum(np.multiply(sigma_abs, n), axis = 1)

    return tau # tau(lambda)


def evaluateChord(point, args): # Function to be multiprocessed

    phi, rho, orbphase = point

    delta_phi, delta_rho, xArray, wavelengthArray, architectureDict, fundamentalsDict, scenarioDict, speciesDict, gridsDict, outputDict, sigmaLookupDict, Fstar_function, FstarIntegrated = args

    if checkBlock(phi, rho, orbphase, architectureDict, fundamentalsDict):

        tau = np.inf * np.ones_like(wavelengthArray)

    else:

        tau = calculateOpticalDepth(phi, rho, orbphase, xArray, wavelengthArray, fundamentalsDict, architectureDict, scenarioDict, speciesDict, gridsDict, sigmaLookupDict)

    Fstar = getFstar(phi, rho, wavelengthArray, architectureDict, fundamentalsDict, Fstar_function)

    singleChord = rho * Fstar * np.exp(-tau) * delta_phi * delta_rho

    if outputDict['recordTau']:

        return np.float32(singleChord / FstarIntegrated), np.float32(tau)

    return np.float32(singleChord / FstarIntegrated), None

def prepareArguments(fundamentalsDict, architectureDict, scenarioDict, speciesDict, gridsDict, outputDict, startTime):

    xArray = constructAxis(gridsDict, architectureDict, 'x')

    wavelengthArray = constructAxis(gridsDict, architectureDict, 'wavelength')

    phiArray = constructAxis(gridsDict, architectureDict, 'phi')
    rhoArray = constructAxis(gridsDict, architectureDict, 'rho')
    orbphaseArray = constructAxis(gridsDict, architectureDict, 'orbphase')

    phiGrid, rhoGrid, orbphaseGrid = np.meshgrid(phiArray, rhoArray, orbphaseArray, indexing = 'ij')

    GRID = np.stack((phiGrid.flatten(), rhoGrid.flatten(), orbphaseGrid.flatten()), axis = -1)

    delta_rho = architectureDict['R_star'] / float(gridsDict['rho_steps'])
    delta_phi = 2 * np.pi / float(gridsDict['phi_steps'])

    Fstar_function = getPHOENIX_output(wavelengthArray, fundamentalsDict, architectureDict)

    print('\nPHOENIX spectrum downloaded (if required):', datetime.now() - startTime)

    FstarIntegrated = getFstarIntegrated(wavelengthArray, fundamentalsDict, architectureDict, gridsDict, Fstar_function) # Stellar flux as a function of wavelength, integrated over the stellar disk

    print('Integrated stellar flux calculated:', datetime.now() - startTime)

    sigmaLookupDict = gasprop.createLookupAbsorption(xArray, wavelengthArray, GRID, fundamentalsDict, architectureDict, scenarioDict, speciesDict)

    print('Lookup absorption cross section calculated (if required):', datetime.now() - startTime)

    args = (delta_phi, delta_rho, xArray, wavelengthArray, architectureDict, fundamentalsDict, scenarioDict, speciesDict, gridsDict, outputDict, sigmaLookupDict, Fstar_function, FstarIntegrated)

    return GRID, args