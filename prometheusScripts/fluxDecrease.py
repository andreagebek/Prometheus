# coding=utf-8
"""
Calculate the light curve of a transiting exoplanet.
Created on 15. July 2021 by Andrea Gebek.
"""

import numpy as np
from scipy.interpolate import interp1d
from datetime import datetime
import prometheusScripts.geometryHandler as geom
import prometheusScripts.gasProperties as gasprop
import prometheusScripts.stellarSpectrum as stellar

def constructAxis(gridsDict, axisName): # For spatial axes, return both the grid and the spacing (difference between grid points)

    if axisName == 'x':

        x_midpoint = gridsDict['x_midpoint']
        x_border = gridsDict['x_border']
        x_steps = gridsDict['x_steps']
        x_axis = np.linspace(x_midpoint - x_border, x_midpoint + x_border, int(x_steps), endpoint = False, retstep = True)
        list(x_axis)[0] += x_border / float(x_steps)

        return x_axis

    elif axisName == 'rho':   

        rho_steps = gridsDict['rho_steps']
        upper_rho = gridsDict['upper_rho']
        rho_axis = np.linspace(0, upper_rho, int(rho_steps), endpoint = False, retstep = True)
        list(rho_axis)[0] += 0.5 * upper_rho / float(rho_steps)

        return rho_axis

    elif axisName == 'phi':

        phi_steps = gridsDict['phi_steps']
        phi_axis = np.linspace(0, 2 * np.pi, int(phi_steps), endpoint = False, retstep = True)
        list(phi_axis)[0] += np.pi / float(phi_steps)

        return phi_axis

    elif axisName == 'orbphase':

        orbphase_border = gridsDict['orbphase_border']
        orbphase_axis = np.linspace(-orbphase_border, orbphase_border, int(gridsDict['orbphase_steps']))

        return orbphase_axis  

def constructWavelengthGrid(gridsDict, scenarioDict, speciesDict):

    lower_w = gridsDict['lower_w']
    upper_w = gridsDict['upper_w']
    widthHighRes = gridsDict['widthHighRes']
    resolutionLow = gridsDict['resolutionLow']
    resolutionHigh = gridsDict['resolutionHigh']

    linesList = []

    for key_scenario in scenarioDict.keys():
        for key_species in speciesDict[key_scenario].keys():

            lines_w = gasprop.readLineList(key_species, np.array([lower_w, upper_w]))[0] # Read in center wavelengths of all absorption lines
            
            linesList.extend(lines_w)

    peaks = np.sort(np.unique(linesList))
    diff = np.concatenate(([np.inf], np.diff(peaks), [np.inf])) # Difference between each absorption line center. Add infinity at beginning and end such that
    # when calculating the high-resolution bands we keep the lines with the lowest/highest wavlengths.

    HighResBorders = ([], [])

    for idx, peak in enumerate(peaks):

        if diff[idx] > widthHighRes:

            HighResBorders[0].append(peak - widthHighRes / 2.)

        if diff[idx + 1] > widthHighRes:

            HighResBorders[1].append(peak + widthHighRes / 2.)

    grid = []

    for idx in range(len(HighResBorders[0])):

        grid.append(np.arange(HighResBorders[0][idx], HighResBorders[1][idx], resolutionHigh))

        if idx == 0 and lower_w < HighResBorders[0][0]:
            grid.append(np.arange(lower_w, HighResBorders[0][0], resolutionLow))

        elif idx == len(HighResBorders[0]) - 1 and upper_w > HighResBorders[1][-1]:
            grid.append(np.arange(HighResBorders[1][-1], upper_w, resolutionLow))
            grid.append(np.arange(HighResBorders[1][idx - 1], HighResBorders[0][idx], resolutionLow))
            
        else:
            grid.append(np.arange(HighResBorders[1][idx - 1], HighResBorders[0][idx], resolutionLow))

    wavelength = np.sort(np.concatenate(grid))

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
    # Return F_star(wavelength, phi, rho) integrated over the stellar disk and over the upper part of the stellar disk (if upper_rho is not the stellar radius)

    R_star = architectureDict['R_star']

    if not fundamentalsDict['CLV_variations'] and not fundamentalsDict['RM_effect']:

        FstarIntegrated = np.pi * R_star**2 * np.ones_like(wavelengthArray)

        FstarUpper = np.pi * (R_star**2 - gridsDict['upper_rho']**2) * np.ones_like(wavelengthArray)


    elif fundamentalsDict['CLV_variations'] and not fundamentalsDict['RM_effect']:

        u1 = architectureDict['u1']
        u2 = architectureDict['u2']
        rho_u = gridsDict['upper_rho']

        FstarIntegrated = np.pi * R_star**2 * (1. - u1 / 3. - u2 / 6.) * np.ones_like(wavelengthArray)


        upperTerm = 0.5 * (-u2 * R_star**2 - u1 * R_star**2 + R_star**2)
        term1 = -4. * R_star**2 * u1 * (1. - rho_u**2 / R_star**2)**1.5
        term2 = R_star**2 * u2 * (6 * rho_u**2 / R_star**2 + 8. * (1. - rho_u**2 / R_star**2)**1.5 - 3. * (R_star**2 - rho_u**2)**2 / R_star**4)
        lowerTerm = 1. / 12. * (term1 - term2 - 6. * u1 * rho_u**2 + 6. * rho_u**2)
        FstarUpper = 2. * np.pi * (upperTerm - lowerTerm) * np.ones_like(wavelengthArray)
    

    else: # RM_effect = True

        phiArray, delta_phi = constructAxis(gridsDict, 'phi')
        rhoArray, delta_rho = constructAxis(gridsDict, 'rho')

        phiGrid, rhoGrid = np.meshgrid(phiArray, rhoArray, indexing = 'ij')

        GRID = np.stack((phiGrid.flatten(), rhoGrid.flatten()), axis = -1)

        FstarIntegrated = np.zeros_like(wavelengthArray)

        for point in GRID:

            phi = point[0]
            rho = point[1]
        
            Fstar = calculateRM(phi, rho, wavelengthArray, architectureDict, PHOENIX_output)

            if fundamentalsDict['CLV_variations']:

                Fstar *= calculateCLV(rho, R_star, architectureDict['u1'], architectureDict['u2'])
            
            FstarIntegrated += Fstar * delta_phi * delta_rho * rho

            FstarUpper = np.zeros_like(wavelengthArray) # To calculate this with the RM-effect we would need the spatial grid to run over the entire stellar disk

    return FstarIntegrated, FstarUpper

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

        n = gasprop.getNumberDensity(phi, rho, orbphase, xArray, key_scenario, specificScenarioDict, architectureDict)     

        sigma_abs = gasprop.getAbsorptionCrossSection(phi, rho, orbphase, xArray, wavelengthArray, key_scenario, fundamentalsDict, specificScenarioDict, architectureDict, speciesDict, sigmaLookupDict)

        tau += delta_x * np.sum(np.multiply(sigma_abs, n), axis = 1)

    return tau # tau(lambda)


def evaluateChord(point, args): # Function to be multiprocessed

    phi, rho, orbphase = point

    delta_x, delta_phi, delta_rho, xArray, wavelengthArray, architectureDict, fundamentalsDict, scenarioDict, speciesDict, gridsDict, outputDict, sigmaLookupDict, Fstar_function = args

    if checkBlock(phi, rho, orbphase, architectureDict, fundamentalsDict):

        tau = np.inf * np.ones_like(wavelengthArray)

    else:

        tau = calculateOpticalDepth(phi, rho, orbphase, xArray, wavelengthArray, fundamentalsDict, architectureDict, scenarioDict, speciesDict, gridsDict, sigmaLookupDict)

    Fstar = getFstar(phi, rho, wavelengthArray, architectureDict, fundamentalsDict, Fstar_function)

    singleChord = rho * Fstar * np.exp(-tau) * delta_phi * delta_rho

    if outputDict['recordTau']:

        return singleChord, tau

    return singleChord, None

def prepareArguments(fundamentalsDict, architectureDict, scenarioDict, speciesDict, gridsDict, outputDict, startTime, verbose):

    xArray, delta_x = constructAxis(gridsDict, 'x')

    wavelengthArray = constructWavelengthGrid(gridsDict, scenarioDict, speciesDict)

    phiArray, delta_phi = constructAxis(gridsDict, 'phi')
    rhoArray, delta_rho = constructAxis(gridsDict, 'rho')
    orbphaseArray = constructAxis(gridsDict, 'orbphase')

    phiGrid, rhoGrid, orbphaseGrid = np.meshgrid(phiArray, rhoArray, orbphaseArray, indexing = 'ij')

    GRID = np.stack((phiGrid.flatten(), rhoGrid.flatten(), orbphaseGrid.flatten()), axis = -1)

    Fstar_function = getPHOENIX_output(wavelengthArray, fundamentalsDict, architectureDict)

    if verbose:
        print('\nPHOENIX spectrum downloaded (if required):', datetime.now() - startTime)

    FstarIntegrated, FstarUpper = getFstarIntegrated(wavelengthArray, fundamentalsDict, architectureDict, gridsDict, Fstar_function) # Stellar flux as a function of wavelength, integrated over the stellar disk

    if verbose:
        print('Integrated stellar flux calculated:', datetime.now() - startTime)

    sigmaLookupDict = gasprop.createLookupAbsorption(xArray, wavelengthArray, GRID, fundamentalsDict, architectureDict, scenarioDict, speciesDict)

    if verbose:
        print('Lookup absorption cross section calculated:', datetime.now() - startTime)

    args = (delta_x, delta_phi, delta_rho, xArray, wavelengthArray, architectureDict, fundamentalsDict, scenarioDict, speciesDict, gridsDict, outputDict, sigmaLookupDict, Fstar_function)

    return GRID, args, FstarIntegrated, FstarUpper


def sumOverChords(RESULTS_R, gridsDict, FstarIntegrated, FstarUpper):

    N_phi = int(gridsDict['phi_steps'])
    N_rho = int(gridsDict['rho_steps'])
    N_orbphase = int(gridsDict['orbphase_steps'])
    N_lambda = len(RESULTS_R[0])

    ChordArray = np.array(RESULTS_R).reshape((N_phi * N_rho, N_orbphase, N_lambda))
    ChordSum = np.sum(ChordArray, axis = 0) # Sum over rho and phi coordinates

    R = ((ChordSum + np.tile(FstarUpper, (N_orbphase, 1))) / np.tile(FstarIntegrated, (N_orbphase, 1))).T

    index_max = np.unravel_index(np.argmin(R, axis = None), R.shape)

    return R.flatten(), index_max

def getTau(RESULTS_TAU, architectureDict, gridsDict, index_max):

    N_phi = int(gridsDict['phi_steps'])
    N_rho = int(gridsDict['rho_steps'])
    N_orbphase = int(gridsDict['orbphase_steps'])
    N_lambda = len(RESULTS_TAU[0])

    rho_axis = constructAxis(gridsDict, 'rho')[0]

    phi_axis = constructAxis(gridsDict, 'phi')[0]
    
    phi, rho = np.meshgrid(phi_axis, rho_axis, indexing = 'ij')
    phi = phi.flatten()
    rho = rho.flatten()

    TAUarray = np.array(RESULTS_TAU).reshape((N_phi, N_rho, N_orbphase, N_lambda))
    tauDisk = TAUarray[:, :, index_max[1], index_max[0]].flatten()

    return phi, rho, tauDisk