"""
Calculate the light curve of a transiting exoplanet.
Created on 15. July 2021 by Andrea Gebek.
"""

import numpy as np
from scipy.interpolate import interp1d
import sys
import os
import datetime
SCRIPTPATH = os.path.realpath(__file__)
GITPATH = os.path.dirname(os.path.dirname(SCRIPTPATH))
sys.path.append(GITPATH) 
import prometheusScripts.constants as const
import prometheusScripts.geometryHandler as geom
import prometheusScripts.gasProperties as gasprop
import prometheusScripts.stellarSpectrum as stellar

def constructAxis(gridsDict, architectureDict, axisName):

    if axisName == 'x':

        x_border = gridsDict['x_border']
        x_steps = gridsDict['x_steps']
        x_axis = np.linspace(-x_border, x_border, int(x_steps) + 1)[:-1] + x_border / float(x_steps)  

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

def constructSpatialGrid(gridsDict, architectureDict):

    x_axis = constructAxis(gridsDict, architectureDict, 'x')

    rho_axis = constructAxis(gridsDict, architectureDict, 'rho')

    phi_axis = constructAxis(gridsDict, architectureDict, 'phi')

    orbphase_axis = constructAxis(gridsDict, architectureDict, 'orbphase')

    x, phi, rho, orbphase = np.meshgrid(x_axis, phi_axis, rho_axis, orbphase_axis, indexing = 'ij')
    
    return x, phi, rho, orbphase

def calculateCLV(rho, R_star, u1, u2):

    arg = 1. - np.sqrt(1. - rho**2 / R_star**2)

    return 1. - u1 * arg - u2 * arg**2

def calculateRM(wavelength, architectureDict, gridsDict):

    T_starrot = architectureDict['period_starrot']
    R_star = architectureDict['R_star']

    phi_axis = constructAxis(gridsDict, architectureDict, 'phi')
    rho_axis = constructAxis(gridsDict, architectureDict, 'rho')

    PHOENIX_output = stellar.readSpectrum(architectureDict['T_eff'], architectureDict['log_g'], architectureDict['Fe_H'], architectureDict['alpha_Fe'])
    w_star = PHOENIX_output[0]
    v_max = 2. * np.pi * R_star / T_starrot
    w_max = np.max(wavelength * gasprop.calculateDopplerShift(-v_max))
    w_min = np.min(wavelength * gasprop.calculateDopplerShift(v_max))
    SEL_w = np.argwhere((w_star > w_min) * (w_star < w_max))[:, 0]

    SEL = np.concatenate((np.array([np.min(SEL_w) - 1]), SEL_w, np.array([np.max(SEL_w) + 1])))

    w_star = w_star[SEL]
    F_0 = PHOENIX_output[1][SEL]

    v_los = geom.calculateStellarLOSvelocity(architectureDict, phi_axis, rho_axis)
    w_shift = gasprop.calculateDopplerShift(v_los)

    F_function = interp1d(w_star, F_0, kind = 'cubic')
    F_shifted = F_function(np.tensordot(wavelength, 1. / w_shift, axes = 0)) # This contains the shifted flux incident on the exoplanet, instead of shifting the spectra
    # at all positions on the stellar surface just read in the spectrum at different wavelengths depending on the position in the spatial grid (phi and rho)

    return F_shifted

def getStarFactors(wavelength, fundamentalsDict, architectureDict, gridsDict):
    # Return F_star(wavelength, phi, rho) and its integral over the stellar disk

    R_star = architectureDict['R_star']
    phi_steps = gridsDict['phi_steps']
    rho_steps = gridsDict['rho_steps']

    rho_axis = constructAxis(gridsDict, architectureDict, 'rho')

    if not fundamentalsDict['CLV_variations'] and not fundamentalsDict['RM_effect']:

        F_star = np.ones((len(wavelength), int(phi_steps), int(rho_steps)))
        F_star_integrated = np.pi * R_star**2

        return F_star, F_star_integrated

    elif fundamentalsDict['CLV_variations'] and not fundamentalsDict['RM_effect']:

        F_star = np.tile(calculateCLV(rho_axis, R_star, architectureDict['u1'], architectureDict['u2']), (len(wavelength), int(phi_steps), 1))
        F_star_integrated = np.pi * R_star**2 * (1. - architectureDict['u1'] / 3. - architectureDict['u2'] / 6.)

        return F_star, F_star_integrated

    else: # RM_effect = True

        F_star = calculateRM(wavelength, architectureDict, gridsDict)

        if fundamentalsDict['CLV_variations']:

            F_star *= np.tile(calculateCLV(rho_axis, R_star, architectureDict['u1'], architectureDict['u2']), (len(wavelength), int(phi_steps), 1))
        
        delta_rho = R_star / float(rho_steps)
        delta_phi = 2 * np.pi / float(phi_steps)

        F_star_integralphi = delta_phi * np.sum(F_star, axis = 1)
        F_star_integralrho = delta_rho * np.tensordot(rho_axis, F_star_integralphi, axes = [0, 1])
        F_star_integrated = np.tile(F_star_integralrho, (int(gridsDict['orbphase_steps']), 1)).T    
        # Here, the integrated F_star depends on wavelength. Promote it to F_star_integrated(wavelength, orbphase).

        return F_star, F_star_integrated




def calculateOpticalDepth(x, phi, rho, orbphase, wavelength, fundamentalsDict, architectureDict, scenarioDict, speciesDict, gridsDict, outputDict, startTime):

    delta_x = 2 * gridsDict['x_border'] / float(gridsDict['x_steps'])

    tau = 0

    for key_scenario in scenarioDict.keys():
        print('\nOptical depth calculation for the ' + key_scenario + ' scenario starting:', datetime.datetime.now() - startTime)
        specificScenarioDict = scenarioDict[key_scenario]

        n = gasprop.getNumberDensity(x, phi, rho, orbphase, key_scenario, specificScenarioDict, architectureDict, fundamentalsDict)     
        print('Number density calculation for the ' + key_scenario + ' scenario finished:', datetime.datetime.now() - startTime)
        sigma_abs = gasprop.getAbsorptionCrossSection(x, phi, rho, orbphase, wavelength, key_scenario, fundamentalsDict, specificScenarioDict, architectureDict, speciesDict, startTime)
        print('Total absorption cross section calculation for the ' + key_scenario + ' scenario finished:', datetime.datetime.now() - startTime)
        n = np.tile(n, (len(wavelength), 1, 1, 1, 1))
        
        BLOCK = (n == np.inf)

        sigma_abs[BLOCK] = 1.   # If sigma_abs = 0 and n = inf we have an issue for the calculation of tau

        tau += delta_x * np.sum(np.multiply(sigma_abs, n), axis = 1)
        print('Optical depth calculation for the ' + key_scenario + ' scenario finished:', datetime.datetime.now() - startTime, '\n')
    return tau # tau(lambda, phi, rho, t)


def calculateTransitDepth(fundamentalsDict, architectureDict, scenarioDict, speciesDict, gridsDict, outputDict, startTime):
    """Calculate the wavelength-dependent transit depth
    """

    x, phi, rho, orbphase = constructSpatialGrid(gridsDict, architectureDict)

    wavelength = constructAxis(gridsDict, architectureDict, 'wavelength')

    tau = calculateOpticalDepth(x, phi, rho, orbphase, wavelength, fundamentalsDict, architectureDict, scenarioDict, speciesDict, gridsDict, outputDict, startTime)
    print('Total optical depth calculated:', datetime.datetime.now() - startTime)
    singleChord = np.exp(-tau)

    delta_rho = architectureDict['R_star'] / float(gridsDict['rho_steps'])
    delta_phi = 2 * np.pi / float(gridsDict['phi_steps'])

    F_star_untiled, F_star_integrated = getStarFactors(wavelength, fundamentalsDict, architectureDict, gridsDict)
    print('Calculations with the stellar flux finished:', datetime.datetime.now() - startTime)
    F_star_tiled = np.tile(F_star_untiled, (int(gridsDict['orbphase_steps']), 1, 1, 1))
    F_star = np.moveaxis(F_star_tiled, 0, -1) # F_star(lambda, phi, rho, orbphase)

    integral_phi = delta_phi * np.sum(np.multiply(F_star, singleChord), axis = 1) # f(lambda, rho, orbphase)

    sum_over_chords = delta_rho * np.tensordot(rho[0, 0, :, 0], integral_phi, axes = [0, 1]) # Integrate along the rho-axis

    R = sum_over_chords / F_star_integrated # R(lambda, orbphase)

    resultsDict = {'R': R}

    if outputDict['recordTau']:

        argmaxR =  np.unravel_index(np.argmin(R, axis = None), sum_over_chords.shape)

        tauDisk = tau[argmaxR[0], :, :, argmaxR[1]]
        
        resultsDict['tauDisk'] = tauDisk

    if outputDict['benchmark']:

        R_benchmark = calculateBarometricBenchmark(x, phi, rho, orbphase, wavelength, fundamentalsDict, architectureDict, scenarioDict['barometric'], speciesDict)
        print('Benchmark transit spectrum calculated:', datetime.datetime.now() - startTime)

        resultsDict['R_benchmark'] = R_benchmark

    return resultsDict


def calculateBarometricBenchmark(x, phi, rho, orbphase, wavelength, fundamentalsDict, architectureDict, specificScenarioDict, speciesDict):
    T = specificScenarioDict['T']
    P_0 = specificScenarioDict['P_0']
    mu = specificScenarioDict['mu']
    R_0 = architectureDict['R_0']
    M_p = architectureDict['M_p']
    R_star = architectureDict['R_star']

    H = const.k_B * T * R_0**2 / (const.G * mu * M_p)
    g = const.G * M_p / R_0**2

    sigma_abs = gasprop.getAbsorptionCrossSection(x, phi, rho, orbphase, wavelength, 'barometric', fundamentalsDict, specificScenarioDict, architectureDict, speciesDict)[:, 0, 0, 0, :]
    # Absorption cross section for the benchmark cannot depend on the spatial position, i.e. sigma_abs(lambda, t)

    kappa = sigma_abs / mu
    R_lambda = R_0 + H * (const.euler_mascheroni + np.log(P_0 * kappa / g * np.sqrt(2 * np.pi * R_0 / H)))
    benchmark_spectrum = (R_star**2 - R_lambda**2) / R_star**2

    return benchmark_spectrum

