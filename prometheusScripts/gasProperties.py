"""
This file stores various functions related
to the properties of the gas (e.g. number densities, velocities,
absorption cross sections).
Created on 19. October 2021 by Andrea Gebek.
"""

import numpy as np
import sys
from scipy.special import erf, voigt_profile
from scipy.interpolate import interp1d, RegularGridInterpolator
import os
import h5py
SCRIPTPATH = os.path.realpath(__file__)
GITPATH = os.path.dirname(os.path.dirname(SCRIPTPATH))
sys.path.append(GITPATH)
import prometheusScripts.constants as const
import prometheusScripts.geometryHandler as geom
import datetime

startTime = datetime.datetime.now()

"""
Number density functions
"""

def getNumberDensity(phi, rho, orbphase, xArray, key_scenario, specificScenarioDict, architectureDict, fundamentalsDict):

    if key_scenario == 'barometric':

        r_fromP = geom.getDistanceFromPlanet(architectureDict, phi, rho, orbphase, xArray)
        n = getBarometric_n(r_fromP, specificScenarioDict, architectureDict)

    elif key_scenario == 'hydrostatic':

        r_fromP = geom.getDistanceFromPlanet(architectureDict, phi, rho, orbphase, xArray)
        n = getHydrostatic_n(r_fromP, specificScenarioDict, architectureDict)

    elif key_scenario == 'escaping':

        r_fromP = geom.getDistanceFromPlanet(architectureDict, phi, rho, orbphase, xArray)
        n = getEscaping_n(r_fromP, specificScenarioDict, architectureDict)

    elif key_scenario == 'exomoon':

        r_fromMoon = geom.getDistanceFromMoon(architectureDict, phi, rho, orbphase, xArray)
        n = getExomoon_n(r_fromMoon, specificScenarioDict, architectureDict)

    elif key_scenario == 'torus':

        a, z = geom.getTorusCoords(architectureDict, phi, rho, orbphase, xArray)
        n = getTorus_n(a, z, specificScenarioDict, architectureDict)

    else:

        print('This key for a number density structure does not exist. The program exits now.')
        sys.exit(1)
    
    return n


def getBarometric_n(r, specificScenarioDict, architectureDict):
    """Barometric law to benchmark the analytical solution
    """
    T = specificScenarioDict['T']
    P_0 = specificScenarioDict['P_0']
    mu = specificScenarioDict['mu']
    M_p = architectureDict['M_p']
    R_0 = architectureDict['R_0']

    n_0 = P_0 / (const.k_B * T)
    H = const.k_B * T * R_0**2 / (const.G * mu * M_p)
    n = n_0 * np.exp((R_0 - r) / H)

    return n 


def getHydrostatic_n(r, specificScenarioDict, architectureDict):
    """Hydrostatic equilibrium with constant temperature and
    mean molecular weight
    """
    T = specificScenarioDict['T']
    P_0 = specificScenarioDict['P_0']
    mu = specificScenarioDict['mu']
    M_p = architectureDict['M_p']
    R_0 = architectureDict['R_0']

    n_0 = P_0 / (const.k_B * T)
    Jeans_0 = const.G * mu * M_p / (const.k_B * T * R_0)
    Jeans = const.G * mu * M_p / (const.k_B * T * r) * np.heaviside(r - R_0, 1) 
    # When the number density is evaluated inside the planet the Jeans parameter gets very large, leading to overflowing n
    n = n_0 * np.exp(Jeans - Jeans_0) 

    return n    



def getEscaping_n(r, specificScenarioDict, architectureDict):
    """Power law with index q_esc, for a hydrodynamically
    escaping atmosphere
    """
    q_esc = specificScenarioDict['q_esc']
    R_0 = architectureDict['R_0']

    if 'T' in specificScenarioDict.keys():

        T = specificScenarioDict['T']
        P_0 = specificScenarioDict['P_0']
        n_0_esc = P_0 / (const.k_B * T)
    
    else:

        n_0_esc = (q_esc - 3.) / (4. * np.pi * R_0**3)

    n = n_0_esc * (R_0 / r)**q_esc

    return n


def getExomoon_n(r, specificScenarioDict, architectureDict):
    """Scaled to observed number density profile at Io (Burger 2001)
    """
    q_moon = specificScenarioDict['q_moon']
    R_moon = architectureDict['R_moon']

    n_0_exomoon = (q_moon - 3.) / (4 * np.pi * R_moon**3)
    n = n_0_exomoon * (R_moon / r)**(q_moon)

    return n

def getTorus_n(a, z, specificScenarioDict, architectureDict):
    """ Toroidal atmosphere (Johnson & Huggins 2006)
    """
    a_torus = specificScenarioDict['a_torus']
    v_ej = specificScenarioDict['v_ej']
    M_p = architectureDict['M_p']

    v_orbit = np.sqrt(const.G * M_p / a_torus)
    H_torus = a_torus * v_ej / v_orbit

    n_a = np.exp(- ((a - a_torus) / (4 * H_torus))**2)
    n_z = np.exp(- (z / H_torus)**2)

    term1 = 8. * H_torus**2 * np.exp(-a_torus**2 / (16 * H_torus**2))
    term2 = 2. * np.sqrt(np.pi) * a_torus * H_torus * (erf(a_torus / (4 * H_torus)) + 1.)

    n_0_torus = 1. / (2. * np.pi**1.5 * H_torus * (term1 + term2))

    n = n_0_torus * np.multiply(n_a, n_z) 
    
    return n


"""
Velocity functions (for Doppler shifts)
"""


def getLOSVelocity(phi, rho, orbphase, xArray, fundamentalsDict, key_scenario, specificScenarioDict, architectureDict):

    v_los = np.zeros(np.shape(xArray))

    if fundamentalsDict['DopplerOrbitalMotion']:

        v_los += geom.getBodyLOSvelocity(architectureDict, orbphase, key_scenario)

    if fundamentalsDict['DopplerPlanetRotation']:

        v_los += geom.getPlanetRotationLOSvelocity(architectureDict, phi, rho, orbphase, key_scenario)

    if key_scenario == 'escaping' and specificScenarioDict['RadialWind']:

        v_los += getEscaping_v(phi, rho, orbphase, xArray, specificScenarioDict, architectureDict)

    elif key_scenario == 'exomoon' and specificScenarioDict['RadialWind']:

        v_los += getExomoon_v(phi, rho, orbphase, xArray, specificScenarioDict, architectureDict)

    return v_los


def getEscaping_v(phi, rho, orbphase, xArray, specificScenarioDict, architectureDict):

    v_0 = specificScenarioDict['vRadial_0']
    q_esc = specificScenarioDict['q_esc']
    R_0 = architectureDict['R_0']

    r = geom.getDistanceFromPlanet(architectureDict, phi, rho, orbphase, xArray)

    v_los = v_0 * (r / R_0) ** (q_esc - 2.) * xArray / r

    return v_los

def getExomoon_v(phi, rho, orbphase, xArray, specificScenarioDict, architectureDict):

    v_0 = specificScenarioDict['vRadial_0']
    q_moon = specificScenarioDict['q_moon']
    R_moon = architectureDict['R_moon']

    r = geom.getDistanceFromMoon(architectureDict, phi, rho, orbphase, xArray)
    x_moon = geom.getMoonPosition(architectureDict, orbphase)[0]

    v_los = v_0 * (r / R_moon)**(q_moon - 2.) * (xArray - x_moon) / r

    return v_los


"""
Calculate absorption cross sections
"""

def readLineList(key_species, wavelength):

    LineList = np.loadtxt(GITPATH + '/Resources/LineList.txt', dtype = str, usecols = (0, 1, 2, 3, 4), skiprows = 1)

    line_wavelength = np.array([x[1:-1] for x in LineList[:, 2]])
    line_A = np.array([x[1:-1] for x in LineList[:, 3]])
    line_f = np.array([x[1:-1] for x in LineList[:, 4]])
    
    SEL_COMPLETE = (line_wavelength != '') * (line_A != '') * (line_f != '') 

    SEL_SPECIES = (LineList[:, 0] == const.speciesInfoDict[key_species][0]) * (LineList[:, 1] == const.speciesInfoDict[key_species][1]) # Select element and ionization state

    line_wavelength = line_wavelength[SEL_SPECIES * SEL_COMPLETE].astype(np.float) * 1e-8 # In cm
    line_gamma = line_A[SEL_SPECIES * SEL_COMPLETE].astype(np.float) / (4. * np.pi) # HWHM of Lorentzian, assumed to be A / (4 * pi) (see Draine 2013 p57)
    line_f = line_f[SEL_SPECIES * SEL_COMPLETE].astype(np.float)

    SEL_WAVELENGTH = (line_wavelength > 0.95 * min(wavelength)) * (line_wavelength < 1.05 * max(wavelength))

    return line_wavelength[SEL_WAVELENGTH], line_gamma[SEL_WAVELENGTH], line_f[SEL_WAVELENGTH]

def calculateLineAbsorption(wavelength, line_wavelength, line_gamma, line_f, specificSpeciesDict):

    chi = specificSpeciesDict['chi']
    sigma_v = specificSpeciesDict['sigma_v']

    sigma_abs_species = np.zeros_like(wavelength)

    for idx in range(len(line_wavelength)):

        lineProfile = voigt_profile(const.c / wavelength - const.c / line_wavelength[idx], sigma_v / line_wavelength[idx], line_gamma[idx]) # Calculation in frequency space
        sigma_abs_species += np.pi * (const.e)**2 / (const.m_e * const.c) * line_f[idx] * chi * lineProfile # Draine 2013

    return sigma_abs_species


def calculateDopplerShift(v_los):
    # If v_los is positive, receiver and source are moving towards each other

    if np.any(np.abs(v_los) >= const.c):
        print('\nCRITICAL ERROR: Velocities of the gas in the simulation exceed speed of light. PROMETHEUS exits now.')
        sys.exit()

    beta = v_los / const.c
    shift = np.sqrt((1. - beta) / (1. + beta))

    return shift

def createLookupAbsorption(xArray, wavelengthArray, GRID, fundamentalsDict, architectureDict, scenarioDict, speciesDict):

    if fundamentalsDict['ExactSigmaAbs']:

        return None

    else:

        sigmaLookupDict = {}

        for key_scenario in scenarioDict.keys():

            v_los = []

            for point in GRID:

                phi = point[0]
                rho = point[1]
                orbphase = point[2]
            
                v_los.append(getLOSVelocity(phi, rho, orbphase, xArray, fundamentalsDict, key_scenario, scenarioDict[key_scenario], architectureDict))

            v_los_max = np.max(np.abs(v_los))

            w_min = np.min(wavelengthArray) * calculateDopplerShift(v_los_max)
            w_max = np.max(wavelengthArray) * calculateDopplerShift(-v_los_max)
            wavelengthHighRes = np.arange(w_min, w_max, fundamentalsDict['LookupResolution'])

            sigmaHighRes = np.zeros_like(wavelengthHighRes)

            for key_species in speciesDict[key_scenario].keys():

                if 'sigma_v' in speciesDict[key_scenario][key_species].keys():

                    line_wavelength, line_gamma, line_f = readLineList(key_species, wavelengthArray)

                    sigmaHighRes += calculateLineAbsorption(wavelengthHighRes, line_wavelength, line_gamma, line_f, speciesDict[key_scenario][key_species])


            if 'RayleighScatt' in scenarioDict[key_scenario].keys():

                if scenarioDict[key_scenario]['RayleighScatt']:

                    sigmaHighRes += 8.49e-45 / wavelengthHighRes**4 # FROM WHERE IS THIS? DEPENDENCY ON H2 mixing ratio?
            
            sigmaLookupFunction = interp1d(wavelengthHighRes, sigmaHighRes, kind = 'cubic')
            sigmaLookupDict[key_scenario] = (sigmaLookupFunction, wavelengthHighRes)

    return sigmaLookupDict
 

def getAbsorptionCrossSection(phi, rho, orbphase, xArray, wavelengthArray, key_scenario, fundamentalsDict, specificScenarioDict, architectureDict, speciesDict, sigmaLookupDict):
    # Note that this absorption cross section is already multiplied by either the mixing ratio or the total number of absorbing atoms,
    # hence it doesn't reflect the bare absorption cross section per particle

    v_los = getLOSVelocity(phi, rho, orbphase, xArray, fundamentalsDict, key_scenario, specificScenarioDict, architectureDict)

    wavelengthShiftFactor = calculateDopplerShift(-v_los)
    wavelengthShifted = np.tensordot(wavelengthArray, wavelengthShiftFactor, axes = 0)

    if fundamentalsDict['ExactSigmaAbs']:

        sigma_abs = 0

        for key_species in speciesDict[key_scenario].keys():

            line_wavelength, line_gamma, line_f = readLineList(key_species, wavelengthArray)

            sigma_abs += calculateLineAbsorption(wavelengthShifted, line_wavelength, line_gamma, line_f, speciesDict[key_scenario][key_species])            
            
        if 'RayleighScatt' in specificScenarioDict.keys():

            if specificScenarioDict['RayleighScatt']:

                sigma_abs += 8.49e-45 / wavelengthShifted**4 # FROM WHERE IS THIS? DEPENDENCY ON H2 mixing ratio?

    else:

        sigmaLookupFunction, wavelengthHighRes = sigmaLookupDict[key_scenario]

        sigmaLookupUnclipped = sigmaLookupFunction(np.clip(wavelengthShifted, np.min(wavelengthHighRes), np.max(wavelengthHighRes))) # Evaluate the lookup absorption cross 
        #section at Doppler-shifted wavelengths. Clip because of rounding & fitting errors

        if np.any(sigmaLookupUnclipped < 0.):
            print('\nWARNING: The absorption cross section is smaller than zero somewhere. This is probably due to an insufficient lookup resolution. \
The absorption cross section will be clipped to positive values.\n')

        sigma_abs = np.clip(sigmaLookupUnclipped, 0., a_max = None)
    
    return sigma_abs