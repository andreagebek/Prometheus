"""
Functions to calculate the positions/velocities of exoplanet/exomoon.
Created on 18. October 2021 by Andrea Gebek.
"""

import numpy as np
import eliteScripts.constants as const

def getPlanetPosition(architectureDict, orbphase):

    x_p_fromStar = architectureDict['a_p'] * np.cos(orbphase)    # The actual x_p is always 0
    y_p = architectureDict['a_p'] * np.sin(orbphase)

    return x_p_fromStar, y_p

def getOrbphaseMoon(architectureDict, orbphase):

    orbphase_moon = architectureDict['starting_orbphase_moon'] + orbphase * np.sqrt((architectureDict['a_p']**3 * architectureDict['M_p']) / (architectureDict['a_moon']**3 * architectureDict['M_star']))

    return orbphase_moon

def getMoonPosition(architectureDict, orbphase):

    orbphase_moon = getOrbphaseMoon(architectureDict, orbphase)

    x_moon = architectureDict['a_moon'] * np.cos(orbphase_moon)     # Planet is always at x = 0
    y_moon = getPlanetPosition(architectureDict, orbphase)[1] + architectureDict['a_moon'] * np.sin(orbphase_moon)

    return x_moon, y_moon


def getBodyLOSvelocity(architectureDict, orbphase, key_scenario):

    v_los = -np.sin(orbphase) * np.sqrt(const.G * architectureDict['M_star'] / architectureDict['a_p'])

    if key_scenario == 'exomoon':

        orbphase_moon = getOrbphaseMoon(architectureDict, orbphase)
        v_los += - np.sin(orbphase_moon) * np.sqrt(const.G * architectureDict['M_p'] / architectureDict['a_moon'])

    return v_los

def getPlanetRotationLOSvelocity(architectureDict, phi, rho, orbphase, key_scenario):

    if key_scenario == 'barometric' or key_scenario == 'hydrostatic' or key_scenario == 'escaping':

        R_0 = architectureDict['R_0']
        period_planetrot = architectureDict['period_planetrot']

        y_p = getPlanetPosition(architectureDict, orbphase)[1]

        v_los = 2 * np.pi / period_planetrot * (rho * np.sin(phi) - y_p)

    else:

        v_los = np.zeros(np.shape(phi))
    
    return v_los

def getCartesianFromCylinder(x, phi, rho):

    y = rho * np.sin(phi)
    z = rho * np.cos(phi)

    return x, y, z


def getDistanceFromPlanet(architectureDict, x, phi, rho, orbphase):

    x, y, z = getCartesianFromCylinder(x, phi, rho)
    
    y_p = getPlanetPosition(architectureDict, orbphase)[1]

    r_fromPlanet = np.sqrt(x**2 + (y - y_p)**2 + z**2)

    return r_fromPlanet

def getTorusCoords(architectureDict, x, phi, rho, orbphase):

    x, y, z = getCartesianFromCylinder(x, phi, rho)
    
    y_p = getPlanetPosition(architectureDict, orbphase)[1]

    a = np.sqrt(x**2 + (y - y_p)**2)

    return a, z


def getDistanceFromMoon(architectureDict, x, phi, rho, orbphase):

    x, y, z = getCartesianFromCylinder(x, phi, rho)

    x_moon, y_moon = getMoonPosition(architectureDict, orbphase)
    
    r_fromMoon = np.sqrt((x - x_moon)**2 + (y - y_moon)**2 + z**2)

    return r_fromMoon

