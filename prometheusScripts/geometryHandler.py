"""
Functions to calculate the positions/velocities of exoplanet/exomoon.
Created on 18. October 2021 by Andrea Gebek.
"""

import numpy as np
import sys
import os
SCRIPTPATH = os.path.realpath(__file__)
GITPATH = os.path.dirname(os.path.dirname(SCRIPTPATH))
sys.path.append(GITPATH)
import prometheusScripts.constants as const

def getPlanetPosition(architectureDict, orbphase): # Circular orbit

    x_p = architectureDict['a_p'] * np.cos(orbphase)
    y_p = architectureDict['a_p'] * np.sin(orbphase)

    return x_p, y_p

def getOrbphaseMoon(architectureDict, orbphase):

    orbphase_moon = architectureDict['starting_orbphase_moon'] + orbphase * np.sqrt((architectureDict['a_p']**3 * architectureDict['M_p']) / (architectureDict['a_moon']**3 * architectureDict['M_star']))

    return orbphase_moon

def getMoonPosition(architectureDict, orbphase):

    orbphase_moon = getOrbphaseMoon(architectureDict, orbphase)

    x_p, y_p = getPlanetPosition(architectureDict, orbphase)

    x_moon = x_p + architectureDict['a_moon'] * np.cos(orbphase_moon)
    y_moon = y_p + architectureDict['a_moon'] * np.sin(orbphase_moon)

    return x_moon, y_moon


def getBodyLOSvelocity(architectureDict, orbphase, key_scenario):

    v_los = -np.sin(orbphase) * np.sqrt(const.G * architectureDict['M_star'] / architectureDict['a_p'])

    if key_scenario == 'exomoon':

        orbphase_moon = getOrbphaseMoon(architectureDict, orbphase)
        v_los += - np.sin(orbphase_moon) * np.sqrt(const.G * architectureDict['M_p'] / architectureDict['a_moon'])

    return v_los

def getPlanetRotationLOSvelocity(architectureDict, phi, rho, orbphase, key_scenario):

    if key_scenario == 'barometric' or key_scenario == 'hydrostatic' or key_scenario == 'escaping':

        period_planetrot = architectureDict['period_planetrot']

        y_p = getPlanetPosition(architectureDict, orbphase)[1]

        v_los = 2 * np.pi / period_planetrot * (rho * np.sin(phi) - y_p)

    else:

        v_los = 0.
    
    return v_los

def getCartesianFromCylinder(phi, rho):

    y = rho * np.sin(phi)
    z = rho * np.cos(phi)

    return y, z


def getDistanceFromPlanet(architectureDict, phi, rho, orbphase, xArray):

    y, z = getCartesianFromCylinder(phi, rho)
    
    x_p, y_p = getPlanetPosition(architectureDict, orbphase)

    r_fromPlanet = np.sqrt((xArray - x_p)**2 + (y - y_p)**2 + z**2)

    return r_fromPlanet

def getTorusCoords(architectureDict, phi, rho, orbphase, xArray):

    y, z = getCartesianFromCylinder(phi, rho)
    
    x_p, y_p = getPlanetPosition(architectureDict, orbphase)

    a = np.sqrt((xArray - x_p)**2 + (y - y_p)**2)

    return a, z


def getDistanceFromMoon(architectureDict, phi, rho, orbphase, xArray):

    y, z = getCartesianFromCylinder(phi, rho)

    x_moon, y_moon = getMoonPosition(architectureDict, orbphase)
    
    r_fromMoon = np.sqrt((xArray - x_moon)**2 + (y - y_moon)**2 + z**2)

    return r_fromMoon


def calculateStellarLOSvelocity(architectureDict, phi, rho):
    
    i_starrot = architectureDict['inclination_starrot']
    phi_starrot = architectureDict['azimuth_starrot']
    T_starrot = architectureDict['period_starrot']
    R_star = architectureDict['R_star']

    dir_omega = np.array([-np.sin(i_starrot) * np.cos(phi_starrot), -np.sin(i_starrot) * np.sin(phi_starrot), np.cos(i_starrot)])
    omega = 2. * np.pi / T_starrot * dir_omega # Angular velocity vector of the stellar rotation
    r_surface = np.array([np.tensordot(np.ones(len(phi)), np.sqrt(R_star**2 - rho**2), axes = 0), 
    np.tensordot(np.sin(phi), rho, axes = 0), np.tensordot(np.cos(phi), rho, axes = 0)])
    # Vector to the surface of the star
    v_los = np.cross(omega, r_surface, axisb = 0)[:, :, 0] # The line-of-sight velocity is the one along the x-axis

    return v_los