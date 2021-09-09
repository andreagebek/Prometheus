"""
This file stores functions to calculate 
velocities of the absorbers, for the simulation
of a lightcurve.
Created on 1. September 2021 by Andrea Gebek.
"""

from constants import *
import numpy as np


def v_escaping(ArchitectureDict, EscapingDict, gridgrid):
    # Radially escaping wind in the escaping scenario
    
    v_0 = EscapingDict['vRadial_0']
    q_esc = EscapingDict['q_esc']
    R_0 = ArchitectureDict['R_0']



    xx, phiphi, rhorho, yy_pp = gridgrid

    r_fromP = np.sqrt(xx**2 + (rhorho * np.sin(phiphi) - yy_pp)**2 + (rhorho * np.cos(phiphi))**2)

    v = v_0 * (r_fromP / R_0) ** (q_esc - 2.) * xx / r_fromP # Constant outwards mass flux

    return v

def v_exomoon(orbphase, ArchitectureDict, ExomoonDict, gridgrid):
    # Radially escaping wind in the exomoon scenario
    v_0 = ExomoonDict['vRadial_0']

    R_0 = ArchitectureDict['R_moon']
    a_moon = ArchitectureDict['a_moon']
    starting_orbphase_moon = ArchitectureDict['starting_orbphase_moon']
    a_p = ArchitectureDict['a_p']
    M_p = ArchitectureDict['M_p']
    M_s = ArchitectureDict['M_star']

    orbphase_moon = starting_orbphase_moon + orbphase * np.sqrt((a_p**3 * M_p) / (a_moon**3 * M_s)) # Assuming Keplerian orbits of massless points

    xx, phiphi, rhorho, yy_pp = gridgrid

    x_moonFrame = xx - a_moon * np.cos(orbphase_moon)
    y_moonFrame = np.sin(phiphi) * rhorho - yy_pp - a_moon * np.sin(orbphase_moon)
    z_moonFrame = np.cos(phiphi) * rhorho

    r_fromMoon = np.sqrt(x_moonFrame**2 + y_moonFrame**2 + z_moonFrame**2)

    v = v_0 * (r_fromMoon / R_0) ** 1.34 * x_moonFrame / r_fromMoon # Constant outwards mass flux

    return v

def v_OrbitalMotion(orbphase, key_scenario, x_steps, phi_steps, z_steps, ArchitectureDict, sigma_dim):
    # Calculate the velocity of a static atmosphere/exosphere around
    # the moving planet/exomoon

    M_s = ArchitectureDict['M_star']
    a_p = ArchitectureDict['a_p']

    v = -np.sin(orbphase) * np.sqrt(G * M_s / a_p)

    if key_scenario == 'exomoon':
        M_p = ArchitectureDict['M_p']
        a_moon = ArchitectureDict['a_moon']
        starting_orbphase_moon = ArchitectureDict['starting_orbphase_moon']

        orbphase_moon = starting_orbphase_moon + orbphase * np.sqrt((a_p**3 * M_p) / (a_moon**3 * M_s))

        v += - np.sin(orbphase_moon) * np.sqrt(G * M_p / a_moon)
    
    if sigma_dim == 5:

        v = np.tile(v, (int(x_steps), int(phi_steps), int(z_steps), 1))

    return v

def v_PlanetRotation(orbphase, key_scenario, x_steps, phi_steps, z_steps, ArchitectureDict, gridgrid):
    # Calculate the velocity of a static rotating atmosphere

    if key_scenario == 'barometric' or key_scenario == 'hydrostatic' or key_scenario == 'escaping':
        R_0 = ArchitectureDict['R_0']
        period_planetrot = ArchitectureDict['period_planetrot']

        xx, phiphi, rhorho, yy_pp = gridgrid

        v = 2 * np.pi / period_planetrot * (rhorho * np.sin(phiphi) - yy_pp)
    
    else:
        v = np.zeros((int(x_steps), int(phi_steps), int(z_steps), len(orbphase)))

    return v


def v_total(orbphase, key_scenario, x_steps, phi_steps, z_steps, ArchitectureDict, ScenarioDict, DopplerOrbitalMotion, DopplerPlanetRotation, sigma_dim, gridgrid):

    if sigma_dim == 2:
        v = np.zeros(len(orbphase))

    elif sigma_dim == 5: 
        v = np.zeros((int(x_steps), int(phi_steps), int(z_steps), len(orbphase)))

    if DopplerOrbitalMotion:

        v += v_OrbitalMotion(orbphase, key_scenario, x_steps, phi_steps, z_steps, ArchitectureDict, sigma_dim)

    if DopplerPlanetRotation:

        v += v_PlanetRotation(orbphase, key_scenario, x_steps, phi_steps, z_steps, ArchitectureDict, gridgrid)

    if key_scenario == 'escaping' and ScenarioDict['escaping']['RadialWind']:

        v += v_escaping(ArchitectureDict, ScenarioDict['escaping'], gridgrid)

    if key_scenario == 'exomoon' and ScenarioDict['exomoon']['RadialWind']:

        v += v_exomoon(orbphase, ArchitectureDict, ScenarioDict['exomoon'], gridgrid)


    return v

