"""
This file stores functions to calculate 
velocities of the absorbers.
Created on 1. September 2021 by Andrea Gebek.
"""

from constants import *

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


def v_total(orbphase, key_scenario, x_steps, phi_steps, z_steps, ArchitectureDict, DopplerOrbitalMotion, DopplerPlanetRotation, sigma_dim, gridgrid):

    if sigma_dim == 2:
        v = np.zeros(len(orbphase))

    elif sigma_dim == 5: 
        v = np.zeros((int(x_steps), int(phi_steps), int(z_steps), len(orbphase)))

    if DopplerOrbitalMotion:

        v += v_OrbitalMotion(orbphase, key_scenario, x_steps, phi_steps, z_steps, ArchitectureDict, sigma_dim)

    if DopplerPlanetRotation: # We always have sigma_dim == 5 if DopplerPlanetRotation == True

        v += v_PlanetRotation(orbphase, key_scenario, x_steps, phi_steps, z_steps, ArchitectureDict, gridgrid)

    return v

