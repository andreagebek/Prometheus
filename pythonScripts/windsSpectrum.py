"""
This file stores functions to calculate 
velocities of the absorbers, for the simulation
of a transmission spectrum.
Created on 8. September 2021 by Andrea Gebek.
"""

from constants import *
import numpy as np

def v_escaping(ArchitectureDict, EscapingDict, gridgrid, sigma_dim, sphericalSymmetry):
    # Radially escaping wind in the escaping scenario
    
    v_0 = EscapingDict['vRadial_0']
    q_esc = EscapingDict['q_esc']
    R_0 = ArchitectureDict['R_0']

    if sigma_dim == 4:

        xx, phiphi, rhorho = gridgrid
        r = np.sqrt(xx**2 + rhorho**2)


    elif sigma_dim == 3:

        if sphericalSymmetry:

            xx, rhorho = gridgrid
        
        else:

            xx, phiphi, rhorho = gridgrid
            xx = xx[:, 0, :]
            rhorho = rhorho[:, 0, :]
        
        r = np.sqrt(xx**2 + rhorho**2)

    v = v_0 * (r / R_0) ** (q_esc - 2.) * xx / r # Constant outwards mass flux

    return v

def v_exomoon(ArchitectureDict, ExomoonDict, gridgrid, sigma_dim, sphericalSymmetry):
    # Radially escaping wind in the exomoon scenario
    v_0 = ExomoonDict['vRadial_0']
    R_0 = ArchitectureDict['R_moon']

    if sigma_dim == 4:

        a_moon = ArchitectureDict['a_moon']
        z_moon = ArchitectureDict['z_moon']
        starting_orbphase_moon = ArchitectureDict['starting_orbphase_moon']

        x_moonFrame = xx - np.sqrt(a_moon**2 - z_moon**2) * np.cos(starting_orbphase_moon)
        y_moonFrame = np.sin(phiphi) * rhorho - np.sqrt(a_moon**2 - z_moon**2) * np.sin(starting_orbphase_moon)
        z_moonFrame = np.cos(phiphi) * rhorho - z_moon

        xx, phiphi, rhorho = gridgrid
        r = np.sqrt(x_moonFrame**2 + y_moonFrame**2 + z_moonFrame**2)

        v = v_0 * (r / R_0) ** 1.34 * x_moonFrame / r # Constant outwards mass flux

    elif sigma_dim == 3:

        if sphericalSymmetry:

            xx, rhorho = gridgrid
        
        else:

            xx, phiphi, rhorho = gridgrid
            xx = xx[:, 0, :]
            rhorho = rhorho[:, 0, :]
        
        r = np.sqrt(xx**2 + rhorho**2)

        v = v_0 * (r / R_0) ** 1.34 * xx / r # Constant outwards mass flux

    return v


def v_PlanetRotation(key_scenario, x_steps, phi_steps, z_steps, ArchitectureDict, gridgrid):
    # Calculate the velocity of a static rotating atmosphere

    if key_scenario == 'barometric' or key_scenario == 'hydrostatic' or key_scenario == 'escaping':
        R_0 = ArchitectureDict['R_0']
        period_planetrot = ArchitectureDict['period_planetrot']

        xx, phiphi, rhorho = gridgrid

        v = 2 * np.pi / period_planetrot * rhorho * np.sin(phiphi)
    
    else:
        v = np.zeros((int(x_steps), int(phi_steps), int(z_steps)))

    return v


def v_total(key_scenario, x_steps, z_steps, ArchitectureDict, ScenarioDict, DopplerPlanetRotation, sigma_dim, gridgrid, sphericalSymmetry = True, phi_steps = None):

    if sigma_dim == 3:
        v = np.zeros((int(x_steps), int(z_steps)))

    elif sigma_dim == 4:    # Only possible when the entire simulation is NOT spherically symmetric
        v = np.zeros((int(x_steps), int(phi_steps), int(z_steps)))

    if DopplerPlanetRotation:

        v += v_PlanetRotation(key_scenario, x_steps, phi_steps, z_steps, ArchitectureDict, gridgrid)

    if key_scenario == 'escaping' and ScenarioDict['escaping']['RadialWind']:

        v += v_escaping(ArchitectureDict, ScenarioDict['escaping'], gridgrid, sigma_dim, sphericalSymmetry)

    if key_scenario == 'exomoon' and ScenarioDict['exomoon']['RadialWind']:

        v += v_exomoon(ArchitectureDict, ScenarioDict['exomoon'], gridgrid, sigma_dim, sphericalSymmetry)

    return v

