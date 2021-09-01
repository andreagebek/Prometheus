"""
This file stores functions to calculate 
velocities of the absorbers.
Created on 1. September 2021 by Andrea Gebek.
"""

from constants import *

def v_intrinsic(orbphase, key_scenario, phi_steps, z_steps, params):
    # Calculate the velocity of a static atmosphere/exosphere around
    # the moving planet/exomoon

    M_s = params['M_star']
    a_p = params['a_p']

    v_intr = -np.sin(orbphase) * np.sqrt(G * M_s / a_p)

    if key_scenario == 'exomoon':

        M_p = params['M_p']
        a_moon = params['a_moon']
        v_intr += -np.sin(orbphase_moon) * np.sqrt(G * M_p / a_moon)

    v_intr = np.tile(v_intr, (int(phi_steps), int(z_steps), 1))

    return v_intr

