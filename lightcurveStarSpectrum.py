"""
Calculate the light curve of a transiting exoplanet,
taking into account the host star spectrum (e.g. for the RM-effect).
Created on 10. August 2021 by Andrea Gebek.
"""

import numpy as np
import json
import sys
from scipy.interpolate import interp1d
from datetime import datetime
from constants import *
from n_profiles import *
from sigma_abs import *
from stellarspectrum import *
import matplotlib.pyplot as plt



startTime = datetime.now()


"""
Calculate the host star spectrum on the spatial grid
"""

def CLV(rho, u1, u2):
    arg = 1. - np.sqrt(1. - rho**2 / R_s**2)
    return 1. - u1 * arg - u2 * arg**2


def Doppler(v):
    # If v is positive, receiver and source are moving towards each other
    beta = v / c
    shift = np.sqrt((1. - beta) / (1. + beta))

    return shift



def F_star(wavelength, phi, rho):
    PHOENIX_output = read_spectrum(T_eff, log_g, Fe_H, alpha_Fe)
    w_star = PHOENIX_output[0]
    v_max = 2. * np.pi * R_s / period_starrot
    w_max = np.max(wavelength * Doppler(-v_max))
    w_min = np.min(wavelength * Doppler(v_max))
    SEL_w = np.argwhere((w_star > w_min) * (w_star < w_max))[:, 0]

    SEL = np.concatenate((np.array([np.min(SEL_w) - 1]), SEL_w, np.array([np.max(SEL_w) + 1])))

    w_star = w_star[SEL]
    F_0 = PHOENIX_output[1][SEL]

    dir_omega = np.array([-np.sin(inclination_starrot) * np.cos(azimuth_starrot), -np.sin(inclination_starrot) * np.sin(azimuth_starrot), np.cos(inclination_starrot)])
    omega = 2. * np.pi / period_starrot * dir_omega # Angular velocity vector of the stellar rotation
    r_surface = np.array([np.tensordot(np.sqrt(R_s**2 - rho**2), np.ones(len(phi)), axes = 0), np.tensordot(rho, np.sin(phi), axes = 0), np.tensordot(rho, np.cos(phi), axes = 0)]) 
    # Vector to the surface of the star
    v_los = np.cross(omega, r_surface, axisb = 0)[:, :, 0] # The line-of-sight velocity is the one along the x-axis
    w_shift = Doppler(v_los)

    F_function = interp1d(w_star, F_0, kind = 'cubic')
    F_shifted = F_function(np.tensordot(wavelength, 1. / w_shift, axes = 0)) # This contains the shifted flux incident on the exoplanet, instead of shifting the spectra
    # at all positions on the stellar surface just read in the spectrum at different wavelengths depending on the position in the spatial grid (phi and rho) 

    if CLV_variations:
        CLV_array = np.repeat(CLV(rho, u1, u2), len(wavelength) * len(phi)).reshape(len(rho), len(wavelength), len(phi))
        
        return (F_shifted * CLV_array.swapaxes(0, 1)).swapaxes(1, 2)
    
    else:
        return F_shifted.swapaxes(1, 2)



"""
Calculation of the theoretical transit depth in two steps,
first calculating the (wavelength-dependent) optical depth along the 
chord at various impact parameters, then integrate over all possible 
impact parameters to obtain the (wavelength-dependent) transit depth
"""



def optical_depth(wavelength, phi, rho, x_p, y_p):

    x = np.linspace(-x_border, x_border, int(x_steps) + 1, dtype = np.dtype('f4'))[:-1] + x_border / float(x_steps)
    delta_x = 2 * x_border / float(x_steps)
    xx, phiphi, rhorho, yy_pp = np.meshgrid(x, phi, rho, y_p)

    r_fromP = np.sqrt(xx**2 + (rhorho * np.sin(phiphi) - yy_pp)**2 + (rhorho * np.cos(phiphi))**2)
    
    tau = 0
    for key_scenario in species_dict.keys():

        number_density = number_density_dict[key_scenario]        

        if key_scenario == 'barometric' or key_scenario == 'hydrostatic' or key_scenario == 'escaping':
        
            n = number_density(r_fromP, scenario_dict[key_scenario])

        elif key_scenario == 'exomoon':

            current_orbphase_moon = orbphase_moon + orbphase * np.sqrt((a_p**3 * M_p) / (a_moon**3 * M_s)) # Assuming Keplerian orbits of massless points

            x_moonFrame = xx - a_moon * np.cos(current_orbphase_moon)
            y_moonFrame = np.sin(phiphi) * rhorho - yy_pp - a_moon * np.sin(current_orbphase_moon)
            z_moonFrame = np.cos(phiphi) * rhorho

            r_fromMoon = np.sqrt(x_moonFrame**2 + y_moonFrame**2 + z_moonFrame**2)

            n = number_density(r_fromMoon, scenario_dict[key_scenario])

        elif key_scenario == 'torus':

            aa = np.sqrt(xx**2 + (np.sin(phiphi) * rhorho - yy_pp)**2)
            zz = np.cos(phiphi) * rhorho

            n = number_density(aa, zz, scenario_dict[key_scenario])

        yy = rhorho * np.sin(phiphi)
        zz = np.cos(phiphi) * rhorho
        blockingPlanet = (np.sqrt((yy - yy_pp)**2 + zz**2) < R_0)
        n = np.where(blockingPlanet, np.inf, n)

        if ExomoonSource: # Correct for chords which are blocked by the exomoon

            current_orbphase_moon = orbphase_moon + orbphase * np.sqrt((a_p**3 * M_p) / (a_moon**3 * M_s)) # Assuming Keplerian orbits of massless points
            y_moon = yy_pp + a_moon * np.sin(current_orbphase_moon)

            blockingMoon = ((yy - y_moon)**2 + zz**2 < R_moon**2) 

            n = np.where(blockingMoon, np.inf, n)
        
        behindStar = (xx + x_p < 0)
        n = np.where(behindStar, 0, n)

        N = delta_x * np.sum(n, axis = 1)

        sigma = 0
        for key_species in species_dict[key_scenario].keys():

            chi = species_dict[key_scenario][key_species]['chi'] # Mixing ratio OR number of absorbing atoms

            T_abs = species_dict[key_scenario][key_species]['T_abs']

            sigma += absorption_cross_section(wavelength, chi, T_abs, key_species, lines_dict)

        if 'RayleighScatt' in scenario_dict[key_scenario].keys():
            if scenario_dict[key_scenario]['RayleighScatt']:
                sigma += rayleigh_scattering(wavelength)

        tau += np.tensordot(sigma, N, axes = 0)

    return tau


def transit_depth(wavelength, orbphase):
    """Calculate the wavelength-dependent transit depth
    """

    phi = np.linspace(0, 2 * np.pi, int(phi_steps) + 1, dtype = np.dtype('f4'))[:-1] + np.pi / float(phi_steps)
    rho = np.linspace(0, R_s, int(z_steps) + 1, dtype = np.dtype('f4'))[:-1] + 0.5 * R_s / float(z_steps)

    x_p = a_p * np.cos(orbphase)
    y_p = np.array(a_p * np.sin(orbphase), dtype = np.dtype('f4'))
    
    single_chord = np.exp(-optical_depth(wavelength, phi, rho, x_p, y_p))

    delta_rho = R_s / float(z_steps)
    delta_phi = 2 * np.pi / float(phi_steps)

    F_input = np.moveaxis(np.tile(F_star(wavelength, phi, rho), (len(orbphase), 1, 1, 1)), 0, 3)
    # Make the stellar input spectrum four-dimensional: wavelength, phi, rho, orbphase

    integral_phi = delta_phi * np.sum(np.multiply(F_input, single_chord), axis = 1)

    sum_over_chords = delta_rho * np.tensordot(rho, integral_phi, axes = [0, 1])

    denominator_integralphi = delta_phi * np.sum(F_star(wavelength, phi, rho), axis = 1)
    denominator_integralrho = delta_rho * np.tensordot(rho, denominator_integralphi, axes = [0, 1])
    denominator = np.tile(denominator_integralrho, (len(orbphase), 1)).T
    
    return sum_over_chords / denominator




""".
Read in values from the setup file
"""

paramsFilename = sys.argv[1]

with open('../' + paramsFilename + '.txt') as file:
    param = json.load(file)

ExomoonSource = param['ExomoonSource']
CLV_variations = param['CLV_variations']


architecture_dict = param['Architecture']

R_s = architecture_dict['R_star']
M_s = architecture_dict['M_star']
R_0 = architecture_dict['R_0']
M_p = architecture_dict['M_p']
a_p = architecture_dict['a_p']
T_eff = architecture_dict['T_eff']
log_g = architecture_dict['log_g']
Fe_H = architecture_dict['Fe_H']
alpha_Fe = architecture_dict['alpha_Fe']
period_starrot = architecture_dict['period_starrot']
inclination_starrot = architecture_dict['inclination_starrot']
azimuth_starrot = architecture_dict['azimuth_starrot']

if ExomoonSource:
    R_moon = architecture_dict['R_moon']
    a_moon = architecture_dict['a_moon']
    orbphase_moon = architecture_dict['orbphase_moon']

if CLV_variations:
    u1 = architecture_dict['u1']
    u2 = architecture_dict['u2']


scenario_dict = param['Scenarios']
lines_dict = param['Lines']
species_dict = param['Species']

grids_dict = param['Grids']

wavelength = np.arange(grids_dict['lower_w'], grids_dict['upper_w'], grids_dict['resolution'])
orbphase_border = grids_dict['orbphase_border']
orbphase = np.linspace(-orbphase_border, orbphase_border, int(grids_dict['orbphase_steps']))
x_border = grids_dict['x_border']
x_steps = grids_dict['x_steps']
z_steps = grids_dict['z_steps']
phi_steps = grids_dict['phi_steps']

output_dict = param['Output']

number_density_dict = {'barometric': barometric, 'hydrostatic': hydrostatic, 'escaping': escaping, 'exomoon': exomoon, 'torus': torus}


"""
Execute the code and save the output as a .txt file
"""

header = 'Orbital phase'

for w in wavelength:
    header += ', Light curve at ' + str(w * 1e8) + ' Å'

lightcurve = transit_depth(wavelength, orbphase)
np.savetxt('../' + paramsFilename + '_lightcurve.txt', np.c_[orbphase / (2 * np.pi), lightcurve.T], header = header)

elapsed_time = datetime.now() - startTime

print("PROMETHEUS finished, yay! Elapsed time is:", elapsed_time)

print("The maximal flux decrease due to atmospheric/exospheric absorption in percent is:", np.abs(np.round(100 * (1 - np.min(lightcurve)), 5)))

print("The minimal flux decrease due to atmospheric/exospheric absorption in percent is:", np.abs(np.round(100 * (1 - np.max(lightcurve)), 5)))