"""
Calculate the transit spectrum for a spherically symmetric atmosphere/exosphere.
Created on 3. June 2021 by Andrea Gebek.
"""

import numpy as np
import json
import sys
from datetime import datetime
from constants import *
from n_profiles import *
from sigma_abs import *
from windsSpectrum import *


startTime = datetime.now()

"""
Doppler shift
"""

def Doppler(v):
    # If v is positive, receiver and source are moving towards each other
    beta = v / c
    shift = np.sqrt((1. - beta) / (1. + beta))

    return shift



"""
Calculation of the theoretical transit depth in two steps,
first calculating the (wavelength-dependent) optical depth along the 
chord at various impact parameters, then integrate over all possible 
impact parameters to obtain the (wavelength-dependent) transit depth
"""

def optical_depth(wavelength, z):   
    """Calculate the optical depth for all absorbers along the chord 
    """

    x = np.linspace(-x_border, x_border, int(x_steps) + 1)[:-1] + x_border / float(x_steps)

    delta_x = 2. * x_border / float(x_steps) 

    gridgrid = np.meshgrid(x, z, indexing = 'ij')
    xx, zz = gridgrid

    r = np.sqrt(xx**2 + zz**2)

    tau = 0
    for key_scenario in species_dict.keys():

        number_density = number_density_dict[key_scenario]
        n = number_density(r, scenario_dict[key_scenario])

        if DopplerPlanetRotation or RadialWinds:

            v_los = v_total(key_scenario, x_steps, z_steps, architecture_dict, scenario_dict, DopplerPlanetRotation, sigma_dim, gridgrid)
            w_shift = Doppler(-v_los)
            wavelength_shifted = np.tensordot(wavelength, w_shift, axes = 0)
        
        else:
            wavelength_shifted = wavelength       

        sigma = 0
        for key_species in species_dict[key_scenario].keys():

            chi = species_dict[key_scenario][key_species]['chi'] # Mixing ratio OR number of absorbing atoms

            T_abs = species_dict[key_scenario][key_species]['T_abs']

            sigma += absorption_cross_section(wavelength_shifted, chi, T_abs, key_species, lines_dict)

        
        if 'RayleighScatt' in scenario_dict[key_scenario].keys():
            if scenario_dict[key_scenario]['RayleighScatt']:
                sigma += rayleigh_scattering(wavelength_shifted)

        if sigma_dim == 1:  # sigma(lambda)

            N = delta_x * np.sum(n, axis = 0)
            tau += np.tensordot(sigma, N, axes = 0)

        elif sigma_dim == 3:    # sigma(lambda, x, rho)

            n = np.tile(n, (len(wavelength), 1, 1))

            tau += delta_x * np.sum(np.multiply(sigma, n), axis = 1)
    
    return tau



def transit_depth(wavelength):
    """Calculate the wavelength-dependent transit depth
    """

    if ExomoonSource:
        starting_z = R_moon
    
    else:
        starting_z = R_0

    z = np.linspace(starting_z, R_s, int(z_steps) + 1)[:-1] + 0.5 * (R_s - starting_z) / float(z_steps)

    single_chord = np.exp(-optical_depth(wavelength, z))
        
    delta_z = (R_s - starting_z) / float(z_steps)

    sum_over_chords = delta_z * np.tensordot(z, single_chord, axes = [0, 1])

    return 2 * sum_over_chords / R_s**2


"""
Additional Functionality
"""

def BarometricBenchmark(params):
    T = params['T']
    P_0 = params['P_0']
    mu = params['mu']

    H = k_B * T * R_0**2 / (G * mu * M_p)
    g = G * M_p / R_0**2

    sigma = 0
    for key_species in species_dict['barometric'].keys():

        chi = species_dict['barometric'][key_species]['chi']
        T_abs = species_dict['barometric'][key_species]['T_abs']
    
        sigma += absorption_cross_section(wavelength, chi, T_abs, key_species, lines_dict)
    
    if scenario_dict['barometric']['RayleighScatt']:
            sigma += rayleigh_scattering(wavelength)

    kappa = sigma / mu
    R_lambda = R_0 + H * (euler_mascheroni + np.log(P_0 * kappa / g * np.sqrt(2 * np.pi * R_0 / H)))
    benchmark_spectrum = (R_s**2 - R_lambda**2) / R_s**2

    return benchmark_spectrum


""".
Read in values from the setup file
"""

paramsFilename = sys.argv[1]

with open('../' + paramsFilename + '.txt') as file:
    param = json.load(file)

ExomoonSource = param['ExomoonSource']
DopplerPlanetRotation = param['DopplerPlanetRotation']
RadialWinds = param['RadialWinds']
sigma_dim = param['sigma_dim']

architecture_dict = param['Architecture']

R_s = architecture_dict['R_star']
R_0 = architecture_dict['R_0']
M_p = architecture_dict['M_p']

if ExomoonSource:
    R_moon = architecture_dict['R_moon']

scenario_dict = param['Scenarios']
lines_dict = param['Lines']
species_dict = param['Species']

grids_dict = param['Grids']

wavelength = np.arange(grids_dict['lower_w'], grids_dict['upper_w'], grids_dict['resolution'])
x_border = grids_dict['x_border']
x_steps = grids_dict['x_steps']
z_steps = grids_dict['z_steps']


output_dict = param['Output']

benchmark = output_dict['benchmark']
record_tau = output_dict['record_tau']


number_density_dict = {'barometric': barometric, 'hydrostatic': hydrostatic, 'escaping': escaping, 'exomoon': exomoon}


"""
Execute the code and save the output as a .txt file
"""

spectrum = transit_depth(wavelength)
np.savetxt('../' + paramsFilename + '_spectrum.txt', np.c_[wavelength, spectrum], header = 'Wavelength [Angstrom], Transit depth spectrum')

elapsed_time = datetime.now() - startTime

print("PROMETHEUS finished, yay! Elapsed time is:", elapsed_time)

print("The maximal flux decrease due to atmospheric/exospheric absorption in percent is:", np.round(100*(1-min(spectrum)), 5))

print("The minimal flux decrease due to atmospheric/exospheric absorption in percent is:", np.round(100*(1-max(spectrum)), 5))

if benchmark:
    benchmark_spectrum = BarometricBenchmark(scenario_dict['barometric'])

    np.savetxt('../' + paramsFilename + '_spectrum.txt', np.c_[wavelength, spectrum, benchmark_spectrum], 
    header = 'Wavelength [Angstrom], Transit depth spectrum, Transit depth benchmark')
    
    print("The maximal BENCHMARK flux decrease due to atmospheric/exospheric absorption in percent is:", np.round(100*(1-min(benchmark_spectrum)), 5))

    print("The minimal BENCHMARK flux decrease due to atmospheric/exospheric absorption in percent is:", np.round(100*(1-max(benchmark_spectrum)), 5))


if record_tau:
    idx_maxabs = np.argmin(spectrum)

    if ExomoonSource:
        starting_z = R_moon
    
    else:
        starting_z = R_0

    z = np.linspace(starting_z, R_s, int(z_steps) + 1)[:-1] + 0.5 * (R_s - starting_z) / float(z_steps)

    tau = optical_depth(wavelength, z)[idx_maxabs, :]

    np.savetxt('../' + paramsFilename + '_tau.txt', np.array([z, tau]).T, header = 'z [cm], tau')        
    

    
    