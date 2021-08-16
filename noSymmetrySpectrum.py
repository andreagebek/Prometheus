"""
Calculate the transit spectrum for a nonsymmetric atmosphere/exosphere in full 3D.
Created on 3. June 2021 by Andrea Gebek.
"""

import numpy as np
import json
import sys
from datetime import datetime
from constants import *
from n_profiles import *
from sigma_abs import *

startTime = datetime.now()

"""
Helper function to calculate the area of two intersecting circles
"""


def A_intersect(r1, r2, d):

    if np.abs(d) > (r1 + r2):
        A = 0
    
    elif np.abs(d) < np.abs(r1 - r2):
        A = np.pi * np.min((r1, r2))**2
    
    else:
        term1 = (r1**2 - r2**2 + d**2) / (2. * np.abs(d))
        term2 = (r2**2 - r1**2 + d**2) / (2. * np.abs(d))
        term3 = r1**2 * np.arccos(term1 / r1) - term1 * np.sqrt(r1**2 - term1**2)
        term4 = r2**2 * np.arccos(term2 / r2) - term2 * np.sqrt(r2**2 - term2**2)

        A = term3 + term4
  
    return A


"""
Calculation of the theoretical transit depth in two steps,
first calculating the (wavelength-dependent) optical depth along the 
chord at various impact parameters, then integrate over all possible 
impact parameters to obtain the (wavelength-dependent) transit depth
"""




def optical_depth(wavelength, phi, rho):

    x = np.linspace(-x_border, x_border, int(x_steps) + 1)[:-1] + x_border / float(x_steps)
    delta_x = 2 * x_border / float(x_steps)
    xx, phiphi, rhorho = np.meshgrid(x, phi, rho)

    tau = 0
    for key_scenario in species_dict.keys():

        number_density = number_density_dict[key_scenario]

        if key_scenario == 'barometric' or key_scenario == 'hydrostatic' or key_scenario == 'escaping':
        
            r = np.sqrt(xx**2 + rhorho**2)

            n = number_density(r, scenario_dict[key_scenario])

        elif key_scenario == 'exomoon':
            x_moonFrame = xx - np.sqrt(a_moon**2 - z_moon**2) * np.cos(orbphase_moon)
            y_moonFrame = np.sin(phiphi) * rhorho - np.sqrt(a_moon**2 - z_moon**2) * np.sin(orbphase_moon)
            z_moonFrame = np.cos(phiphi) * rhorho - z_moon

            r = np.sqrt(x_moonFrame**2 + y_moonFrame**2 + z_moonFrame**2)

            n = number_density(r, scenario_dict[key_scenario])

        elif key_scenario == 'torus':

            aa = np.sqrt(xx**2 + (np.sin(phiphi) * rhorho)**2)
            zz = np.cos(phiphi) * rhorho
            
            n = number_density(aa, zz, scenario_dict[key_scenario])
        

        if ExomoonSource: # Correct for chords which are blocked by the off-center exomoon 

            y_moon = np.sqrt(a_moon**2 - z_moon**2) * np.sin(orbphase_moon)
            yy = np.sin(phiphi) * rhorho
            zz = np.cos(phiphi) * rhorho

            blockingMoon = ((yy - y_moon)**2 + (zz - z_moon)**2 < R_moon**2)
            
            n = np.where(blockingMoon, 0, n)

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





def transit_depth(wavelength):
    """Calculate the wavelength-dependent transit depth
    """

    phi = np.linspace(0, 2 * np.pi, int(phi_steps) + 1)[:-1] + np.pi / float(phi_steps)
    rho = np.linspace(R_0, R_s, int(z_steps) + 1)[:-1] + 0.5 * (R_s - R_0) / float(z_steps)
    
    single_chord = np.exp(-optical_depth(wavelength, phi, rho))

    delta_rho = (R_s - R_0) / float(z_steps)
    delta_phi = 2 * np.pi / float(phi_steps)

    integral_phi = delta_phi * np.sum(single_chord, axis = 1)

    sum_over_chords = delta_rho * np.tensordot(rho, integral_phi, axes = [0, 1])

    if ExomoonSource:
        y_moon = np.sqrt(a_moon**2 - z_moon**2) * np.sin(orbphase_moon)
        A_occ_moon = np.pi * R_moon**2 - A_intersect(R_moon, R_0, np.sqrt(y_moon**2 + z_moon**2)) # This is the area of the moon which is not blocked by the exoplanet
        sum_over_chords -= A_occ_moon

    return sum_over_chords / (np.pi * R_s**2)

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



architecture_dict = param['Architecture']

R_s = architecture_dict['R_star']
R_0 = architecture_dict['R_0']
M_p = architecture_dict['M_p']

if ExomoonSource:
    R_moon = architecture_dict['R_moon']
    a_moon = architecture_dict['a_moon']
    orbphase_moon = architecture_dict['orbphase_moon']
    z_moon = architecture_dict['z_moon']


scenario_dict = param['Scenarios']
lines_dict = param['Lines']
species_dict = param['Species']

grids_dict = param['Grids']

wavelength = np.linspace(grids_dict['lower_w'], grids_dict['upper_w'], 1 + int((grids_dict['upper_w'] - grids_dict['lower_w']) / grids_dict['resolution']))
x_border = grids_dict['x_border']
x_steps = grids_dict['x_steps']
z_steps = grids_dict['z_steps']
phi_steps = grids_dict['phi_steps']

output_dict = param['Output']

benchmark = output_dict['benchmark']
record_tau = output_dict['record_tau']


number_density_dict = {'barometric': barometric, 'hydrostatic': hydrostatic, 'escaping': escaping, 'exomoon': exomoon, 'torus': torus}


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
    w_maxabs = wavelength[idx_maxabs]

    phi = np.linspace(0, 2 * np.pi, int(phi_steps) + 1)[:-1] + np.pi / float(phi_steps)
    rho = np.linspace(R_0, R_s, int(z_steps) + 1)[:-1] + 0.5 * (R_s - R_0) / float(z_steps)

    phiphi, rhorho = np.meshgrid(phi, rho, indexing = 'ij')
    
    phiphi = phiphi.reshape(int(phi_steps * z_steps))
    rhorho = rhorho.reshape(int(phi_steps * z_steps))

    tau = optical_depth(w_maxabs, phi, rho).reshape(int(phi_steps * z_steps))

    np.savetxt('../' + paramsFilename + '_tau.txt', np.array([phiphi, rhorho, tau]).T, header = 'phi grid [rad], rho grid [cm], tau')   
    