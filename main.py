"""
Run the main radiative transfer code.
Created on 3. June 2021 by Andrea Gebek.
"""

import numpy as np
import json
from scipy.special import wofz
from datetime import datetime
from constants import *


startTime = datetime.now()

"""
Number density functions
"""

def barometric(r, params):
    """Barometric law (to benchmark the numerical solution against analytical
    results)
    """
    T = params[0]
    P_0 = params[1]
    mu = params[2]

    n_0 = P_0 / (k_B * T)
    H = k_B * T * R_0**2 / (G * mu * M_p)
    n = n_0 * np.exp((R_0 - r) / H)
    return n 


def hydrostatic(r, params):
    """Hydrostatic equilibrium with constant temperature and
    mean molecular weight
    """
    T = params[0]
    P_0 = params[1]
    mu = params[2]

    n_0 = P_0 / (k_B * T)
    Jeans_0 = G * mu * M_p / (k_B * T * R_0)
    Jeans = G * mu * M_p / (k_B * T * r)
    n = n_0 * np.exp(Jeans - Jeans_0) 

    return n    



def escaping(r, params):
    """Power law with index q, for a hydrodynamically
    escaping atmosphere
    """
    q = params[0]
    norm_esc = params[1]

    if norm_esc == 'pressure':
        T = params[2]
        P_0 = params[3]
        n_0_esc = P_0 / (k_B * T)
    
    else:
        n_0_esc = (q - 3.) / (4. * np.pi * R_0**3)

    n = n_0_esc * (R_0 / r)**q

    return n


def exomoon(r, params):
    """Scaled to observed number density profile at Io (Burger 2001)
    """
    R_moon = params[0]

    n_0_exomoon = 0.34 / (4 * np.pi * R_moon**3)
    n = n_0_exomoon * (r / R_moon)**(-3.34)

    return n


def central_torus_density(N):
    """ Calculates the number density in the center of the torus
    """
    v_orbit = np.sqrt(G * M_p / a_torus)
    H_torus = a_torus * v_ej / v_orbit
    
    z_grid = np.linspace(0, 4 * H_torus, 100)
    delta_z = 4. * H_torus / 100.
    
    a_grid = np.linspace(a_torus, a_torus + 16 * H_torus, 100)
    delta_a = 16. * H_torus / 100.
    
    n_a = 2 * np.pi * a_grid * 4 * np.exp(-((a_grid - a_torus)/(4 * H_torus))**2)
    n_z = np.exp(-(z_grid / H_torus)**2)    
    n = np.tensordot(n_a, n_z, axes = 0) * delta_z * delta_a

    return N / np.sum(n)

def torus(a, z):
    """ Toroidal atmosphere (Johnson & Huggins 2006)
    """
    n_0_torus = central_torus_density(N)
    v_orbit = np.sqrt(G * M_p / a_torus)
    H_torus = a_torus * v_ej / v_orbit
    n_a = np.exp(- ((a - a_torus) / (4 * H_torus))**2)
    n_z = np.exp(- (z / H_torus)**2)
    n = n_0_torus * np.multiply(n_a, n_z) 
    
    return n


"""
Calculation of the absorption cross section from the line list
"""



def line_absorption(wavelength, T, absorber_mass,
                  line_wavelength, line_intensity, line_HWHM):
    """Voigt or Lorentz profile
    """
    x = 1 / wavelength - 1 / line_wavelength # Frequency difference, converted to wavenumber
 
    if T > 0:
        sigma_v = 9.12 * np.sqrt(T * amu / (1e4 * absorber_mass)) * 1e5 / (c * line_wavelength)
        #Doppler broadening, converted to wavenumber
        argument = (x + 1j * line_HWHM / c) / (sigma_v * np.sqrt(2)) 
        phi_line = np.real(wofz(argument)) / (sigma_v * np.sqrt(2 * np.pi))
    
    else:
        phi_line = (line_HWHM / (c * np.pi)) / (x**2 + line_HWHM**2 / c**2) # Lorentz broadening, converted to wavenumber

    return phi_line * line_intensity



def absorption_cross_section(wavelength, chi, T, species):
    """Compute the total absorption cross section from all absorbing lines
    """

    sigma = 0

    for value in param['Lines'].values():

        if value[4] != species:
            continue

        line_intensity = value[0]
        line_wavelength = value[1]
        line_HWHM = value[2]
        line_mass = value[3]

        sigma = sigma + line_absorption(wavelength, T, line_mass, line_wavelength,
                                        line_intensity, line_HWHM) * chi

    return sigma


def rayleigh_scattering(wavelength):
    sigma = 8.49e-45 / wavelength**4 # FROM WHERE IS THIS? DEPENDENCY ON H2 mixing ratio?
    return sigma


"""
Calculation of the theoretical transit depth in two steps,
first calculating the (wavelength-dependent) optical depth along the 
chord at various impact parameters, then integrate over all possible 
impact parameters to obtain the (wavelength-dependent) transit depth
"""

def optical_depth(z, wavelength):   
    """Calculate the optical depth for all absorbers along the chord 
    """

    x = np.linspace(0, R_s, int(x_steps) + 1)[:-1] + 0.5 * R_s / float(x_steps)

    delta_x = R_s / float(x_steps)

    xx, zz = np.meshgrid(x, z)

    r = np.sqrt(xx**2 + zz**2)

    tau = 0
    for idx_scenario, key_scenario in enumerate(scenario_dict.keys()):
        number_density = number_density_dict[key_scenario]
        n = number_density(r, scenario_dict[key_scenario])

        N = 2 * delta_x * np.sum(n, axis = 1)

        sigma = 0
        for key_species in species_dict.keys():

            chi = species_dict[key_species][idx_scenario][1] # Mixing ratio OR number of absorbing atoms

            if not check('evaporative', key_scenario, scenario_dict[key_scenario]):
                T = species_dict[key_species][idx_scenario][2]


            else:

                if 'therm_broad' in scenario_dict[key_scenario]:
                    v_mean = species_dict[key_species][idx_scenario][2]
                    absorber_mass = species_dict[key_species][idx_scenario][3]
                    T = v_mean**2 * absorber_mass * np.pi / (8. * k_B)
                
                else:
                    T = 0   # No thermal broadening

            sigma += absorption_cross_section(wavelength, chi, T, key_species)

        
        if 'rayleigh_scatt' in scenario_dict[key_scenario]:
            sigma += rayleigh_scattering(wavelength)

        tau += np.tensordot(sigma, N, axes = 0)

    return tau


def transit_depth(wavelength):
    """Calculate the wavelength-dependent transit depth
    """

    
    z = np.linspace(R_0, R_s, int(z_steps) + 1)[:-1] + 0.5 * (R_s - R_0) / float(z_steps)

    single_chord = np.exp(-optical_depth(z, wavelength))
        
    delta_z = (R_s - R_0) / float(z_steps)

    sum_over_chords = delta_z * np.tensordot(z, single_chord, axes = [0, 1])
    
    return 2. / (R_s**2 - R_0**2) * sum_over_chords


"""
Additional Functionality
"""

def BarometricBenchmark(params):
    T = params[0]
    P_0 = params[1]
    mu = params[2]

    H = k_B * T * R_0**2 / (G * mu * M_p)
    g = G * M_p / R_0**2

    sigma = 0
    for key_species in param['Species'].keys():
        for parameter_list in param['Species'][key_species]:
            if parameter_list[0] == 'barometric':
                chi = parameter_list[1]
                T = parameter_list[2]
            else:
                continue
    
        sigma += absorption_cross_section(wavelength, chi, T, key_species)

    kappa = sigma / mu
    R = R_0 + H * (euler_mascheroni + np.log(P_0 * kappa / g * np.sqrt(2 * np.pi * R_0 / H)))
    benchmark_spectrum = 1 - (R**2 - R_0**2) / (R_s**2 - R_0**2)

    return benchmark_spectrum



""".
Read in values from the setup file
"""

with open('../settings.txt') as file:
    param = json.load(file)



R_s = param['R_star']
R_0 = param['R_0']
M_p = param['M_p']

scenario_dict = param['Scenarios']
lines_dict = param['Lines']
species_dict = param['Species']


wavelength = np.linspace(param['lower_w'], param['upper_w'], 1 + int((param['upper_w'] - param['lower_w']) / param['resolution']))
x_steps = param['x_steps']
z_steps = param['z_steps']

benchmark_ON = param['benchmark']


number_density_dict = {'barometric': barometric, 'hydrostatic': hydrostatic, 'escaping': escaping, 'exomoon': exomoon}




"""
Execute the code and save the output as a .txt file
"""

spectrum = transit_depth(wavelength)
np.savetxt('../test.txt', np.c_[wavelength, spectrum])

elapsed_time = datetime.now() - startTime

print("PROMETHEUS finished, yay! Elapsed time is:")
print(elapsed_time)
print("The maximal flux decrease due to atmospheric/exospheric absorption in percent is:")
print(np.round(100*(1-min(spectrum)), 5))
print("The minimal flux decrease due to atmospheric/exospheric absorption in percent is:")
print(np.round(100*(1-max(spectrum)), 5))

if benchmark_ON == 'yes':
    benchmark_spectrum = BarometricBenchmark(param['Scenarios']['barometric'])

    np.savetxt('../test.txt', np.c_[wavelength, spectrum, benchmark_spectrum])
    
    print("The maximal BENCHMARK flux decrease due to atmospheric/exospheric absorption in percent is:")
    print(np.round(100*(1-min(benchmark_spectrum)), 5))
    print("The minimal BENCHMARK flux decrease due to atmospheric/exospheric absorption in percent is:")
    print(np.round(100*(1-max(benchmark_spectrum)), 5))