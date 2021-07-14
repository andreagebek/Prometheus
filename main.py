"""
Run the main radiative transfer code.
Created on 3. June 2021 by Andrea Gebek.
"""

import numpy as np
import json
import sys
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
    T = params['T']
    P_0 = params['P_0']
    mu = params['mu']

    n_0 = P_0 / (k_B * T)
    H = k_B * T * R_0**2 / (G * mu * M_p)
    n = n_0 * np.exp((R_0 - r) / H)
    return n 


def hydrostatic(r, params):
    """Hydrostatic equilibrium with constant temperature and
    mean molecular weight
    """
    T = params['T']
    P_0 = params['P_0']
    mu = params['mu']

    n_0 = P_0 / (k_B * T)
    Jeans_0 = G * mu * M_p / (k_B * T * R_0)
    Jeans = G * mu * M_p / (k_B * T * r)
    n = n_0 * np.exp(Jeans - Jeans_0) 

    return n    



def escaping(r, params):
    """Power law with index q, for a hydrodynamically
    escaping atmosphere
    """
    q_esc = params['q_esc']

    if 'T' in params.keys():
        T = params['T']
        P_0 = params['P_0']
        n_0_esc = P_0 / (k_B * T)
    
    else:
        n_0_esc = (q_esc - 3.) / (4. * np.pi * R_0**3)

    n = n_0_esc * (R_0 / r)**q_esc

    return n


def exomoon(r, params):
    """Scaled to observed number density profile at Io (Burger 2001)
    """
    R_moon = params['R_moon']

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

    x = np.linspace(a_p - upper_x, a_p + upper_x, int(x_steps) + 1)[:-1] + upper_x / float(x_steps)

    delta_x = 2 * upper_x / float(x_steps)

    xx, zz = np.meshgrid(x, z)

    r = np.sqrt((xx-a_p)**2 + zz**2)

    tau = 0
    for key_scenario in species_dict.keys():

        number_density = number_density_dict[key_scenario]
        n = number_density(r, scenario_dict[key_scenario])

        N = delta_x * np.sum(n, axis = 1)

        sigma = 0
        for key_species in species_dict[key_scenario].keys():

            chi = species_dict[key_scenario][key_species]['chi'] # Mixing ratio OR number of absorbing atoms

            T_abs = species_dict[key_scenario][key_species]['T_abs']

            sigma += absorption_cross_section(wavelength, chi, T_abs, key_species)

        
        if 'RayleighScatt' in scenario_dict[key_scenario].keys():
            if scenario_dict[key_scenario]['RayleighScatt']:
                sigma += rayleigh_scattering(wavelength)

        tau += np.tensordot(sigma, N, axes = 0)
    
    return tau


def optical_depth3D(phi, rho, wavelength):

    x = np.linspace(a_p - upper_x, a_p + upper_x, int(x_steps) + 1)[:-1] + upper_x / float(x_steps)
    delta_x = 2 * upper_x / float(x_steps)
    xx, phiphi, rhorho = np.meshgrid(x, phi, rho)

    tau = 0
    for key_scenario in species_dict.keys():

        number_density = number_density_dict[key_scenario]

        if key_scenario == 'barometric' or key_scenario == 'hydrostatic' or key_scenario == 'escaping':
        
            r = np.sqrt((xx - a_p)**2 + rhorho**2)

        elif key_scenario == 'exomoon':
            x_moonFrame = xx - a_p + np.sqrt(a_moon**2 - z_moon**2) * np.sin(alpha_moon)
            y_moonFrame = np.sin(phiphi) * rhorho - np.sqrt(a_moon**2 - z_moon**2) * np.cos(alpha_moon)
            z_moonFrame = np.cos(phiphi) * rhorho - z_moon

            r = np.sqrt(x_moonFrame**2 + y_moonFrame**2 + z_moonFrame**2)

        n = number_density(r, scenario_dict[key_scenario])

        if ExomoonSource: # Correct for chords which are blocked by the off-center exomoon (but NOT by the exoplanet)

            y_moon = np.sqrt(a_moon**2-z_moon**2) * np.cos(alpha_moon)
            yy = np.sin(phiphi) * rhorho
            zz = np.cos(phiphi) * rhorho

            blockingMoon = ((yy-y_moon)**2 + (zz-z_moon)**2 < R_moon**2) * (rhorho > R_0)
            
            n = np.where(blockingMoon, 0, n)

        N = delta_x * np.sum(n, axis = 1)

        sigma = 0
        for key_species in species_dict[key_scenario].keys():

            chi = species_dict[key_scenario][key_species]['chi'] # Mixing ratio OR number of absorbing atoms

            T_abs = species_dict[key_scenario][key_species]['T_abs']

            sigma += absorption_cross_section(wavelength, chi, T_abs, key_species)

        if 'RayleighScatt' in scenario_dict[key_scenario].keys():
            if scenario_dict[key_scenario]['RayleighScatt']:
                sigma += rayleigh_scattering(wavelength)

        tau += np.tensordot(sigma, N, axes = 0)

    return tau





def transit_depth(wavelength):
    """Calculate the wavelength-dependent transit depth
    """
    if sphericalSymmetry:

        if ExomoonSource:
            starting_z = R_moon
        
        else:
            starting_z = R_0

        z = np.linspace(starting_z, R_s, int(z_steps) + 1)[:-1] + 0.5 * (R_s - starting_z) / float(z_steps)

        single_chord = np.exp(-optical_depth(z, wavelength))
            
        delta_z = (R_s - starting_z) / float(z_steps)

        sum_over_chords = delta_z * np.tensordot(z, single_chord, axes = [0, 1])
    
        return 2 * sum_over_chords / R_s**2

    else:

        phi = np.linspace(0, 2 * np.pi, int(phi_steps) + 1)[:-1] + np.pi / float(phi_steps)
        rho = np.linspace(R_0, R_s, int(z_steps) + 1)[:-1] + 0.5 * (R_s - R_0) / float(z_steps)
        
        single_chord = np.exp(-optical_depth3D(phi, rho, wavelength))

        delta_rho = (R_s - R_0) / float(z_steps)
        delta_phi = 2 * np.pi / float(phi_steps)

        integral_phi = delta_phi * np.sum(single_chord, axis = 1)

        sum_over_chords = delta_rho * np.tensordot(rho, integral_phi, axes = [0, 1])

        if ExomoonSource:
            sum_over_chords -= np.pi * R_moon**2

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
    
        sigma += absorption_cross_section(wavelength, chi, T_abs, key_species)
    
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

sphericalSymmetry = param['sphericalSymmetry']
ExomoonSource = param['ExomoonSource']
ExomoonOffCenter = param['ExomoonOffCenter']


architecture_dict = param['Architecture']

R_s = architecture_dict['R_star']
R_0 = architecture_dict['R_0']
M_p = architecture_dict['M_p']
a_p = architecture_dict['a_p']

if ExomoonSource:
    R_moon = architecture_dict['R_moon']

if ExomoonOffCenter:
    a_moon = architecture_dict['a_moon']
    alpha_moon = architecture_dict['alpha_moon']
    z_moon = architecture_dict['z_moon']


scenario_dict = param['Scenarios']
lines_dict = param['Lines']
species_dict = param['Species']

grids_dict = param['Grids']

wavelength = np.linspace(grids_dict['lower_w'], grids_dict['upper_w'], 1 + int((grids_dict['upper_w'] - grids_dict['lower_w']) / grids_dict['resolution']))
upper_x = grids_dict['upper_x']
x_steps = grids_dict['x_steps']
z_steps = grids_dict['z_steps']

if not sphericalSymmetry:
    phi_steps = grids_dict['phi_steps']

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

print("PROMETHEUS finished, yay! Elapsed time is:")
print(elapsed_time)
print("The maximal flux decrease due to atmospheric/exospheric absorption in percent is:")
print(np.round(100*(1-min(spectrum)), 5))
print("The minimal flux decrease due to atmospheric/exospheric absorption in percent is:")
print(np.round(100*(1-max(spectrum)), 5))

if benchmark:
    benchmark_spectrum = BarometricBenchmark(scenario_dict['barometric'])

    np.savetxt('../' + paramsFilename + '_spectrum.txt', np.c_[wavelength, spectrum, benchmark_spectrum], 
    header = 'Wavelength [Angstrom], Transit depth spectrum, Transit depth benchmark')
    
    print("The maximal BENCHMARK flux decrease due to atmospheric/exospheric absorption in percent is:")
    print(np.round(100*(1-min(benchmark_spectrum)), 5))
    print("The minimal BENCHMARK flux decrease due to atmospheric/exospheric absorption in percent is:")
    print(np.round(100*(1-max(benchmark_spectrum)), 5))

if record_tau:
    idx_maxabs = np.argmin(spectrum)
    w_maxabs = wavelength[idx_maxabs]

    if sphericalSymmetry:

        if ExomoonSource:
            starting_z = R_moon
        
        else:
            starting_z = R_0

        z = np.linspace(starting_z, R_s, int(z_steps) + 1)[:-1] + 0.5 * (R_s - starting_z) / float(z_steps)

        tau = optical_depth(z, w_maxabs)

        np.savetxt('../' + paramsFilename + '_tau.txt', np.array([z, tau]).T, header = 'z [cm], tau')        
    
    else:

        phi = np.linspace(0, 2 * np.pi, int(phi_steps) + 1)[:-1] + np.pi / float(phi_steps)
        rho = np.linspace(R_0, R_s, int(z_steps) + 1)[:-1] + 0.5 * (R_s - R_0) / float(z_steps)

        phiphi, rhorho = np.meshgrid(phi, rho, indexing = 'ij')
        
        phiphi = phiphi.reshape(int(phi_steps * z_steps))
        rhorho = rhorho.reshape(int(phi_steps * z_steps))

        tau = optical_depth3D(phi, rho, w_maxabs).reshape(int(phi_steps * z_steps))

        np.savetxt('../' + paramsFilename + '_tau.txt', np.array([phiphi, rhorho, tau]).T, header = 'phi grid [rad], rho grid [cm], tau')
    
    