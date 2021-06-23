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


def barometric(R):
    """Barometric law (to benchmark the numerical solution against analytical
    results)
    Input: Radial coordinate (2d-array, cm)
    Output: Number density (2d-array, 1/cm^3)
    """
    n_0 = P_0/(k_B*T)
    H = k_B*T*R_0**2/(G*mu*M_p)
    n = n_0*np.exp((R_0-R)/H)
    return n 


def hydrostatic(R):
    """Hydrostatic equilibrium with constant temperature and
    mean molecular weight
    Input: Radial coordinate (2d-array, cm)
    Output: Number density (2d-array, 1/cm^3)
    """
    n_0 = P_0/(k_B*T)
    Jeans_0 = G*mu*M_p/(k_B*T*R_0)
    Jeans = G*mu*M_p/(k_B*T*R)
    n = n_0*np.exp(Jeans-Jeans_0)  
    return n    



def escaping(R):
    """Power law with index q,
    indicative of a hydrodynamically escaping atmosphere
    Input: Radial coordinate (2d-array, cm)
    Output: Number density (2d-array, 1/cm^3)
    """
    n_0_esc = N*(q-3.)/(4.*pi*R_0**3)
    n = n_0_esc*(R_0/R)**q
    return n


def exomoon(R):
    """Scaled to observed number density profile at Io (Burger 2001)
    Input: Radial coordinate (2d-array, cm)
    Output: Number density (2d-array, 1/cm^3)
    """
    n_0_exomoon = N*0.34/(4*pi*R_moon**3)
    n = n_0_exomoon*(R/R_moon)**(-3.34)
    return n


def central_torus_density(N):
    """ Calculates the number density in the center of the torus
    Input: Number of particles 
    Output: Number density (1/cm^3)
    """
    v_orbit = np.sqrt(G*M_p/a_torus)
    H_torus = a_torus*v_ej/v_orbit
    
    z_grid = np.linspace(0,4*H_torus,100)
    delta_z = 4.*H_torus/100.
    
    a_grid = np.linspace(a_torus, a_torus+16*H_torus, 100)
    delta_a = 16.*H_torus/100.
    
    n_a = 2*pi*a_grid*4*np.exp(-((a_grid-a_torus)/(4*H_torus))**2)
    n_z = np.exp(-(z_grid/H_torus)**2)    
    n = np.tensordot(n_a, n_z, axes = 0)*delta_z*delta_a

    return N/np.sum(n)

def torus(a,z):
    """ Toroidal atmosphere (Johnson & Huggins 2006)
    Input: Radial-planar coordinate (cm), vertical coordinate (cm)
    Output: Number density (1/cm^3)
    """
    n_0_torus = central_torus_density(N)
    v_orbit = np.sqrt(G * M_p / a_torus)
    H_torus = a_torus * v_ej / v_orbit
    n_a = np.exp(- ((a - a_torus) / (4 * H_torus))**2)
    n_z = np.exp(- (z / H_torus)**2)
    n = n_0_torus*np.multiply(n_a, n_z) 
    
    return n


"""
Calculation of the absorption cross section from the line list
"""



def line_absorption(wavelength, T, absorber_mass,
                  line_wavelength, line_intensity, line_HWHM):
    """Voigt profile
    Input: Wavelength (array, cm), Temperature (scalar, K), Absorber mass (scalar, g),
    Wavelength of transition (scalar, cm), Line intensity (scalar, cm/particle),
    HWHM of the intrinsic Lorentzian line profile (scalar, 1/s)
    Output: Absorption cross section of specific absorption line (array, cm^2/particle)
    """
    x = 1 / wavelength - 1 / line_wavelength
    sigma_v = 9.12 * np.sqrt(T * amu / (1e4 * absorber_mass)) * 1e5 / (c * line_wavelength)
    #Doppler broadening, converted to wavenumber
    argument = (x + 1j * line_HWHM / c) / (sigma_v * np.sqrt(2)) 
    phi_line = np.real(wofz(argument)) / (sigma_v * np.sqrt(2 * np.pi))
    return phi_line*line_intensity



def absorption_cross_section(wavelength, T, chi):
    """Compute the total absorption cross section from all absorbing lines
    Input: Wavelength (array, cm), Temperature (scalar, K), Mixing ratio (scalar)
    Output: Total absorption cross section (array, cm^2/particle)
    """
    sigma = 0
    for value in param['Lines'].values():

        line_intensity = value[0]
        line_wavelength = value[1]
        line_HWHM = value[2]
        line_mass = value[3]

        sigma = sigma + line_absorption(wavelength,T,line_mass,line_wavelength,
                                        line_intensity,line_HWHM) * chi

    return sigma




"""
Calculation of the theoretical transit depth in two steps,
first calculating the (wavelength-dependent) optical depth along the 
chord at various impact parameters, then integrate over all possible 
impact parameters to obtain the (wavelength-dependent) transit depth
"""

def optical_depth(z, wavelength, yz_grid = None):   
    """Calculate the optical depth for all absorbers along the chord 
    """

    x = np.linspace(0, R_s, int(x_steps) + 1)[:-1] + 0.5 * R_s / float(x_steps)
   
    sigma = absorption_cross_section(wavelength, T, chi)
    
    xx, zz = np.meshgrid(x, z)

    r = np.sqrt(xx**2 + zz**2)

    n = number_density(r)
        
    delta_x = R_s / float(x_steps)

    N = 2 * delta_x * np.sum(n, axis = 1)

    tau = np.tensordot(sigma, N, axes = 0)

    return tau


def transit_depth(wavelength):
    """Calculate the wavelength-dependent transit depth
    """

    
    z = np.linspace(R_0, R_s, int(z_steps) + 1)[:-1] + 0.5 * (R_s - R_0) / float(z_steps)

    single_chord = np.exp(-optical_depth(z, wavelength))
        
    delta_z = (R_s - R_0) / float(z_steps)

    sum_over_chords = delta_z * np.tensordot(z, single_chord, axes = [0, 1])
    
    return 2. / (R_s**2 - R_0**2) * sum_over_chords

""".
Read in values from the setup file
"""

with open('../settings.txt') as file:
    param = json.load(file)



R_s = param['R_star']
R_0 = param['R_0']
M_p = param['M_p']

wavelength = np.linspace(param['lower_w'], param['upper_w'], 1 + int((param['upper_w'] - param['lower_w']) / param['resolution']))
x_steps = param['x_steps']
z_steps = param['z_steps']

number_density_dict = {'barometric': barometric, 'hydrostatic': hydrostatic, 'exomoon': exomoon}

for key_scenario in param['Scenarios'].keys():
    number_density = number_density_dict[key_scenario]
    T = param['Scenarios'][key_scenario][0]
    P_0 = param['Scenarios'][key_scenario][1]
    mu = param['Scenarios'][key_scenario][2]

    chi = param['Species']['sodium'][0][1]


"""
Execute the code and save the output as a .txt file
"""

H = k_B*T*R_0**2/(G*mu*M_p)
g = G*M_p/R_0**2
kappa = absorption_cross_section(wavelength, T, chi) / mu
R = R_0 + H * (euler_mascheroni + np.log(P_0 * kappa / g * np.sqrt(2 * np.pi * R_0 / H)))
benchmark_spectrum = 1 - (R**2 - R_0**2) / (R_s**2 - R_0**2)



spectrum = transit_depth(wavelength)

np.savetxt('../test.txt', np.c_[wavelength, spectrum, benchmark_spectrum])





elapsed_time = datetime.now() - startTime

print("PROMETHEUS finished, yay! Elapsed time is:")
print(elapsed_time)
print("The maximal flux decrease due to atmospheric/exospheric absorption in percent is:")
print(np.round(100*(1-min(spectrum)), 5))
print("The minimal flux decrease due to atmospheric/exospheric absorption in percent is:")
print(np.round(100*(1-max(spectrum)), 5))

print("The maximal BENCHMARK flux decrease due to atmospheric/exospheric absorption in percent is:")
print(np.round(100*(1-min(benchmark_spectrum)), 5))
print("The minimal BENCHMARK flux decrease due to atmospheric/exospheric absorption in percent is:")
print(np.round(100*(1-max(benchmark_spectrum)), 5))