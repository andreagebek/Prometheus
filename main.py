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

""".
Read in values from the setup file
"""

with open('../settings.txt') as file:
    parameters = json.load(file)

print(parameters)
print(parameters['Scenarios'].keys())
stop

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


def hydrostatic_basic(R):
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
    x = 1/wavelength-1/line_wavelength
    sigma_v = 9.12*np.sqrt(T*amu/(1e4*absorber_mass))*1e5*1/(c*line_wavelength)
    #Doppler broadening, converted to wavenumber
    argument = (x+1j*line_HWHM/c)/(sigma_v*np.sqrt(2)) 
    phi_line = np.real(wofz(argument))/(sigma_v*np.sqrt(2*pi))
    return phi_line*line_intensity



def absorption_cross_section(wavelength, T, chi):
    """Compute the total absorption cross section from all absorbing lines
    Input: Wavelength (array, cm), Temperature (scalar, K), Mixing ratio (scalar)
    Output: Total absorption cross section (array, cm^2/particle)
    """
    sigma = 0
    for line in LineList:
        
        if type(absorption_lines) == list:
            line_wavelength = line[0]
            line_intensity = line[1]
            line_HWHM = line[2]
            
        else:
            line_wavelength = 1. / line[3]
            line_intensity = line[4]
            line_HWHM = line[5] / 2.  #line[5] is the Einstein A coefficient

            
        sigma = sigma + line_absorption(wavelength,T,absorber_mass,line_wavelength,
                                        line_intensity,line_HWHM)*chi

    return sigma




"""
Calculation of the theoretical transit depth in two steps,
first calculating the (wavelength-dependent) optical depth along the 
chord at various impact parameters, then integrate over all possible 
impact parameters to obtain the (wavelength-dependent) transit depth
"""

def optical_depth(z, wavelength, yz_grid = None):   
    """Calculate the optical depth for all absorbers along the chord 
    Input: Grid in z-direction (array, cm), Wavelength (array, cm)
    Output: Optical depth (2d-array)
    """
    x = np.linspace(0,R_s,int(x_steps))
    delta_x_i = np.array([R_s/(2*float(x_steps-1))])
    delta_x_middle = np.linspace(R_s/float(x_steps-1),
                                 R_s/float(x_steps-1),
                                  int(x_steps-2))
    delta_x = np.concatenate((delta_x_i, delta_x_middle, delta_x_i))
    x_squared = x**2
    x_ones = np.ones(int(x_steps))
    


    
    sigma = absorption_cross_section(wavelength, T, chi)
    
    if number_density == torus: 
        yz_ones = np.ones((int(y_steps), int(z_steps)))
        x_matrix = np.tensordot(x_squared, yz_ones, axes = 0)
        y_matrix = np.tensordot(x_ones, np.square(yz_grid), axes = 0)
        
        a = np.sqrt(x_matrix + y_matrix)
        
        n = number_density(a, z) # Function of x, y and z


    else:
        z_ones = np.ones(int(len(z)))
        z_squared = z**2  
        x_matrix = np.tensordot(x_squared, z_ones, axes = 0)
        z_matrix = np.tensordot(x_ones, z_squared, axes = 0)
        r = np.sqrt(x_matrix + z_matrix)    
    
        n = number_density(r)
        
        
    N = 2*np.tensordot(delta_x, n, axes = 1)
    tau = np.tensordot(sigma, N, axes = 0)
    
    return tau


def transit_depth(wavelength):
    """Calculate the wavelength-dependent transit depth
    Input: Wavelength (array, cm)
    Output: Transit depth (array)
    """
    
    if number_density == torus:
        starting_z = 1
    else:
        starting_z = R_0
    
    z = np.logspace(np.log10(starting_z), np.log10(R_s), int(z_steps))
    
    q_z = (R_s / starting_z)**(1. / (float(z_steps) - 1.))
    delta_z_i = np.array([starting_z * (q_z - 1.) / 2.])
    delta_z_middle = np.logspace(np.log10(starting_z * (q_z - 1.) * (q_z + 1.) / 2.),
                                 np.log10(starting_z * q_z**(float(z_steps) - 3.) * (q_z + 1.) * (q_z - 1.) / 2.),
                                 int(z_steps - 2.))
    delta_z_f = np.array([starting_z * q_z**(float(z_steps) - 2.) * (q_z - 1.) / 2.])
    delta_z = np.concatenate((delta_z_i, delta_z_middle, delta_z_f))
    
    
    if number_density == torus:
        delta_y = np.array([])
        
        for i in range(len(z)):
            lowest_y = np.real(np.sqrt(R_0**2 - z[i]**2 + 0j))
            highest_y = np.real(np.sqrt(R_s**2 - z[i]**2 + 0j))
            delta_y = np.append(delta_y, (highest_y - lowest_y) / float(y_steps))
            
            if i == 0:
                yz_grid = np.linspace(lowest_y, highest_y, int(y_steps))
                
            else:
                y_grid = np.linspace(lowest_y, highest_y, int(y_steps))
                yz_grid = np.vstack((yz_grid, y_grid))

        single_chord = np.exp(-optical_depth(z, wavelength, yz_grid.T))
        
        integral_over_y = np.sum(single_chord, axis = 1) * delta_y

        sum_over_chords = np.tensordot(integral_over_y, delta_z, axes = (1,0)) * 2. / pi


    
    else:
        z_delta_z = np.multiply(z, delta_z)
        
        single_chord = np.exp(-optical_depth(z, wavelength))
        
        sum_over_chords = np.tensordot(z_delta_z, single_chord, axes = [0, 1])
    
    return 2. / (R_s**2 - R_0**2) * sum_over_chords


"""
Read in the configuration file (config_v3.py) and save the output .txt file
"""



planet_dict = {"WASP-49b": wasp49,
               "WASP-52b": wasp52,
               "WASP-69b": wasp69,
               "WASP-76b": wasp76,
               "WASP-121b": wasp121,
               "HD189733b": HD189,
               "HD209733b": HD209,
               "55 Cancri e": ffcnce,
               "TRAPPIST-1b": T1b,
               "TRAPPIST-1c": T1c}
               

planet_parameters = planet_dict[planet_name]

R_0 = planet_parameters[0]
M_p = planet_parameters[1]
R_s = planet_parameters[2]

if type(absorption_lines) == list:
    LineList = absorption_lines
    absorber_mass = LineList[0][3]
        
else:
    LineList = np.loadtxt(absorption_lines)
    

number_density_dict = {"Barometric": barometric,
               "Hydrostatic-Basic": hydrostatic_basic,
               "Escaping": escaping,
               "Exomoon": exomoon,
               "Torus": torus}

number_density = number_density_dict[number_density_name]

if number_density_name == "Escaping" or number_density_name == "Torus":
    T = v_mean**2 *absorber_mass * pi/ ( 8 * k_B)


if number_density_name == "Exomoon":
    T = v_mean**2 * absorber_mass * pi / (8 * k_B)
    R_0 = R_moon


wavelength = np.arange(lower_wavelength, upper_wavelength+resolution, resolution)

spectrum = transit_depth(wavelength)

np.savetxt(OUTPUT_PATH + name_txt_file, np.c_[wavelength, spectrum])

elapsed_time = datetime.now() - startTime

print("PROMETHEUS finished, yay! Elapsed time is:")
print(elapsed_time)
print("The maximal flux decrease due to atmospheric/exospheric absorption in percent is:")
print(np.round(100*(1-min(spectrum)), 5))
print("The minimal flux decrease due to atmospheric/exospheric absorption in percent is:")
print(np.round(100*(1-max(spectrum)), 5))
