"""
This file stores functions to calculate absorption cross sections.
Created on 15. July 2021 by Andrea Gebek.
"""

import numpy as np
from scipy.special import wofz
from constants import *

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



def absorption_cross_section(wavelength, chi, T, species, lines_dict):
    """Compute the total absorption cross section from all absorbing lines
    """

    sigma = 0

    for value in lines_dict.values():

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