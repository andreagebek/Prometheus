"""
File which stores natural constants in cgs units.
Created on 2. June 2021 by Andrea Gebek.
"""

import numpy as np

e = 4.803e-10   # Elementary charge
m_e = 9.109e-28
c = 2.998e10
G = 6.674*10**(-8)
k_B = 1.381*10**(-16)
amu = 1.661*10**(-24)
R_J = 6.99e9        #Jupiter radius
M_J = 1.898e30      #Jupiter mass
M_E = 5.974e27      #Earth mass
R_sun = 6.96e10     # Solar radius
M_sun = 1.988e33
R_Io = 1.822e8      # Io radius
euler_mascheroni = 0.57721 
AU = 1.496e13   # Conversion of one astronomical unit into cm

"""
Masses for the absorbing atoms/ions (this is not a really elegant way,
but I couldn't find another solution yet)
"""

speciesInfoDict = {'NaI': ['Na', '1', 22.99 * amu], 'KI': ['K', '1', 39.0983 * amu], 'SiI': ['Si', '1', 28.0855 * amu], 'SiII': ['Si', '2', 28.0855 * amu],
'SiIII': ['Si', '3', 28.0855 * amu], 'SiIV': ['Si', '4', 28.0855 * amu], 'MgI': ['Mg', '1', 24.305 * amu], 'MgII': ['Mg', '2', 24.305 * amu]}

"""
Planetary parameters
Format: [Stellar radius (cm), Stellar mass (g), Reference radius (cm), Planetary mass(g), Orbital distance (cm),
Stellar effective temperature (K), Stellar surface gravity (log10(cm/s^2)), Metallicity [Fe/H], Alpha-enhancement [alpha/Fe]]
WASP-49b: Wyttenbach et al. 2017 (Metallicity from Sousa+ 2018, Alpha-enhancement unknown)
HD189733b: Wyttenbach et al. 2015 (T_eff, log_g, and Metallicity from Chavero+ 2019, Alpha-enhancement unknown)
55Cancri-e: Crida et al. 2018 (a_p, T_eff, log_g, and Metallcity from Bourrier+ 2018, Alpha-enhancement unknown) 
"""

planetsDict = {'WASP-49b': [1.038 * R_sun, 1.003 * M_sun, 1.198 * R_J, 0.399 * M_J, 0.03873 * AU, 5600, 4.5, -0.08, 0],
'HD189733b': [0.756 * R_sun, 0.823 * M_sun, 1.138 * R_J, 1.138 * M_J, 0.0312 * AU, 5201, 4.64, -0.02, 0],
'55Cancri-e': [0.98 * R_sun, 1.015 * M_sun, 0.1737 * R_J, 0.02703 * M_J, 0.01544 * AU, 5172, 4.43, 0.35, 0]}


