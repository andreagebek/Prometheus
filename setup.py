"""
Run a question and answer session and store a settings.txt 
file to run prometheus.
Created on 2. June 2021 by Andrea Gebek.
"""

import numpy as np
import sys
import json
sys.path.insert(0, 'pythonScripts') # Import from pythonScripts folder
from constants import *


def read_value(text, lower, upper, unit, round = True, digits = 0, accept_borders = False, break_zero = False):
    
    while True:

        if accept_borders:
            brackl = ' ['
            brackr = '] '
        else:
            brackl = ' ('
            brackr = ') '
        
        if round:
            string = input(text + brackl + '{:.0e}'.format(lower) + ', ' + '{:.0e}'.format(upper) + brackr).replace(' ', '')

        elif digits > 0:
            string = input(text + brackl + '{:.{n}e}'.format(lower, n = digits) + ', ' + '{:.{n}e}'.format(upper, n = digits) + brackr).replace(' ', '')

        else:
            string = input(text + brackl + '{:e}'.format(lower) + ', ' + '{:e}'.format(upper) + brackr).replace(' ', '')

        if string == '':
            print('You actually do have to enter something!')
            continue
        
        value = float(string)

        if break_zero and value == 0:
            return 0
            break

        if value > lower and value < upper:
            return value * unit
            break

        elif value >= lower and value <= upper and accept_borders:
            return value * unit
            break

        else:
            print('The value you entered is not in the appropriate interval.')



def read_str(text, options = 0):

    while True:

        if options == 0:
            string = input(text + ' ').replace(' ', '')

        else:
            string = input(text + ' ' + str(options).replace("'",'') + ' ').replace(' ', '')

        if string == '':
            print('You actually do have to enter something!')
        
        elif options != 0 and string not in options:
            print('Please select a valid option.')

        else:
            if string == 'yes':
                return True
            elif string == 'no':
                return False
            else:
                return string
            break


"""
Overall properties
"""

PlanetarySource = False
ExomoonSource = False
ExomoonOffCenter = False
sphericalSymmetry = True
CLV_variations = False
RM_effect = False
DopplerOrbitalMotion = False
DopplerPlanetRotation = False
RadialWinds = False



"""
Fundamentals
"""
print('\nWelcome to PROMETHEUS! First, read in some parameters related to the architecture of the exoplanetary system.\n')

mode = read_str('Do you want to model a transmission spectrum or a light curve?', ['spectrum', 'lightcurve'])

if mode == 'lightcurve':
    sphericalSymmetry = False

    CLV_variations = read_str('Do you want to take center-to-limb variations into account?', ['yes', 'no'])
    RM_effect = read_str('Do you want to take the Rossiter-McLaughlin-Effect into account (note that this means that you \
have to provide additional information about the host star to specifiy its spectrum)?', ['yes', 'no'])
    DopplerOrbitalMotion = read_str('Do you want to consider the Doppler shifts due to planetary/exomoon orbital motion?', ['yes', 'no'])


system_list = []
for key, value in planets_dict.items():
    system_list.append(key)

system_list.append('0')

systemname = read_str('Enter the name of the exoplanetary system or zero if parameters are entered manually:', system_list)

if systemname == '0':
    R_star = read_value('Enter the radius of the host star in solar radii:', 1e-5, 1e5, R_sun)
    R_0 = read_value('Enter the radius of the exoplanet in Jupiter radii:', 1e-5, 1e5, R_J)
    M_p = read_value('Enter the mass of the exoplanet in Jupiter masses:', 1e-5, 1e3, M_J)

    if mode == 'lightcurve':
        a_p = read_value('Enter the orbital distance between planet and star in AU:', 1e-5, 1e3, AU)
        M_star = read_value('Enter the mass of the host star in solar masses:', 1e-5, 1e10, M_sun)

    if RM_effect:
        T_eff = read_value('Enter the effective temperature of the star in Kelvin:', 2300, 12000, 1, round = False, accept_borders = True)
        log_g = read_value('Enter the logarithmic value of the surface gravity in log10(cm/s^2):', 0, 6, 1, round = False, accept_borders = True)
        Fe_H = read_value('Enter the metallicity of the star [Fe/H]:', -4, 1, 1, round = False, accept_borders = True)

        if Fe_H > -3.5 and Fe_H < 0.25: # The PHOENIX library has varying alpha-enhancement only for this metallicity range
            alpha_Fe = read_value('Enter the alpha-enhancement of the star [alpha/Fe]:', -0.2, 1.2, 1, round = False, accept_borders = True)
        else:
            alpha_Fe = 0

else:
    R_star = planets_dict[systemname][0]
    M_star = planets_dict[systemname][1]
    R_0 = planets_dict[systemname][2]
    M_p = planets_dict[systemname][3]
    a_p = planets_dict[systemname][4]
    T_eff = planets_dict[systemname][5]
    log_g = planets_dict[systemname][6]
    Fe_H = planets_dict[systemname][7]
    alpha_Fe = planets_dict[systemname][8]

architecture_dict = {'R_star': R_star, 'R_0': R_0, 'M_p': M_p}

if mode == 'lightcurve':
    architecture_dict.update({'M_star': M_star, 'a_p': a_p})

if RM_effect:
    period_starrot = read_value('Enter the period of the stellar rotation in days:', 0, 1000, 86400)
    inclination_starrot = read_value('Enter the inclination of the stellar rotation against the planetary orbital plane in degrees:', 0, 180, np.pi / 180., round = False, accept_borders = True)
    azimuth_starrot = read_value('Enter the angle between the line of sight from the observer to the system and the angular momentum vector of the stellar rotation, \
projected onto the orbital plane in degrees', 0, 360, np.pi / 180., round = False, accept_borders = True)

    architecture_dict.update({'T_eff': T_eff, 'log_g': log_g, 'Fe_H': Fe_H, 'alpha_Fe': alpha_Fe, 'period_starrot': period_starrot, 'inclination_starrot': inclination_starrot, 'azimuth_starrot': azimuth_starrot})

if CLV_variations:
    u1 = read_value('Enter the first (linear) coefficient for limb darkening:', -1, 1, 1, accept_borders = True)
    u2 = read_value('Enter the second (quadratic) coefficient for limb darkening:', -1, 1, 1, accept_borders = True)
    architecture_dict.update({'u1': u1, 'u2': u2})


#direction = read_str('Do you want to perform forward or inverse modelling?', ['forward'])


#dishoom_import = read_str('Do you want to import parameters from DISHOOM?', ['no'])


"""
Scenarios for the spatial distribution of the medium
"""

print('\nNow, specifiy the spatial distribution of the absorbing medium, i.e. the structure of the atmosphere/exosphere.\n')

scenario_dict = {}
scenario_list.append('0')

while True:
    scenario_name = read_str('Enter the name of the scenario for the number density profile or \
0 to stop adding absorption sources:', scenario_list)

    if scenario_name == '0':
        break


    if scenario_name == 'barometric' or scenario_name == 'hydrostatic':

        T = read_value('Enter the temperature of the atmosphere in Kelvin:', 1, 1e6, 1)
        P_0 = read_value('Enter the pressure at the reference radius in bar:', 1e-12, 1e6, 1e6)
        mu = read_value('Enter the mean molecular weight of the atmosphere in atomic mass units:', 0.05, 500, amu)
        #press_broad_ON = read_str('Do you want to add pressure broadening in this scenario?', ['no']) NOT HERE!
        RayleighScatt = read_str('Do you want to add Rayleigh scattering in this scenario?', ['yes', 'no'])

        params = {'T': T, 'P_0': P_0, 'mu': mu, 'M_p': M_p, 'R_0': R_0, 'RayleighScatt': RayleighScatt}

        if scenario_name == 'barometric':
            scenario_list.remove('hydrostatic')
        else:
            scenario_list.remove('barometric')

        PlanetarySource = True

    elif scenario_name == 'escaping':

        q_esc = read_value('Enter the power law index for the escaping atmosphere:', 3, 20, 1)
        norm_esc = read_str('Do you want to normalize the number density profile via (total) pressure or via number \
of absorbing atoms at the base of the wind?', ['pressure', 'number'])

        if norm_esc == 'pressure':

            P_0 = read_value('Enter the pressure at the base of the wind in bar:', 1e-15, 1e3, 1e6)
            T = read_value('Enter the temperature of the escaping wind in Kelvin:', 1, 1e6, 1)


            params = {'q_esc': q_esc, 'P_0': P_0, 'T': T, 'R_0': R_0}
        
        else:

            params = {'q_esc': q_esc, 'R_0': R_0}

        RadialWind = read_str('Do you want to add radially escaping winds in the escaping scenario?', ['yes', 'no'])
        params['RadialWind'] = RadialWind

        if RadialWind:

            vRadial_0 = read_value('Enter the velocity of the radially escaping wind at the reference radius in km/s:', 1e-3, 1e3, 1e5)
            params['vRadial_0'] = vRadial_0
            RadialWinds = True

        PlanetarySource = True

    elif scenario_name == 'exomoon':

        R_moon = read_value('Enter the radius of the moon in Io radii:', 1e-3, 1e3, R_Io)

        RadialWind = read_str('Do you want to add radially escaping winds in the escaping scenario?', ['yes', 'no'])

        params = {'R_moon': R_moon, 'RadialWind': RadialWind}

        if RadialWind:

            vRadial_0 = read_value('Enter the velocity of the radially escaping wind at the exomoon radius in km/s:', 1e-3, 1e3, 1e5)   
            params['vRadial_0'] = vRadial_0
            RadialWinds = True

        ExomoonSource = True


    elif scenario_name == 'torus':

        a_torus = read_value('Enter the distance between the center of the torus and the center of the exoplanet in planetary radii:', 1, 1000, R_0)
        v_ej = read_value('Enter the ejection velocity (which sets the torus scale height) of the particles from the torus in km/s:', 1e-2, 1e3, 1e5)
 
        params = {'a_torus': a_torus, 'v_ej': v_ej, 'M_p': M_p}

        sphericalSymmetry = False

    if 'T' in params.keys():

        thermBroad = read_str('Do you want to add thermal line broadening in this scenario?', ['yes', 'no'])

    else:

        thermBroad = read_str('Do you want to add pseudo-thermal line broadening via mean velocity of the absorber in this scenario?', ['yes', 'no'])

    params['thermBroad'] =  thermBroad


    scenario_dict[scenario_name] = params

    scenario_list.remove(scenario_name)



"""
If there is an exomoon, it is possible to take both planet and exomoon into account, requiring additional 
system architecture parameters for the exomoon.
"""



if ExomoonSource:

    architecture_dict['R_moon'] = R_moon

    if not PlanetarySource and sphericalSymmetry:
        
        ExomoonOffCenter = not(read_str('Do you want to make the approximation that the exomoon sits at the planetary center, neglecting the planet itself?', ['yes', 'no']))
    
    else:

        ExomoonOffCenter = True # Both a planetary source and an exomoon are present, or we are considering a light curve
    
if ExomoonOffCenter:

    print('\nYou need to specify some additional parameters related to the system architecture of the exoplanet-exomoon system.\n')

    R_moon = scenario_dict['exomoon']['R_moon']

    a_moon = read_value('Enter the orbital distance between the exomoon and the planet in planetary radii (measured from the centers of the bodies):', 1, a_p / R_0, R_0, round = False)
    starting_orbphase_moon = read_value('Enter the orbital phase of the moon when the planet is transiting. A moon orbital phase of 0 corresponds to the moon sitting \
between the planet and the observer, 0.25 means that the exomoon is located to the right of the planet when viewed from the observer.', -0.5, 0.5, 2. * np.pi, accept_borders = True)

    architecture_dict.update({'R_moon': R_moon, 'a_moon': a_moon, 'starting_orbphase_moon': starting_orbphase_moon})

    if mode == 'spectrum':
        z_moon = read_value('Enter the elevation (z-coordinate) of the exomoon with respect to the planet orbital plane in planetary radii:', 0, a_moon / R_0, R_0, 
        round = False, accept_borders = True)

        architecture_dict['z_moon'] = z_moon
    """
    elif mode == 'lightcurve':
        i_moon = read_value('Enter the inclination of the moon orbital plane with respect to the planet orbital plane in degrees:', 0, 90, 2. * np.pi / 360., accept_borders = True)
        
        if i_moon > 0:
            iphase_moon = read_value('Enter the orbital phase of the moon when it ascends through the planet orbital plane:', 0)
    """

    sphericalSymmetry = False


if len(scenario_dict) == 0:
    print('You have not added any absorption sources! Your loss. PROMETHEUS exits now.')
    sys.exit()

"""
If there is a planetary source, specifiy if planetary rotation is taken into account.
"""

if PlanetarySource:

    DopplerPlanetRotation = read_str('Do you want to consider the Doppler shifts due to planetary rotation? This only applies to planetary sources \
(barometric, hydrostatic, escaping).', ['yes', 'no'])

    if DopplerPlanetRotation:

        period_planetrot = read_value('Enter the period for the planetary rotation in days:', 0, 1000, 86400)
        architecture_dict['period_planetrot'] = period_planetrot

        sphericalSymmetry = False

# Now, calculate the dimension required for the absorption cross section with all the available information

if mode == 'lightcurve':
    if DopplerPlanetRotation or RadialWinds:
        sigma_dim = 5   # sigma(lambda, x, phi, rho, t)
    elif DopplerOrbitalMotion:
        sigma_dim = 2   # sigma(lambda, t)
    else:
        sigma_dim = 1   # sigma(lambda)

else:
    if (RadialWinds and ExomoonOffCenter) or DopplerPlanetRotation:
        sigma_dim = 4   # sigma(lambda, x, phi, rho)
    elif RadialWinds:
        sigma_dim = 3   # sigma(lambda, x, rho)
    else:
        sigma_dim = 1   # sigma(lambda)

"""
Specify the absorption lines and species-related parameters.
"""

print('\nSpecify the absorption lines and the parameters related to the abundance and other properties of the absorbing species.\n')

lines_dict = {}
species_dict = {}


absorptionlines_list = []
for key, value in absorptionlines_dict.items():
    absorptionlines_list.append(key)
absorptionlines_list.append('0')


while True:
    line = read_str('Enter the name of the line you want to consider, or 0 to move on:', absorptionlines_list)

    if line == '0':
        break

    else:
        lines_dict[line] = absorptionlines_dict[line]
        absorptionlines_list.remove(line)

for key_scenario in scenario_dict.keys():

    species_dict[key_scenario] = {}

    for value in lines_dict.values():

        element = value[4]
        species_dict[key_scenario][element] = {}

for key_scenario in species_dict.keys():

    for key_species in species_dict[key_scenario].keys():

        if 'P_0' in scenario_dict[key_scenario].keys():

            chi = read_value('Enter the mixing ratio of ' + key_species + ' in the ' + key_scenario + ' scenario:', 
            1e-15, 1, 1)

            params = {'chi': chi}
                
        else:

            N_abs = read_value('Enter the number of absorbing ' + key_species + ' atoms in the ' + key_scenario + ' scenario:', 
            1e10, 1e50, 1)

            params = {'chi': N_abs}

        if scenario_dict[key_scenario]['thermBroad']:

            if 'T' in scenario_dict[key_scenario].keys():

                params['T_abs'] = scenario_dict[key_scenario]['T']
            
            else:

                v_mean = read_value('Enter the mean (thermal) velocity of ' + key_species + ' in the ' + key_scenario + ' scenario in km/s:', 1e-3, 1e5, 1e5)
                absorber_mass = speciesMass_dict[key_species]

                T_abs = v_mean**2 * absorber_mass * np.pi / (8. * k_B)

                params['T_abs'] = T_abs
        
        else:

            params['T_abs'] = 0

        species_dict[key_scenario][key_species] = params


"""
Performance parameters
"""

print('\nAlmost done! Specify the discretization parameters for the wavelength and spatial grids.\n')


lower_w = read_value('Enter the lower wavelength border in Angstrom:', 1e-3, 1e12, 1e-8)
upper_w = read_value('Enter the upper wavelength border in Angstrom:', lower_w * 1e8, 1e12, 1e-8, round = False)
resolution = read_value('Enter the resolution of the wavelength grid in Angstrom:', 1e-6, (upper_w - lower_w) * 1e8 / 2., 1e-8, round = False)

x_border = read_value('Enter the half chord length (x-direction) for the numerical integration along the x-axis in planetary radii:', 0, a_p / R_0, R_0, round = False)
x_steps = read_value('Enter the number of bins for the spatial discretization along the chord (x-direction):', 2, 1e6, 1)

z_steps = read_value('Enter the number of bins for the spatial discretization in z-direction:', 2, 1e6, 1)


grids_dict = {'lower_w': lower_w, 'upper_w': upper_w, 'resolution': resolution, 'x_border': x_border, 'x_steps': x_steps, 'z_steps': z_steps}

if not sphericalSymmetry:
    phi_steps = read_value('Enter the number of bins for the spatial discretization for the polar coordinate (phi-direction):', 2, 1e4, 1)
    grids_dict['phi_steps'] = phi_steps

if mode == 'lightcurve':
    orbphase_border = read_value('Enter the orbital phase at which the light curve calculation starts and stops:', 0, 0.5, 2 * np.pi, accept_borders = True)
    orbphase_steps = read_value('Enter the number of bins for the orbital phase discretization:', 2, 1e4, 1)

    grids_dict.update({'orbphase_border': orbphase_border, 'orbphase_steps': orbphase_steps})

"""
Additional output
"""

print('\nFinally, specify the output.\n')

output_dict = {}

paramsFilename = read_str('How do you want to name the txt file containing the parameters of this session (enter the file name without the .txt ending)?')

if 'barometric' in scenario_dict.keys() and mode == 'spectrum':
    benchmark = read_str('Do you want to record the analytical benchmark for the barometric scenario?', ['yes', 'no'])
    output_dict['benchmark'] = benchmark

else:
    output_dict['benchmark'] = False

if mode == 'spectrum':
    record_tau = read_str('Do you want to record the optical depth for all chords of the spatial grid at the wavelength with the largest flux decrease?', ['yes', 'no'])

elif mode == 'lightcurve':
    record_tau = False

output_dict['record_tau'] = record_tau

"""
Write parameter dictionary and store it as json file
"""

print('\n\nAll parameters are stored! To run PROMETHEUS, type <python main.py ' + paramsFilename + '>.\n')

parameters = {'Architecture': architecture_dict, 'Scenarios': scenario_dict, 'Lines': lines_dict, 'Species': species_dict, 'Grids': grids_dict, 'Output': output_dict,
'sphericalSymmetry': sphericalSymmetry, 'ExomoonSource': ExomoonSource, 'CLV_variations': CLV_variations, 'RM_effect': RM_effect, 
'DopplerOrbitalMotion': DopplerOrbitalMotion, 'DopplerPlanetRotation': DopplerPlanetRotation, 'RadialWinds': RadialWinds, 'sigma_dim': sigma_dim, 'mode': mode}


with open('../' + paramsFilename + '.txt', 'w') as outfile:
    json.dump(parameters, outfile)