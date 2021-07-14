"""
Run a question and answer session and store a settings.txt 
file to run prometheus.
Created on 2. June 2021 by Andrea Gebek.
"""

import numpy as np
import sys
import json
from constants import *


def read_value(text, lower, upper, unit, round = True, digits = 0, accept_borders = False):
    
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



"""
Fundamentals
"""
print('\nWelcome to PROMETHEUS! First, read in some parameters related to the architecture of the exoplanetary system.\n')

system_list = []
for key, value in planets_dict.items():
    system_list.append(key)

system_list.append('0')

systemname = read_str('Enter the name of the exoplanetary system or zero if parameters are entered manually:', system_list)

if systemname == '0':
    R_star = read_value('Enter the radius of the host star in solar radii:', 1e-5, 1e5, R_sun)
    R_0 = read_value('Enter the radius of the exoplanet in Jupiter radii:', 1e-5, 1e5, R_J)
    M_p = read_value('Enter the mass of the exoplanet in Jupiter masses:', 1e-5, 1e3, M_J)
    a_p = read_value('Enter the orbital distance between planet and star in AU:', 1e-5, 1e3, AU)

else:
    R_star = planets_dict[systemname][0]
    R_0 = planets_dict[systemname][1]
    M_p = planets_dict[systemname][2]
    a_p = planets_dict[systemname][3]

architecture_dict = {'R_star': R_star, 'R_0': R_0, 'M_p': M_p, 'a_p': a_p}

"""
direction = read_str('Do you want to perform forward or inverse modelling?', ['forward'])

mode = read_str('Do you want to model a transmission spectrum or a light curve?', ['spectrum'])

dishoom_import = read_str('Do you want to import parameters from DISHOOM?', ['no'])
"""


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
            
        #winds_ON = read_str('Do you want to add winds in this scenario?', ['no'])

        params = {'T': T, 'P_0': P_0, 'mu': mu, 'RayleighScatt': RayleighScatt}

        if scenario_name == 'barometric':
            scenario_list.remove('hydrostatic')
        else:
            scenario_list.remove('barometric')

        PlanetarySource = True

    if scenario_name == 'escaping':

        q_esc = read_value('Enter the power law index for the escaping atmosphere:', 3, 20, 1)
        norm_esc = read_str('Do you want to normalize the number density profile via (total) pressure or via number \
of absorbing atoms at the base of the wind?', ['pressure', 'number'])

        if norm_esc == 'pressure':

            T = read_value('Enter the temperature of the escaping wind in Kelvin:', 1, 1e6, 1)
            P_0 = read_value('Enter the pressure at the base of the wind in bar:', 1e-15, 1e3, 1e6)

            params = {'q_esc': q_esc, 'T': T, 'P_0': P_0}
        
        else:

            params = {'q_esc': q_esc}

        PlanetarySource = True

    if scenario_name == 'exomoon':

        R_moon = read_value('Enter the radius of the moon in Io radii:', 1e-3, 1e3, R_Io)

        params = {'R_moon': R_moon}
    
        ExomoonSource = True

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

    if not PlanetarySource:
        
        ExomoonOffCenter = not(read_str('Do you want to make the approximation that the exomoon sits at the planetary center, neglecting the planet itself?', ['yes', 'no']))
    
    else:

        ExomoonOffCenter = True # Both a planetary source and an exomoon are present
    
if ExomoonOffCenter:

    print('\nYou need to specify some additional parameters related to the system architecture of the exoplanet-exomoon system.\n')

    R_moon = scenario_dict['exomoon']['R_moon']

    a_moon = read_value('Enter the orbital distance between the exomoon and the planet in planetary radii (measured from the centers of the bodies):', 1, a_p / R_0, R_0, round = False)
    alpha_moon = read_value('Enter the phase angle of the exomoon in degrees (0 is to the right when the exoplanetary system is viewed from the observers \
perspective, 90 corresponds to the exomoon sitting between the planet and the host star):', 0, 360, 2. * np.pi / 360., round = False, digits = 2, accept_borders = True)
    z_moon = read_value('Enter the elevation (z-coordinate) of the exomoon with respect to the star-planet plane in planetary radii:', 0, a_moon / R_0, R_0, round = False, accept_borders = True)

    architecture_dict.update({'R_moon': R_moon, 'a_moon': a_moon, 'alpha_moon': alpha_moon, 'z_moon': z_moon})

    sphericalSymmetry = False


if len(scenario_dict) == 0:
    print('You have not added any absorption sources! Your loss. PROMETHEUS exits now.')
    sys.exit()



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

upper_x = read_value('Enter the half chord length (x-direction) for the numerical integration along the x-axis in planetary radii:', 0, a_p / R_0, R_0, round = False)
x_steps = read_value('Enter the steps for the spatial discretization along the chord (x-direction):', 2, 1e6, 1)

z_steps = read_value('Enter the steps for the spatial discretization in z-direction:', 2, 1e6, 1)

grids_dict = {'lower_w': lower_w, 'upper_w': upper_w, 'resolution': resolution, 'upper_x': upper_x, 'x_steps': x_steps, 'z_steps': z_steps}

if not sphericalSymmetry:
    phi_steps = read_value('Enter the steps for the spatial discretization for the polar coordinate (phi-direction):', 2, 1e4, 1)
    grids_dict['phi_steps'] = phi_steps

"""
Additional output
"""

print('\nFinally, specify the output.\n')

output_dict = {}

paramsFilename = read_str('How do you want to name the txt file containing the parameters of this session (enter the file name without the .txt ending)?')

if 'barometric' in scenario_dict.keys():
    benchmark = read_str('Do you want to record the analytical benchmark for the barometric scenario?', ['yes', 'no'])
    output_dict['benchmark'] = benchmark

else:
    output_dict['benchmark'] = False

record_tau = read_str('Do you want to record the optical depth for all chords of the spatial grid at the wavelength with the largest flux decrease?', ['yes', 'no'])

output_dict['record_tau'] = record_tau

"""
Write parameter dictionary and store it as json file
"""

print('All parameters are stored! To run PROMETHEUS, type <python main.py filename> and replace <filename> with the name you specified for the parameter txt file.')

parameters = {'Architecture': architecture_dict, 'Scenarios': scenario_dict, 'Lines': lines_dict, 'Species': species_dict, 'Grids': grids_dict, 'Output': output_dict,
'sphericalSymmetry': sphericalSymmetry, 'ExomoonSource': ExomoonSource, 'ExomoonOffCenter': ExomoonOffCenter}


with open('../' + paramsFilename + '.txt', 'w') as outfile:
    json.dump(parameters, outfile)