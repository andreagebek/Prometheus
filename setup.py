"""
Run a question and answer session and store a settings.txt 
file to run prometheus.
Created on 2. June 2021 by Andrea Gebek.
"""

import numpy as np
import sys
import json
from constants import *


def read_value(text, lower, upper, unit, round = True):
    
    while True:
        
        if round:
            string = input(text + ' [' + '{:.0e}'.format(lower) + ', ' + '{:.0e}'.format(upper) + '] ').replace(' ', '')
        
        else:
            string = input(text + ' [' + '{:e}'.format(lower) + ', ' + '{:e}'.format(upper) + '] ').replace(' ', '')

        if string == '':
            print('You actually do have to enter something you scumbag.')
            continue
        
        value = float(string)

        if value > lower and value < upper:
            return value * unit
            break

        else:
            print('The value you entered is not in the appropriate interval.')



def read_str(text, options = 0):

    while True:

        if options == 0:
            string = input(text).replace(' ', '')

        else:
            string = input(text + ' ' + str(options).replace("'",'') + ' ').replace(' ', '')

        if string == '':
            print('You actually do have to enter something you scumbag.')
        
        elif options != 0 and string not in options:
            print('Please select a valid option.')

        else:
            return string
            break


"""
Fundamentals
"""


system_list = []
for key, value in planets_dict.items():
    system_list.append(key)

system_list.append('0')

systemname = read_str('Enter the name of the exoplanetary system or zero if parameters are entered manually:', system_list)

if systemname == '0':
    R_star = read_value('Enter the radius of the host star in solar radii:', 1e-5, 1e5, R_sun)
    R_0 = read_value('Enter the radius of the exoplanet in Jupiter radii:', 1e-5, 1e5, R_J)
    M_p = read_value('Enter the mass of the exoplanet in Jupiter masses:', 1e-5, 1e3, M_J)

else:
    R_star = planets_dict[systemname][0]
    R_0 = planets_dict[systemname][1]
    M_p = planets_dict[systemname][2]



direction = read_str('Do you want to perform forward or inverse modelling?', ['forward'])

mode = read_str('Do you want to model a transmission spectrum or a light curve?', ['spectrum'])

dishoom_import = read_str('Do you want to import parameters from DISHOOM?', ['no'])



"""
Scenarios for the spatial distribution of the medium
"""

scenario_dict = {}
scenario_list.append('0')

while True:
    scenario_name = read_str('Enter the name of the scenario for the number density profile or \
0 to stop adding absorption sources:', scenario_list)

    if scenario_name == '0':
        break


    elif scenario_name == 'barometric' or scenario_name == 'hydrostatic':
        winds = read_str('Add winds to the model?', ['no'])
        T = read_value('Enter the temperature of the atmosphere in Kelvin:', 1, 1e6, 1)
        P_0 = read_value('Enter the pressure at the reference radius in bar:', 1e-12, 1e6, 1e6)
        mu = read_value('Enter the mean molecular weight of the atmosphere in atomic mass units:', 0.05, 500, amu)

        params = [T, P_0, mu]

    elif scenario_name == 'exomoon':
        vel_dist = read_str('Enter the velocity distribution of the atoms:', ['maxwell_boltzmann'])
        R_moon = read_value('Enter the radius of the moon in Io radii:', 1e-3, 1e3, R_Io)

        params = [R_moon]

    
    scenario_dict[scenario_name] = params

    scenario_list.remove(scenario_name)

if len(scenario_dict) == 0:
    print('You have not added any absorption sources you scrub.')
    sys.exit()


"""
Absorbers
"""


while True:
    absorber_mode = read_str('How do you want to specify the absorption cross section?', ['BuiltinLine']) #BuiltinCont, ImportLineList

    if absorber_mode == 'BuiltinLine':
        
        absorptionlines_list = []
        for key, value in absorptionlines_dict.items():
            absorptionlines_list.append(key)
        absorptionlines_list.append('0')
        lines_dict = {}        
        
        while True:
            line = read_str('Enter the name of the line you want to consider, or 0 to move on:', absorptionlines_list)

            if line == '0':
                break

            else:
                lines_dict[line] = absorptionlines_dict[line]
                absorptionlines_list.remove(line)

        if len(lines_dict) == 0:
            print('You have not added any absorption lines you scrub.')
            sys.exit()

        species_dict = {}
        for value in lines_dict.values():
            element = value[4]
            species_dict[element] = 0
        
        for key_species in species_dict.keys():
            params = []
            for key_scenario in scenario_dict.keys():

                if key_scenario == 'barometric' or key_scenario == 'hydrostatic':
                    chi = read_value('Enter the mixing ratio of ' + key_species + ' in the ' + key_scenario + ' scenario:', 
                    1e-15, 1, 1)
                    params.append([key_scenario, chi])

                elif key_scenario == 'exomoon':
                    N = read_value('Enter the number of absorbing ' + key_species + ' atoms in the ' + key_scenario + ' scenario:', 
                    1e10, 1e50, 1)
                    v_mean = read_value('Enter the mean (thermal) velocity of ' + key_species + ' in the ' + key_scenario + '\
 scenario in km/s:', 1e-3, 1e5, 1e5)
                    params.append([key_scenario, N, v_mean])                   

            species_dict[key_species] = params
                    

    break


"""
Performance parameters
"""

if mode == 'spectrum':
    lower_w = read_value('Enter the lower wavelength border in Angstrom:', 1e-3, 1e12, 1e-8)
    upper_w = read_value('Enter the upper wavelength border in Angstrom:', lower_w * 1e8, 1e12, 1e-8,
    round = False)
    resolution = read_value('Eneter the resolution of the wavelength grid in Angstrom:', 1e-6,
    (upper_w - lower_w) * 1e8 / 2., 1e-8, round = False)

x_steps = read_value('Enter the steps for the spatial discretization along the chord:', 2, 1e6, 1)
z_steps = read_value('Enter the steps for the spatial discretization in z-direction:', 2, 1e6, 1)


"""
Write parameter dictionary and store it as json file
"""


parameters = {'R_star': R_star, 'R_0': R_0, 'M_p': M_p, 'direction': direction, 'mode': mode, 'dishoom_import': dishoom_import,
'Scenarios': scenario_dict, 'Lines': lines_dict, 'Species': species_dict,
'lower_w': lower_w, 'upper_w': upper_w, 'resolution': resolution, 'x_steps': x_steps, 'z_steps': z_steps}


with open('../settings.txt', 'w') as outfile:
    json.dump(parameters, outfile)

