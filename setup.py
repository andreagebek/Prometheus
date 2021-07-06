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
            string = input(text).replace(' ', '')

        else:
            string = input(text + ' ' + str(options).replace("'",'') + ' ').replace(' ', '')

        if string == '':
            print('You actually do have to enter something!')
        
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
    a_p = read_value('Enter the orbital distance between planet and star in AU:', 1e-5, 1e3, AU)

else:
    R_star = planets_dict[systemname][0]
    R_0 = planets_dict[systemname][1]
    M_p = planets_dict[systemname][2]
    a_p = planets_dict[systemname][3]

R_moon = 0
a_moon = 0
alpha_moon = 0
z_moon = 0


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


    if check('atmosphere', scenario_name):

        T = read_value('Enter the temperature of the atmosphere in Kelvin:', 1, 1e6, 1)
        P_0 = read_value('Enter the pressure at the reference radius in bar:', 1e-12, 1e6, 1e6)
        mu = read_value('Enter the mean molecular weight of the atmosphere in atomic mass units:', 0.05, 500, amu)
        #press_broad_ON = read_str('Do you want to add pressure broadening in this scenario?', ['no']) NOT HERE!
        rayleigh_scatt_ON = read_str('Do you want to add Rayleigh scattering in this scenario?', ['yes', 'no'])
            
        #winds_ON = read_str('Do you want to add winds in this scenario?', ['no'])

        params = [T, P_0, mu]

        if rayleigh_scatt_ON == 'yes':
            params.append('rayleigh_scatt')

    if scenario_name == 'escaping':

        q_esc = read_value('Enter the power law index for the escaping atmosphere:', 3, 20, 1)
        norm_esc = read_str('Do you want to normalize the number density profile via (total) pressure or via number \
of absorbing atoms at the base of the wind?', ['pressure', 'number'])

        if norm_esc == 'pressure':

            T = read_value('Enter the temperature of the escaping wind in Kelvin:', 1, 1e6, 1)
            P_0 = read_value('Enter the pressure at the base of the wind in bar:', 1e-15, 1e3, 1e6)

            params = [q_esc, norm_esc, T, P_0]
        
        else:

            params = [q_esc, norm_esc]

    if scenario_name == 'exomoon':

        R_moon = read_value('Enter the radius of the moon in Io radii:', 1e-3, 1e3, R_Io)

        params = [R_moon, 'center_exomoon_off']

    if check('evaporative', scenario_name, params):

        therm_broad_ON = read_str('Do you want to add pseudo-thermal broadening (via mean velocity of absorber) in this scenario?', ['yes', 'no'])

        if therm_broad_ON == 'yes':
            params.append('therm_broad')

    scenario_dict[scenario_name] = params

    scenario_list.remove(scenario_name)


if check('exomoon_center', scenario_dict.keys()):
    exomoon_center_ON = read_str('Do you want to make the approximation that the exomoon sits at the planetary center, neglecting the planet itself?', ['yes', 'no'])
    
    if exomoon_center_ON == 'yes':
        scenario_dict['exomoon'][1] = 'center_exomoon_on'
    
    else:
        a_moon = read_value('Enter the orbital distance between the exomoon and the planet in planetary radii:', 1, a_p / R_0, R_0, round = False)
        alpha_moon = read_value('Enter the phase angle of the exomoon in degrees (0 is to the left when the exoplanetary system is viewed from the observers \
perspective, 90 corresponds to the exomoon sitting between the planet and the observer):', 0, 360, 2. * np.pi / 360., round = False, digits = 2, accept_borders = True)
        z_moon = read_value('Enter the elevation (z-coordinate) of the exomoon with respect to the star-planet plane in planetary radii:', 0, a_moon / R_0, R_0, round = False, accept_borders = True)



if len(scenario_dict) == 0:
    print('You have not added any absorption sources! Your loss. PROMETHEUS exits now.')
    sys.exit()

"""
Absorbers
"""
lines_dict = {}
species_dict = {}

while True:
    if check('only_rayleigh_scatt', scenario_dict):
        add_lines_ON = read_str('Do you want to add absorption lines on top of Rayleigh scattering?', ['yes', 'no'])

        if add_lines_ON == 'no':
            break

    absorber_mode = read_str('How do you want to specify the absorption cross section for the absorption lines?', ['BuiltinLine']) #ImportLineList

    if absorber_mode == 'BuiltinLine':
        
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

        if len(lines_dict) == 0:
            print('You have not added any source of absorption! Your loss. PROMETHEUS exits now.')
            sys.exit()


        for value in lines_dict.values():
            element = value[4]
            species_dict[element] = 0
        
        for key_species in species_dict.keys():
            params = []
            for key_scenario in scenario_dict.keys():

                if not check('evaporative', key_scenario, scenario_dict[key_scenario]):
                    chi = read_value('Enter the mixing ratio of ' + key_species + ' in the ' + key_scenario + ' scenario:', 
                    1e-15, 1, 1)

                    if check('atmosphere', key_scenario):
                        T = scenario_dict[key_scenario][0]
                    if key_scenario == 'escaping':
                        T = scenario_dict[key_scenario][2]

                    params.append([key_scenario, chi, T])
                        
                else:
                    N = read_value('Enter the number of absorbing ' + key_species + ' atoms in the ' + key_scenario + ' scenario:', 
                    1e10, 1e50, 1)

                    if 'therm_broad' in scenario_dict[key_scenario]:

                        v_mean = read_value('Enter the mean (thermal) velocity of ' + key_species + ' in the ' + key_scenario + '\
 scenario in km/s:', 1e-3, 1e5, 1e5)
                        absorber_mass = speciesMass_dict[key_species]

                        params.append([key_scenario, N, v_mean, absorber_mass])

                    else:

                        params.append([key_scenario, N])                   

            species_dict[key_species] = params

    break

"""
Performance parameters
"""

if mode == 'spectrum':
    lower_w = read_value('Enter the lower wavelength border in Angstrom:', 1e-3, 1e12, 1e-8)
    upper_w = read_value('Enter the upper wavelength border in Angstrom:', lower_w * 1e8, 1e12, 1e-8,
    round = False)
    resolution = read_value('Enter the resolution of the wavelength grid in Angstrom:', 1e-6,
    (upper_w - lower_w) * 1e8 / 2., 1e-8, round = False)

upper_x = read_value('Enter the half chord length (x-direction) for the numerical integration along the x-axis in planetary radii:', 0, a_p / R_0, R_0, round = False)
x_steps = read_value('Enter the steps for the spatial discretization along the chord (x-direction):', 2, 1e6, 1)

if check('spherical_symmetry', scenario_dict):
    phi_steps = 0
else:
    phi_steps = read_value('Enter the steps for the spatial discretization for the polar coordinate (phi-direction):', 2, 1e4, 1)

z_steps = read_value('Enter the steps for the spatial discretization in z-direction:', 2, 1e6, 1)


"""
Additional output
"""

if 'barometric' in scenario_dict.keys():
    benchmark = read_str('Do you want to record the analytical benchmark for the barometric scenario?', ['yes', 'no'])
else:
    benchmark = 'no'

"""
Write parameter dictionary and store it as json file
"""


parameters = {'R_star': R_star, 'R_0': R_0, 'M_p': M_p, 'a_p': a_p, 
'R_moon': R_moon, 'a_moon': a_moon, 'alpha_moon': alpha_moon, 'z_moon': z_moon,
'direction': direction, 'mode': mode, 'dishoom_import': dishoom_import,
'Scenarios': scenario_dict, 'Lines': lines_dict, 'Species': species_dict,
'lower_w': lower_w, 'upper_w': upper_w, 'resolution': resolution, 
'upper_x': upper_x, 'x_steps': x_steps, 'phi_steps': phi_steps, 'z_steps': z_steps, 'benchmark': benchmark}


with open('../settings.txt', 'w') as outfile:
    json.dump(parameters, outfile)

