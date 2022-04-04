# coding=utf-8
"""
Run a question and answer session and store a settings.txt 
file to run PROMETHEUS.
Created on 2. June 2021 by Andrea Gebek.
"""

import numpy as np
import json
import sys
import os
SCRIPTPATH = os.path.realpath(__file__)
GITPATH = os.path.dirname(os.path.dirname(SCRIPTPATH))
PARENTPATH = os.path.dirname(GITPATH)
sys.path.append(GITPATH)

import prometheusScripts.constants as const


def read_value(text, lower, upper, unit, roundBorders = True, digits = 0, acceptLowerBorder = False, acceptUpperBorder = False):
    
    while True:

        if acceptLowerBorder:
            brackl = ' ['
        else:
            brackl = ' ('

        if acceptUpperBorder:
            brackr = '] '
        else:
            brackr = ') '

        
        if roundBorders:
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

        elif value == lower and acceptLowerBorder:
            return value * unit
            break

        elif value == upper and acceptUpperBorder:
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
Fundamentals
"""

fundamentalsDict = {'PlanetarySource': False, 'ExomoonSource': False, 'DopplerPlanetRotation': False}

print(r"""
 *******  *******           ,/MMM8&&&.         ****     **** ******** ********** **      ** ******** **     **  ********
/**////**/**////**     _...MMMMM88&&&&..._    /**/**   **/**/**///// /////**/// /**     /**/**///// /**    /** **////// 
/**   /**/**   /**   .:'''MMMMM88&&&&&&''':.  /**//** ** /**/**          /**    /**     /**/**      /**    /**/**       
/******* /*******   :     MMMMM88&&&&&&     : /** //***  /**/*******     /**    /**********/******* /**    /**/*********
/**////  /**///**    ':...MMMMM88&&&&&&....:  /**  //*   /**/**////      /**    /**//////**/**////  /**    /**////////**
/**      /**  //**      `''MMMMM88&&&&'''`    /**   /    /**/**          /**    /**     /**/**      /**    /**       /**
/**      /**   //**         'MMM8&&&'         /**        /**/********    /**    /**     /**/********//*******  ******** 
//       //     //                            //         // ////////     //     //      // ////////  ///////  ////////  
""")

print('\nWelcome to PROMETHEUS! First, define if you want to make some fundamental simplifications.\n')


fundamentalsDict['CLV_variations'] = read_str('Do you want to take center-to-limb variations into account?', ['yes', 'no'])
fundamentalsDict['RM_effect'] = read_str('Do you want to take the Rossiter-McLaughlin-Effect into account (note that this means that you \
have to provide additional information about the host star to specifiy its spectrum)?', ['yes', 'no'])
fundamentalsDict['DopplerOrbitalMotion'] = read_str('Do you want to consider the Doppler shifts due to planetary/exomoon orbital motion?', ['yes', 'no'])
fundamentalsDict['ExactSigmaAbs'] = read_str('Do you want to calculate the absorption cross section exactly (NOT RECOMMENDED, makes the code very slow)?', ['yes', 'no'])

if not fundamentalsDict['ExactSigmaAbs']:
    fundamentalsDict['LookupResolution'] = read_value('Enter the resolution for the lookup table in Angstrom. Note that this effectively corresponds to a velocity resolution, \
to resolve a velocity of 1 km/s at 6000 Angstrom a lookup table resolution of 0.02 Angstrom is needed.', 1e-5, 100, 1e-8)

"""
Scenarios for the spatial distribution of the medium
"""

print('\nNow, specifiy the spatial distribution of the absorbing medium, i.e. the structure of the atmosphere/exosphere.\n')

scenarioDict = {}
PossibleScenarios = ['barometric', 'hydrostatic', 'escaping', 'exomoon', 'torus', 'plasma', '0'] # Possible scenarios

while True:

    scenario_name = read_str('Enter the name of the scenario for the number density profile or \
0 to stop adding absorption sources:', PossibleScenarios)

    if scenario_name == '0':
        break

    params = {}

    if scenario_name == 'barometric' or scenario_name == 'hydrostatic':

        fundamentalsDict['PlanetarySource'] = True

        params['T'] = read_value('Enter the temperature of the atmosphere in Kelvin:', 1, 1e6, 1)
        params['P_0'] = read_value('Enter the pressure at the reference radius in bar:', 1e-12, 1e6, 1e6)
        params['mu'] = read_value('Enter the mean molecular weight of the atmosphere in atomic mass units:', 0.05, 500, const.amu)
        params['RayleighScatt'] = read_str('Do you want to add Rayleigh scattering in this scenario?', ['yes', 'no'])

        if scenario_name == 'barometric':
            PossibleScenarios.remove('hydrostatic')
        else:
            PossibleScenarios.remove('barometric')



    elif scenario_name == 'escaping':

        fundamentalsDict['PlanetarySource'] = True

        params['q_esc'] = read_value('Enter the power law index for the escaping atmosphere:', 3, 20, 1)

        if read_str('Do you want to normalize the number density profile via (total) pressure or via number \
of absorbing atoms at the base of the wind?', ['pressure', 'number']) == 'pressure':

            params['P_0'] = read_value('Enter the pressure at the base of the wind in bar:', 1e-15, 1e3, 1e6)
            params['T'] = read_value('Enter the temperature of the escaping wind in Kelvin:', 1, 1e6, 1)

        params['RadialWind'] = read_str('Do you want to add radially escaping winds in the escaping scenario?', ['yes', 'no'])


        if params['RadialWind']:

            params['vRadial_0'] = read_value('Enter the velocity of the radially escaping wind at the reference radius in km/s:', 1e-3, 1e3, 1e5)


    elif scenario_name == 'exomoon':

        params['q_moon'] = read_value('Enter the power law index for the exomoon exosphere:', 3, 20, 1)

        params['RadialWind'] = read_str('Do you want to add radially escaping winds in the exomoon scenario?', ['yes', 'no'])

        if params['RadialWind']:

            params['vRadial_0'] = read_value('Enter the velocity of the radially escaping wind at the exomoon radius in km/s:', 1e-3, 1e3, 1e5)   

        fundamentalsDict['ExomoonSource'] = True

    elif scenario_name == 'torus':

        params['a_torus'] = read_value('Enter the distance between the center of the torus and the center of the exoplanet in Jovian radii:', 0.01, 1000, const.R_J)
        params['v_ej'] = read_value('Enter the ejection velocity (which sets the torus scale height) of the particles from the torus in km/s:', 1e-2, 1e3, 1e5)

    elif scenario_name == 'plasma':
        params['amitisFilename'] = read_str('Enter the name of the AMITIS output file (located in "amitis_outputs" folder, without the .5h ending:')

    if 'T' in params.keys():

        params['thermBroad'] = read_str('Do you want to add thermal line broadening in this scenario?', ['yes', 'no'])

    else:

        params['thermBroad'] = read_str('Do you want to add pseudo-thermal line broadening via mean velocity of the absorber in this scenario?', ['yes', 'no'])


    scenarioDict[scenario_name] = params

    PossibleScenarios.remove(scenario_name)


if len(scenarioDict) == 0:
    print('You have not added any absorption sources! Your loss. PROMETHEUS exits now.')
    sys.exit()


"""
Architecture
"""

print('\nProvide parameters related to the architecture of the system.\n')

PossibleSystems = ['0']

for key in const.planetsDict.keys():
    PossibleSystems.append(key)


architectureDict = {}
systemname = read_str('Enter the name of the exoplanetary system or zero if parameters are entered manually:', PossibleSystems)

if systemname == '0':

    architectureDict['R_star'] = read_value('Enter the radius of the host star in solar radii:', 1e-5, 1e5, const.R_sun)
    architectureDict['R_0'] = read_value('Enter the radius of the exoplanet in Jupiter radii:', 1e-5, 1e5, const.R_J)
    architectureDict['M_p'] = read_value('Enter the mass of the exoplanet in Jupiter masses:', 1e-5, 1e3, const.M_J)
    architectureDict['a_p'] = read_value('Enter the orbital distance between planet and star in AU:', 1e-5, 1e3, const.AU)
    architectureDict['M_star'] = read_value('Enter the mass of the host star in solar masses:', 1e-5, 1e10, const.M_sun)

    if fundamentalsDict['RM_effect']:

        architectureDict['T_eff'] = read_value('Enter the effective temperature of the star in Kelvin:', 2300, 12000, 1, roundBorders = False, acceptLowerBorder = True, acceptUpperBorder = True)
        architectureDict['log_g'] = read_value('Enter the logarithmic value of the surface gravity in log10(cm/s^2):', 0, 6, 1, roundBorders = False, acceptLowerBorder = True, acceptUpperBorder = True)
        architectureDict['Fe_H'] = read_value('Enter the metallicity of the star [Fe/H]:', -4, 1, 1, roundBorders = False, acceptLowerBorder = True, acceptUpperBorder = True)

        if architectureDict['Fe_H'] > -3.5 and architectureDict['Fe_H'] < 0.25: # The PHOENIX library has varying alpha-enhancement only for this metallicity range
            architectureDict['alpha_Fe'] = read_value('Enter the alpha-enhancement of the star [alpha/Fe]:', -0.2, 1.2, 1, roundBorders = False, acceptLowerBorder = True, acceptUpperBorder = True)
        else:
            architectureDict['alpha_Fe'] = 0

else:
    architectureDict['R_star'] = const.planetsDict[systemname][0]
    architectureDict['M_star'] = const.planetsDict[systemname][1]
    architectureDict['R_0'] = const.planetsDict[systemname][2]
    architectureDict['M_p'] = const.planetsDict[systemname][3]
    architectureDict['a_p'] = const.planetsDict[systemname][4]
    architectureDict['T_eff'] = const.planetsDict[systemname][5]
    architectureDict['log_g'] = const.planetsDict[systemname][6]
    architectureDict['Fe_H'] = const.planetsDict[systemname][7]
    architectureDict['alpha_Fe'] = const.planetsDict[systemname][8]


if fundamentalsDict['RM_effect']:

    architectureDict['period_starrot'] = read_value('Enter the period of the stellar rotation in days:', 0, 1000, 86400)
    architectureDict['inclination_starrot'] = read_value('Enter the inclination of the stellar rotation against the planetary orbital plane in degrees:', 0, 180, np.pi / 180., roundBorders = False, acceptLowerBorder = True)
    architectureDict['azimuth_starrot'] = read_value('Enter the angle between the line of sight from the observer to the system and the angular momentum vector of the stellar rotation \
projected onto the orbital plane in degrees', 0, 360, np.pi / 180., roundBorders = False, acceptLowerBorder = True)


if fundamentalsDict['CLV_variations']:

    architectureDict['u1'] = read_value('Enter the first (linear) coefficient for limb darkening:', -1, 1, 1, acceptLowerBorder = True, acceptUpperBorder = True)
    architectureDict['u2'] = read_value('Enter the second (quadratic) coefficient for limb darkening:', -1, 1, 1, acceptLowerBorder = True, acceptUpperBorder = True)


if fundamentalsDict['ExomoonSource']:

    architectureDict['R_moon'] = read_value('Enter the radius of the moon in Io radii:', 1e-3, 1e3, const.R_Io)

    architectureDict['a_moon'] = read_value('Enter the orbital distance between the exomoon and the planet in planetary radii (measured from the centers of the bodies):', 
1, architectureDict['a_p'] / architectureDict['R_0'], architectureDict['R_0'], roundBorders = False)

    architectureDict['starting_orbphase_moon'] = read_value('Enter the orbital phase of the moon when the planet is transiting. A moon orbital phase of 0 corresponds to the moon sitting \
between the planet and the observer, 0.25 means that the exomoon is located to the right of the planet when viewed from the observer.', -0.5, 0.5, 2. * np.pi,
acceptLowerBorder = True, acceptUpperBorder = True)

    """

        i_moon = read_value('Enter the inclination of the moon orbital plane with respect to the planet orbital plane in degrees:', 0, 90, 2. * np.pi / 360., accept_borders = True)
        
        if i_moon > 0:
            iphase_moon = read_value('Enter the orbital phase of the moon when it ascends through the planet orbital plane:', 0)
    """

if fundamentalsDict['PlanetarySource']:

    fundamentalsDict['DopplerPlanetRotation'] = read_str('Do you want to consider the Doppler shifts due to planetary rotation? This only applies to planetary sources \
(barometric, hydrostatic, escaping).', ['yes', 'no'])

    if fundamentalsDict['DopplerPlanetRotation']:

        architectureDict['period_planetrot'] = read_value('Enter the period for the planetary rotation in days:', 0, 1000, 86400)


"""
Specify the absorption species.
"""

print('\nSpecify the absorbing species and their abundances.\n')

speciesDict = {}


for key_scenario in scenarioDict.keys():

    PossibleAbsorbers = list(const.speciesInfoDict.keys())
    PossibleAbsorbers.append('0')

    speciesDict[key_scenario] = {}

    if 'T' in scenarioDict[key_scenario].keys(): # Scenarios that incorporate a temperature

        while True:

            params = {}
            key_species = read_str('Enter the name of the absorbing species you want to consider for the ' + key_scenario + ' scenario, or 0 to move on:', PossibleAbsorbers)

            if key_species == '0':
                break

            else:

                params['chi'] = read_value('Enter the mixing ratio of ' + key_species + ' in the ' + key_scenario + ' scenario:', 
                0, 1, 1, acceptLowerBorder = True, acceptUpperBorder = True)

            if scenarioDict[key_scenario]['thermBroad']: # Velocity dispersion only for atoms/ions

                params['sigma_v'] = np.sqrt(const.k_B * scenarioDict[key_scenario]['T'] / const.speciesInfoDict[key_species][2])

            else:

                params['sigma_v'] = 0

            speciesDict[key_scenario][key_species] = params
            PossibleAbsorbers.remove(key_species)

    else: # Evaporative scenarios (no temperature)

        while True:

            params = {}
            key_species = read_str('Enter the name of the absorbing species you want to consider for the ' + key_scenario + ' scenario, or 0 to move on:', PossibleAbsorbers)
            
            if key_species == '0':
                break

            else:

                params['chi'] = read_value('Enter the number of absorbing ' + key_species + ' atoms/ions in the ' + key_scenario + ' scenario:', 
                0, 1e50, 1, acceptLowerBorder = True, acceptUpperBorder = True)


            if scenarioDict[key_scenario]['thermBroad']:

                params['sigma_v'] = read_value('Enter the pseudo-thermal velocity dispersion (sigma_v = sqrt(k_B * T / m))  of ' + key_species + ' in the ' + key_scenario + ' \
scenario in km/s:', 1e-3, 1e5, 1e5)

            else:

                params['sigma_v'] = 0

            speciesDict[key_scenario][key_species] = params
            PossibleAbsorbers.remove(key_species)



"""
Grid parameters
"""

print('\nAlmost done! Specify the discretization parameters for the wavelength and spatial grids. Default values give an expected runtime of ~10s.\n')

gridsDict = {}

gridsDict['lower_w'] = read_value('Enter the lower wavelength border in Angstrom:', 500, 55000, 1e-8, roundBorders = False)
gridsDict['upper_w'] = read_value('Enter the upper wavelength border in Angstrom:', gridsDict['lower_w'] * 1e8, 55000, 1e-8, roundBorders = False)
gridsDict['resolution'] = read_value('Enter the resolution of the wavelength grid in Angstrom \
(default: ' + str((gridsDict['upper_w'] - gridsDict['lower_w']) * 1e8 / 20.) + '):', 1e-6, (gridsDict['upper_w'] - gridsDict['lower_w']) * 1e8 / 2., 1e-8, roundBorders = False)

gridsDict['x_midpoint'] = read_value('Enter the midpoint of the grid in x-direction in planetary orbital radii:', architectureDict['R_star'] / architectureDict['a_p'], 10.,
architectureDict['a_p'], roundBorders = False)
gridsDict['x_border'] = read_value('Enter the half chord length (x-direction) for the numerical integration along the x-axis in planetary radii:',
0., (gridsDict['x_midpoint'] - architectureDict['R_star']) / architectureDict['R_0'], architectureDict['R_0'], roundBorders = False)
gridsDict['x_steps'] = read_value('Enter the number of bins for the spatial discretization along the chord (x-direction, default: 20):', 2, 1e6, 1)

gridsDict['phi_steps'] = read_value('Enter the number of bins for the spatial discretization for the polar coordinate (phi-direction, default: 20):', 1, 1e4, 1, acceptLowerBorder = True)
gridsDict['rho_steps'] = read_value('Enter the number of bins for the spatial discretization in z-direction (default: 20):', 2, 1e6, 1)

gridsDict['orbphase_border'] = read_value('Enter the orbital phase at which the light curve calculation starts and stops:', 0, 0.5, 2 * np.pi, acceptLowerBorder = True, acceptUpperBorder = True)
gridsDict['orbphase_steps'] = read_value('Enter the number of bins for the orbital phase discretization (default: 20):', 1, 1e4, 1, acceptLowerBorder = True)

"""
Additional output
"""

print('\nFinally, specify some output options.\n')

outputDict = {}

paramsFilename = read_str('How do you want to name the txt file containing the parameters of this session (enter the file name without the .txt ending)?')

outputDict['recordTau'] = read_str('Do you want to record the optical depth for all chords of the spatial grid at the wavelength and orbital phase with the largest flux decrease?', ['yes', 'no'])


"""
Write parameter dictionary and store it as json file
"""

print('\n\nAll parameters are stored! To run PROMETHEUS, type <python main.py ' + paramsFilename + '>.\n\n')

print(r"""
 *******  *******           ,/MMM8&&&.         ****     **** ******** ********** **      ** ******** **     **  ********
/**////**/**////**     _...MMMMM88&&&&..._    /**/**   **/**/**///// /////**/// /**     /**/**///// /**    /** **////// 
/**   /**/**   /**   .:'''MMMMM88&&&&&&''':.  /**//** ** /**/**          /**    /**     /**/**      /**    /**/**       
/******* /*******   :     MMMMM88&&&&&&     : /** //***  /**/*******     /**    /**********/******* /**    /**/*********
/**////  /**///**    ':...MMMMM88&&&&&&....:  /**  //*   /**/**////      /**    /**//////**/**////  /**    /**////////**
/**      /**  //**      `''MMMMM88&&&&'''`    /**   /    /**/**          /**    /**     /**/**      /**    /**       /**
/**      /**   //**         'MMM8&&&'         /**        /**/********    /**    /**     /**/********//*******  ******** 
//       //     //                            //         // ////////     //     //      // ////////  ///////  ////////  
""")

parameters = {'Fundamentals': fundamentalsDict, 'Architecture': architectureDict, 'Scenarios': scenarioDict, 'Species': speciesDict, 'Grids': gridsDict, 'Output': outputDict}


with open(PARENTPATH +'/setupFiles/' + paramsFilename + '.txt', 'w') as outfile:
    json.dump(parameters, outfile)