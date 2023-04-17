# coding=utf-8
"""
Run a question and answer session and store a settings.txt 
file to run PROMETHEUS.
Created on 2. June 2021 by Andrea Gebek.
"""

import numpy as np
import json
import sys
import pythonScripts.constants as const
import pythonScripts.celestialBodies as bodies

class NumericalQuestion:

    def __init__(self, questionText, lowerBorder, upperBorder, unit, roundBorders = True, digits = 0, acceptLowerBorder = True, acceptUpperBorder = True):
        self.questionText = questionText
        self.lowerBorder = lowerBorder
        self.upperBorder = upperBorder
        self.unit = unit
        self.roundBorders = roundBorders
        self.digits = digits
        self.acceptLowerBorder = acceptLowerBorder
        self.acceptUpperBorder = acceptUpperBorder

    def readValue(self):
        
        while True:

            if self.acceptLowerBorder:
                brackl = ' ['
            else:
                brackl = ' ('

            if self.acceptUpperBorder:
                brackr = '] '
            else:
                brackr = ') '

            
            if self.roundBorders:
                string = input(self.questionText + brackl + '{:.0e}'.format(self.lowerBorder) + ', ' + '{:.0e}'.format(self.upperBorder) + brackr).replace(' ', '')

            elif self.digits > 0:
                string = input(self.questionText + brackl + '{:.{n}e}'.format(self.lowerBorder, n = self.digits) + ', ' + '{:.{n}e}'.format(self.upperBorder, n = self.digits) + brackr).replace(' ', '')

            else:
                string = input(self.questionText + brackl + '{:e}'.format(self.lowerBorder) + ', ' + '{:e}'.format(self.upperBorder) + brackr).replace(' ', '')

            if string == '':
                print('You actually do have to enter something!')
                continue
            
            value = float(string)

            if value > self.lowerBorder and value < self.upperBorder:
                return value * self.unit

            elif value == self.lowerBorder and self.acceptLowerBorder:
                return value * self.unit

            elif value == self.upperBorder and self.acceptUpperBorder:
                return value * self.unit

            else:
                print('The value you entered is not in the appropriate interval.')
                continue


class TextQuestion:

    def __init__(self, questionText, boolQuestion = False, options = None):
        self.questionText = questionText
        if boolQuestion:
            self.options = ['yes', 'no']
        else:
            self.options = options

    def readStr(self):

        while True:

            if self.options is None:
                string = input(self.questionText + ' ').replace(' ', '')

            else:
                string = input(self.questionText + ' ' + str(self.options).replace("'",'') + ' ').replace(' ', '')

            if string == '':
                print('You actually do have to enter something!')
                continue
            
            elif self.options is not None and string not in self.options:
                print('Please select a valid option.')
                continue

            else:
                if string == 'yes':
                    return True
                elif string == 'no':
                    return False
                else:
                    return string




"""
Fundamentals
"""

QCLV = TextQuestion('Do you want to take center-to-limb variations into account?', True)
QRM = TextQuestion('Do you want to take the Rossiter-McLaughlin-Effect into account (note that this means that you \
have to provide additional information about the host star to specifiy its spectrum)?', True)
QOrbitalMotion = TextQuestion('Do you want to consider the Doppler shifts due to planetary/exomoon orbital motion?', True)

def setupFundamentals():

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

    fundamentalsDict = {'ExomoonSource': False, 'DopplerPlanetRotation': False}

    fundamentalsDict['CLV_variations'] = QCLV.readStr()
    fundamentalsDict['RM_effect'] = QRM.readStr()
    fundamentalsDict['DopplerOrbitalMotion'] = QOrbitalMotion.readStr()

    return fundamentalsDict


"""
Scenarios for the spatial distribution of the medium
"""

QScenario = TextQuestion('Enter the name of the scenario for the number density profile or \
0 to stop adding absorption sources:', False, ['barometric', 'hydrostatic', 'powerLaw', 'exomoon', 'torus', 'serpens', '0'])
QTemperature = NumericalQuestion('Enter the temperature of the atmosphere in Kelvin:', 1., 1e6, 1.)
QPressure = NumericalQuestion('Enter the pressure at the reference radius in bar:', 1e-15, 1e6, 1e6)
Qmu = NumericalQuestion('Enter the mean molecular weight of the atmosphere in atomic mass units:', 0.05, 500., const.amu)
QpowerLaw = NumericalQuestion('Enter the power law index for the escaping atmosphere:', 3., 20., 1., acceptLowerBorder = False)
QpowerLawNormalization = TextQuestion('Do you want to normalize the number density profile via (total) pressure or via number \
of absorbing atoms at the base of the wind?', False, ['pressure', 'number'])
QpowerLawMoon = NumericalQuestion('Enter the power law index for the exomoon exosphere:', 3., 20., 1.)
QtorusAxis = NumericalQuestion('Enter the distance between the center of the torus and the center of the exoplanet in Jovian radii:', 0.01, 1000, const.R_J)
QtorusEjectionSpeed = NumericalQuestion('Enter the ejection velocity (which sets the torus scale height) of the particles from the torus in km/s:', 1e-2, 1e3, 1e5)
QserpensPath = TextQuestion('Enter the full path to the serpens input txt file, including the filename with the .txt ending:')

def setupScenario(fundamentalsDict):

    print('\nNow, specifiy the spatial distribution of the absorbing medium, i.e. the structure of the atmosphere/exosphere.\n')

    scenarioDict = {}

    while True:

        scenario_name = QScenario.readStr()

        if scenario_name == '0':
            break

        params = {}

        if scenario_name == 'barometric' or scenario_name == 'hydrostatic':

            params['T'] = QTemperature.readValue()
            params['P_0'] = QPressure.readValue()
            params['mu'] = Qmu.readValue()

            if scenario_name == 'barometric':
                QScenario.options.remove('hydrostatic')
            else:
                QScenario.options.remove('barometric')

        
        elif scenario_name == 'powerLaw':

            params['q_esc'] = QpowerLaw.readValue()

            if QpowerLawNormalization.readStr() == 'pressure':

                params['P_0'] = QPressure.readValue()
                params['T'] = QTemperature.readValue()



        elif scenario_name == 'exomoon':

            params['q_moon'] = QpowerLawMoon.readValue()

            fundamentalsDict['ExomoonSource'] = True


        elif scenario_name == 'torus':

            params['a_torus'] = QtorusAxis.readValue()
            params['v_ej'] = QtorusEjectionSpeed.readValue()


        elif scenario_name == 'serpens':

            params['serpensPath'] = QserpensPath.readStr()


        scenarioDict[scenario_name] = params

        QScenario.options.remove(scenario_name)


    if len(scenarioDict) == 0:
        print('You have not added any absorption sources! Your loss. PROMETHEUS exits now.')
        sys.exit()

    return scenarioDict



"""
Architecture
"""

Qsystem = TextQuestion('Enter the name of the exoplanetary system', False, bodies.AvailablePlanets().listPlanetNames())
QRstar = NumericalQuestion('Enter the radius of the host star in solar radii:', 1e-5, 1e5, const.R_sun)
QMstar = NumericalQuestion('Enter the mass of the host star in solar masses:', 1e-5, 1e10, const.M_sun)
QRplanet = NumericalQuestion('Enter the radius of the exoplanet in Jupiter radii:', 1e-5, 1e5, const.R_J)
QMplanet = NumericalQuestion('Enter the mass of the exoplanet in Jupiter masses:', 1e-5, 1e3, const.M_J)
QorbitalRadiusPlanet = NumericalQuestion('Enter the orbital distance between planet and star in AU:', 1e-5, 1e3, const.AU)
QTstar = NumericalQuestion('Enter the effective temperature of the star in Kelvin:', 2300., 12000., 1., roundBorders = False)
QloggStar = NumericalQuestion('Enter the logarithmic value of the surface gravity in log10(cm/s^2):', 0., 6., 1., roundBorders = False)
QmetallicityStar = NumericalQuestion('Enter the metallicity of the star [Fe/H]:', -4., 1., 1., roundBorders = False)
QalphaStar = NumericalQuestion('Enter the alpha-enhancement of the star [alpha/Fe]:', -0.2, 1.2, 1., roundBorders = False)
QRMvsini = NumericalQuestion('Enter the maximum velocity of the stellar rotation (v*sin(i)) in km/s:', 1e-2, 1e4, 1e5)
QRMazimuth = NumericalQuestion('Enter the angle at which the maximum velocity towards the observer lies on the stellar disk \
in degrees (measured in the canonical prometheus coordinate system of the angular coordinate)', 0., 360., np.pi / 180., roundBorders = False)
QCLVu1 = NumericalQuestion('Enter the first (linear) coefficient for limb darkening:', -1., 1., 1.)
QCLVu2 = NumericalQuestion('Enter the second (quadratic) coefficient for limb darkening:', -1., 1., 1.)
QRmoon = NumericalQuestion('Enter the radius of the moon in Io radii:', 1e-3, 1e3, const.R_Io)
QorbphaseMoon = NumericalQuestion('Enter the orbital phase of the moon when the planet is transiting. A moon orbital phase of 0 corresponds to the moon sitting \
between the planet and the observer, 0.25 means that the exomoon is located to the right of the planet when viewed from the observer.', -0.5, 0.5, 2. * np.pi)


def setupArchitecture(fundamentalsDict):

    print('\nProvide parameters related to the architecture of the system.\n')

    architectureDict = {}
    
    planetName = Qsystem.readStr()

    architectureDict['planetName'] = planetName

    planet = bodies.AvailablePlanets().findPlanet(planetName)

    if fundamentalsDict['RM_effect']:

        architectureDict['vsini'] = QRMvsini.readValue()
        architectureDict['azimuth_starrot'] = QRMazimuth.readValue()


    if fundamentalsDict['CLV_variations']:

        architectureDict['u1'] = QCLVu1.readValue()
        architectureDict['u2'] = QCLVu2.readValue()


    if fundamentalsDict['ExomoonSource']:

        architectureDict['R_moon'] = QRmoon.readValue()
        architectureDict['a_moon'] = NumericalQuestion('Enter the orbital distance between the exomoon and the planet in planetary radii (measured from the centers of the bodies):', 
        1, planet.a / planet.R, planet.R, roundBorders = False).readValue()
        architectureDict['starting_orbphase_moon'] = QorbphaseMoon.readValue()

    return architectureDict



"""
Specify the absorption species.
"""        


QabsorberType = TextQuestion('Do you want to add an atomic/ionic or a molecular absorber (the atom option includes ions)?', False, ['atom', 'molecule'])
QstopAbsorbers = TextQuestion('Do you want to add another absorber?', True)
Qmolecule = TextQuestion('Enter the name of the molecular absorber:')
Qsigmav = NumericalQuestion('Enter the pseudo-thermal velocity dispersion (sigma_v = sqrt(k_B * T / m)) in km/s:', 1e-3, 1e5, 1e5)
QMoleculesTemperature = NumericalQuestion('Enter the pseudo-temperature for the molecular absorber in K:', 100., 3400., 1.)

def setupSpecies(scenarioDict):

    print('\nSpecify the absorbing species and their abundances.\n')

    speciesDict = {}


    for key_scenario in scenarioDict.keys():

        PossibleAbsorbers = const.AvailableSpecies().listSpeciesNames()

        speciesDict[key_scenario] = {}
        params = {}

        if 'T' in scenarioDict[key_scenario].keys(): # Scenarios that incorporate a temperature

            while True:

                absorberType = QabsorberType.readStr()

                if absorberType == 'atom':

                    key_species = TextQuestion('Enter the name of the absorbing species you want to consider for the ' + key_scenario + ' scenario:', False, PossibleAbsorbers).readStr()
                    #params['sigma_v'] = np.sqrt(const.k_B * scenarioDict[key_scenario]['T'] / const.AvailableSpecies().findSpecies(key_species).mass) # Velocity dispersion only for atoms/ions
                    PossibleAbsorbers.remove(key_species)
                
                else:
                    key_species = Qmolecule.readStr()

                params['chi'] = NumericalQuestion('Enter the mixing ratio of ' + key_species + ' in the ' + key_scenario + ' scenario:', 0., 1., 1., acceptLowerBorder = False).readValue()
                    
                speciesDict[key_scenario][key_species] = params

                if not QstopAbsorbers.readStr():
                    break

        else: # Evaporative scenarios (no temperature), only one absorbing species

            absorberType = QabsorberType.readStr()

            if absorberType == 'atom':
                key_species = TextQuestion('Enter the name of the absorbing species you want to consider for the ' + key_scenario + ' scenario:', False, PossibleAbsorbers).readStr()
                params['sigma_v'] = Qsigmav.readValue()
            
            else:
                key_species = Qmolecule.readStr()
                params['T'] = QMoleculesTemperature.readValue()

            params['Nparticles'] = NumericalQuestion('Enter the number of ' + key_species + ' particles in the ' + key_scenario + ' scenario:', 0., 1e50, 1., acceptLowerBorder = False).readValue()
  

            speciesDict[key_scenario][key_species] = params

    return speciesDict


"""
Grid parameters
"""

QlowerWavelengthBorder = NumericalQuestion('Enter the lower wavelength border in Angstrom:', 500, 55000, 1e-8)
QxBins = NumericalQuestion('Enter the number of bins for the spatial discretization along the chord (x-direction):', 2., 1e6, 1.)
QphiBins = NumericalQuestion('Enter the number of bins for the spatial discretization for the polar coordinate (phi-direction):', 1., 1e4, 1.)
QrhoBins = NumericalQuestion('Enter the number of bins for the spatial discretization in radial direction (rho-direction):', 2., 1e6, 1.)
QorbphaseBorder = NumericalQuestion('Enter the orbital phase at which the light curve calculation starts and stops:', 0., 0.5, 2. * np.pi)
QorbphaseBins = NumericalQuestion('Enter the number of bins for the orbital phase discretization:', 1., 1e4, 1.)

def setupGrid(architectureDict):

    print('\nAlmost done! Specify the discretization parameters for the wavelength and spatial grids. \n')

    planet = bodies.AvailablePlanets().findPlanet(architectureDict['planetName'])

    gridsDict = {}

    gridsDict['lower_w'] = QlowerWavelengthBorder.readValue()
    gridsDict['upper_w'] = NumericalQuestion('Enter the upper wavelength border in Angstrom:', gridsDict['lower_w'] * 1e8, 55000, 1e-8, roundBorders = False).readValue()
    gridsDict['resolutionLow'] = NumericalQuestion('Enter the resolution of the coarse wavelength grid in Angstrom:', 1e-6, (gridsDict['upper_w'] - gridsDict['lower_w']) * 1e8 / 2., 1e-8, roundBorders = False).readValue()
    gridsDict['widthHighRes'] = NumericalQuestion('Enter the bandwidth (centered on each absorption line) over which to consider a higher resolution in Angstrom:', 1e-6, (gridsDict['upper_w'] - gridsDict['lower_w']) * 1e8 / 2., 1e-8, roundBorders = False).readValue()
    gridsDict['resolutionHigh'] = NumericalQuestion('Enter the resolution of the fine wavelength grid in Angstrom:', 1e-6, gridsDict['widthHighRes'] * 1e8, 1e-8, roundBorders = False).readValue()

    gridsDict['x_midpoint'] = planet.a # Hardcoded default option is the orbital radius of the planetary orbit around the host star
    gridsDict['x_border'] = NumericalQuestion('Enter the half chord length (x-direction) for the numerical integration along the x-axis in planetary radii:',
    0., (gridsDict['x_midpoint'] - planet.hostStar.R) / planet.R, planet.R, roundBorders = False).readValue()
    gridsDict['x_steps'] = QxBins.readValue()

    gridsDict['phi_steps'] = QphiBins.readValue()
    gridsDict['rho_steps'] = QrhoBins.readValue()
    gridsDict['upper_rho'] = planet.hostStar.R # Hardcoded default option is the radius of the host star
    gridsDict['orbphase_border'] = QorbphaseBorder.readValue()
    gridsDict['orbphase_steps'] = QorbphaseBins.readValue()

    return gridsDict



"""
Write parameter dictionary and store it as json file
"""

def createSetupFile(PATH):

    inputFileName = TextQuestion('Write the name of this setup file (without file name ending):').readStr()

    fundamentalsDict = setupFundamentals()
    scenarioDict = setupScenario(fundamentalsDict)
    architectureDict = setupArchitecture(fundamentalsDict)
    speciesDict = setupSpecies(scenarioDict)
    gridsDict = setupGrid(architectureDict)

    parameters = {'Fundamentals': fundamentalsDict, 'Architecture': architectureDict, 'Scenarios': scenarioDict, 'Species': speciesDict, 'Grids': gridsDict}


    with open(PATH +'/setupFiles/' + inputFileName + '.txt', 'w') as outfile:
        json.dump(parameters, outfile)

    print('\n\nAll parameters are stored! To run PROMETHEUS, type <python prometheus.py ' + inputFileName + '>.\n\n')

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
