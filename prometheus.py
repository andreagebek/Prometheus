# coding=utf-8
"""
Main code to run the radiative transfer calculation.
Created on 15. July 2021 by Andrea Gebek.
"""

import sys
import json
import numpy as np
from datetime import datetime
import os
import multiprocessing as mp
from functools import partial
import pythonScripts.setup as setup
import pythonScripts.gasProperties as gasprop
import pythonScripts.celestialBodies as bodies
import pythonScripts.geometryHandler as geom
import pythonScripts.constants as const


PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) # This is the folder that encloses the github Prometheus folder

if __name__ == '__main__':

    N_arguments = len(sys.argv)

    # If main.py is ran with the setup argument, run the setup script

    if 'setup' in sys.argv:

        setup.createSetupFile(PATH)

        sys.exit(0)

    # If there is no setup argument, run the main script

    startTime = datetime.now()

    """
    First, set up the classes from the json input file created during setup
    """

    paramsFilename = sys.argv[1]

    with open(PATH + '/setupFiles/' + paramsFilename + '.txt') as file:
        param = json.load(file)

    fundamentalsDict = param['Fundamentals']
    scenarioDict = param['Scenarios']
    architectureDict = param['Architecture']
    speciesDict = param['Species']
    gridsDict = param['Grids']


    planet = bodies.AvailablePlanets().findPlanet(architectureDict['planetName'])
    
    wavelengthGrid = gasprop.WavelengthGrid(gridsDict['lower_w'], gridsDict['upper_w'], gridsDict['widthHighRes'], gridsDict['resolutionLow'], gridsDict['resolutionHigh'])
    spatialGrid = geom.Grid(gridsDict['x_midpoint'], gridsDict['x_border'], int(gridsDict['x_steps']), gridsDict['upper_rho'], int(gridsDict['rho_steps']),
                            int(gridsDict['phi_steps']), gridsDict['orbphase_border'], int(gridsDict['orbphase_steps']))

    scenarioList = []

    for key_scenario in scenarioDict.keys():
        if key_scenario == 'barometric':
            scenarioList.append(gasprop.BarometricAtmosphere(scenarioDict['barometric']['T'], scenarioDict['barometric']['P_0'], scenarioDict['barometric']['mu'], planet))

        elif key_scenario == 'hydrostatic':
            scenarioList.append(gasprop.HydrostaticAtmosphere(scenarioDict['hydrostatic']['T'], scenarioDict['hydrostatic']['P_0'], scenarioDict['hydrostatic']['mu'], planet))

        elif key_scenario == 'powerLaw':
            if 'P_0' in scenarioDict['powerLaw'].keys():
                scenarioList.append(gasprop.PowerLawAtmosphere(scenarioDict['powerLaw']['T'], scenarioDict['powerLaw']['P_0'], scenarioDict['powerLaw']['q_esc'], planet))
            else:
                key_species = speciesDict['powerLaw'].keys()[0] # Only one absorber
                Nparticles = speciesDict['powerLaw'][key_species]['Nparticles']
                scenarioList.append(gasprop.PowerLawExosphere(Nparticles, scenarioDict['powerLaw']['q_esc'], planet))

        elif key_scenario == 'exomoon':
            moon = bodies.Moon(architectureDict['starting_orbphase_moon'], architectureDict['R_moon'], architectureDict['a_moon'], planet)
            key_species = list(speciesDict['exomoon'].keys())[0] # Only one absorber
            Nparticles = speciesDict['exomoon'][key_species]['Nparticles']
            scenarioList.append(gasprop.MoonExosphere(Nparticles, scenarioDict['exomoon']['q_moon'], moon))
        
        elif key_scenario == 'torus':
            key_species = list(speciesDict['torus'].keys())[0] # Only one absorber
            Nparticles = speciesDict['torus'][key_species]['Nparticles']
            scenarioList.append(gasprop.TorusExosphere(Nparticles, scenarioDict['torus']['a_torus'], scenarioDict['torus']['v_ej'], planet))
        
        elif key_scenario == 'serpens':
            key_species = list(speciesDict['serpens'].keys())[0] # Only one absorber
            Nparticles = speciesDict['serpens'][key_species]['Nparticles']
            scenario = gasprop.SerpensExosphere(scenarioDict['serpens']['serpensPath'], Nparticles, planet, 0.)
            scenario.addInterpolatedDensity(spatialGrid)
            scenarioList.append(scenario) # sigmaSmoothing hardcoded to 0, i.e. no Gaussian smoothing of the serpens density distribution.



    for idx, key_scenario in enumerate(scenarioDict.keys()):

        if 'T' in scenarioDict[key_scenario].keys():
            for key_species in speciesDict[key_scenario].keys():
                absorberDict = speciesDict[key_scenario][key_species]
                if key_species in const.AvailableSpecies().listSpeciesNames(): # Atom/ion
                    scenarioList[idx].addConstituent(key_species, absorberDict['chi'])
                    scenarioList[idx].constituents[-1].addLookupFunctionToConstituent(wavelengthGrid)
                else:
                    scenarioList[idx].addMolecularConstituent(key_species, absorberDict['chi'])
                    scenarioList[idx].constituents[-1].addLookupFunctionToConstituent()
        
        else: # Evaporative scenario
            for key_species in speciesDict[key_scenario].keys():
                absorberDict = speciesDict[key_scenario][key_species]
                if key_species in const.AvailableSpecies().listSpeciesNames(): # Atom/ion
                    scenarioList[idx].addConstituent(key_species, absorberDict['sigma_v'])
                    scenarioList[idx].constituents[-1].addLookupFunctionToConstituent(wavelengthGrid)
                else:
                    scenarioList[idx].addMolecularConstituent(key_species, absorberDict['T'])
                    scenarioList[idx].constituents[-1].addLookupFunctionToConstituent()

    atmos = gasprop.Atmosphere(scenarioList, fundamentalsDict['DopplerOrbitalMotion'])

    main = gasprop.Transit(atmos, wavelengthGrid, spatialGrid)
    main.addWavelength()


    """
    Execute forward model
    """

    R = main.sumOverChords()
    wavelength = wavelengthGrid.constructWavelengthGrid(scenarioList)
    orbphase = spatialGrid.constructOrbphaseAxis()

    """
    Store the result
    """

    orbphaseOutput = np.insert(orbphase / (2. * np.pi), 0, np.nan)
    mainOutput = np.vstack((wavelength, R))
    output = np.vstack((orbphaseOutput, mainOutput.T))

    header = 'Prometheus output file.\nFirst row: Orbital phases [1]\n\
All other rows: Wavelength [cm] (first column), Transit depth R(orbital phase, wavelength) [1] (other columns)'
    
    np.savetxt(PATH + '/output/' + paramsFilename + '.txt', output, header = header)

    elapsedTime = datetime.now() - startTime

    print("\nPROMETHEUS finished, yay! Elapsed time is:", elapsedTime)

    print("The maximal flux decrease due to atmospheric/exospheric absorption in percent is:", np.abs(np.round(100 * (1 - np.min(R)), 5)))

    print("The minimal flux decrease due to atmospheric/exospheric absorption in percent is:", np.abs(np.round(100 * (1 - np.max(R)), 5)))


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