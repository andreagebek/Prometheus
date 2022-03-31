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
SCRIPTPATH = os.path.realpath(__file__)
GITPATH = os.path.dirname(SCRIPTPATH)
PARENTPATH = os.path.dirname(GITPATH)
sys.path.append(GITPATH)
import prometheusScripts.constants as const

if __name__ == '__main__':

    N_arguments = len(sys.argv)

    if N_arguments == 1:

        import prometheusScripts.setup

        sys.exit(0)

    import prometheusScripts.fluxDecrease as flux

    startTime = datetime.now()

    """
    Read in the json parameters and perform the radiative transfer calculation
    """

    paramsFilename = sys.argv[1]

    with open(PARENTPATH + '/setupFiles/' + paramsFilename + '.txt') as file:
        param = json.load(file)

    fundamentalsDict = param['Fundamentals']
    scenarioDict = param['Scenarios']
    architectureDict = param['Architecture']
    speciesDict = param['Species']
    gridsDict = param['Grids']
    outputDict = param['Output']

    GRID, args = flux.prepareArguments(fundamentalsDict, architectureDict, scenarioDict, speciesDict, gridsDict, outputDict, startTime)

    N_cores = mp.cpu_count()

    print('Starting main calculation on', N_cores, 'cores:', datetime.now() - startTime)

    RESULTS = []

    with mp.Pool(processes = N_cores) as pool:

        RESULTS_R, RESULTS_TAU = zip(*pool.map(partial(flux.evaluateChord, args = args), GRID))

    pool.close()
    pool.join()

    print('Main calculation finished:', datetime.now() - startTime)

    Rarray = np.array(RESULTS_R).reshape((int(gridsDict['phi_steps']) * int(gridsDict['rho_steps']), int(gridsDict['orbphase_steps']), len(args[3])))
    Runflattened = np.sum(Rarray, axis = 0).T # Sum over rho and phi coordinates
    index_max = np.unravel_index(np.argmin(Runflattened, axis = None), Runflattened.shape)

    R = Runflattened.flatten()

    """
    Store the output in .txt files
    """

    orbphase_axis = flux.constructAxis(gridsDict, architectureDict, 'orbphase')
    wavelength_axis = flux.constructAxis(gridsDict, architectureDict, 'wavelength') * 1e8 # Conversion from cm to Angstrom

    wavelength, orbphase = np.meshgrid(wavelength_axis, orbphase_axis, indexing = 'ij')
    wavelength = wavelength.flatten()
    orbphase = orbphase.flatten()


    header = 'Wavelength grid (Å), Orbital phase grid [rad], R'

    np.savetxt(PARENTPATH + '/output/' + paramsFilename + '_lightcurve.txt', np.array([wavelength, orbphase, R]).T, header = header)

    
    if outputDict['recordTau']:

        rho_axis = flux.constructAxis(gridsDict, architectureDict, 'rho')

        phi_axis = flux.constructAxis(gridsDict, architectureDict, 'phi')
        
        phi, rho = np.meshgrid(phi_axis, rho_axis, indexing = 'ij')
        phi = phi.flatten()
        rho = rho.flatten()

        TAUarray = np.array(RESULTS_TAU).reshape((int(gridsDict['phi_steps']), int(gridsDict['rho_steps']), int(gridsDict['orbphase_steps']), len(args[3])))
        tauDisk = TAUarray[:, :, index_max[1], index_max[0]].flatten()

        np.savetxt(PARENTPATH + '/output/' + paramsFilename + '_tau.txt', np.array([phi, rho, tauDisk]).T, header = 'phi grid [rad], rho grid [cm], tau')  
    
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