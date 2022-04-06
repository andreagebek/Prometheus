# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 15:22:03 2022

@author: Lucian
"""

import numpy as np
import json
import sys
import os
import h5py

SCRIPTPATH = os.path.realpath(__file__)
GITPATH = os.path.dirname(os.path.dirname(SCRIPTPATH))
PARENTPATH = os.path.dirname(GITPATH)

sys.path.append(GITPATH)

import prometheusScripts.constants as const
from scipy.interpolate import RegularGridInterpolator


def plasmaGrid(phi, rho, orbphase, xArray, key_species, scenarioDict, architectureDict, fundamentalsDict):

    amitisfile = scenarioDict['AmitisPlasma']['AmitisFilename']
    
    f = h5py.File(PARENTPATH + '/amitis_outputs/' + amitisfile + '.h5')
    
    nx    = f.attrs["nx"]             # Number of grid cells along X
    ny    = f.attrs["ny"]             # Number of grid cells along Y
    nz    = f.attrs["nz"]             # Number of grid cells along Z
    xmin  = f.attrs["xmin"]           # Simulation domain settings [m]
    xmax  = f.attrs["xmax"]           # [m]
    ymin  = f.attrs["ymin"]
    ymax  = f.attrs["ymax"]
    zmin  = f.attrs["zmin"]
    zmax  = f.attrs["zmax"]
    
    x_axis = (np.linspace(xmin, xmax, nx, endpoint=False) +
         (np.abs(xmin) + np.abs(xmax))/(2*float(nx)))

    y_axis = (np.linspace(ymin, ymax, ny, endpoint=False) +
             (np.abs(ymin) + np.abs(ymax))/(2*float(ny)))
    
    z_axis = (np.linspace(zmin, zmax, nz, endpoint=False) +
             (np.abs(zmin) + np.abs(zmax))/(2*float(nz)))
    
    
    s = []
    for keys in f.attrs.keys():
        
        if keys.startswith('s') and 'name' in keys:
            
            s.append(f.attrs[keys].decode('UTF-8'))
    try:
        
        key_species_name = "rho" + "%02d"%(s.index(key_species) + 1)
        
    except ValueError:
        
        print(f'\nThe element ({key_species}) you wish to include could not be found in the Amitis file.' +
              ' Please specify the names correctly. PROMETHEUS exits now.')
        
        sys.exit()

    rho = f[key_species_name][:]          # total charge density [C/m^3]
    den = rho/1.6e-19                       # total number density [#/m^3]
    
    pol_func = RegularGridInterpolator((x_axis, y_axis, z_axis), den, fill_value=0)
    
    return pol_func