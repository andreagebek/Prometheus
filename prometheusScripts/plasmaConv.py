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

SCRIPTPATH = os.path.realpath(__file__)
GITPATH = os.path.dirname(os.path.dirname(SCRIPTPATH))
PARENTPATH = os.path.dirname(os.path.dirname(GITPATH))

import prometheusScripts.geometryHandler as geom
import prometheusScripts.fluxDecrease as flux



def plasmaGrid(phi, rho, orbphase, xArray, key_scenario, specificScenarioDict, architectureDict, fundamentalsDict):
    amitisfile = specificScenarioDict['amitisFilename']
    f = h5py.File(amitisfile + '.h5')
    
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
    
    rho_tot = f["rho02"][:] + f["rho03"][:]        # total charge density [C/m^3]
    den_tot = rho_tot/1.6e-19        # total number density [#/m^3]
    
    pol_func = RegularGridInterpolator((x_axis, y_axis, z_axis), den_tot, fill_value=0)