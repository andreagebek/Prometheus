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



def plasmagrid(fundamentalsDict, architectureDict, scenarioDict, speciesDict, gridsDict):
    amitisfile = scenarioDict['amitisFilename']
    f = h5py.File(amitisfile + '.h5')
    