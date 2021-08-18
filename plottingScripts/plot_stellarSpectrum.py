"""
Author: Andrea Gebek
Created on 18.8.2021
Plot the Phoenix input 
spectrum from the host star.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import json
import sys
sys.path.insert(0, '../pythonScripts') # Import from pythonScripts folder
from stellarspectrum import *
from constants import *

matplotlib.rcParams['axes.linewidth'] = 2.5
matplotlib.rcParams['xtick.major.size'] = 10
matplotlib.rcParams['xtick.minor.size'] = 6
matplotlib.rcParams['xtick.major.width'] = 2.5
matplotlib.rcParams['xtick.minor.width'] = 1.5
matplotlib.rcParams['ytick.major.size'] = 10
matplotlib.rcParams['ytick.minor.size'] = 6
matplotlib.rcParams['ytick.major.width'] = 2.5
matplotlib.rcParams['ytick.minor.width'] = 1.5
matplotlib.rcParams['image.origin'] = 'lower'
matplotlib.rcParams.update({'font.size': 26, 'font.weight': 'bold'})

"""
Helper function for the Doppler shift
"""

def Doppler(v):
    # If v is positive, receiver and source are moving towards each other
    beta = v / c
    shift = np.sqrt((1. - beta) / (1. + beta))

    return shift

"""
Read in settings file and calculate stellar spectrum
"""

paramsFilename = sys.argv[1]

with open('../../' + paramsFilename + '.txt') as file:
    param = json.load(file)

R_s = param['Architecture']['R_star']
period_starrot = param['Architecture']['period_starrot']

grids_dict = param['Grids']

wavelength = np.arange(grids_dict['lower_w'], grids_dict['upper_w'], grids_dict['resolution'])

T_eff = param['Architecture']['T_eff']
log_g = param['Architecture']['log_g']
Fe_H = param['Architecture']['Fe_H']
alpha_Fe = param['Architecture']['alpha_Fe']

PHOENIX_output = read_spectrum(T_eff, log_g, Fe_H, alpha_Fe)

w_star = PHOENIX_output[0]
v_max = 2. * np.pi * R_s / period_starrot
w_max = np.max(wavelength * Doppler(-v_max))
w_min = np.min(wavelength * Doppler(v_max))
SEL_w = np.argwhere((w_star > w_min) * (w_star < w_max))[:, 0]

SEL = np.concatenate((np.array([np.min(SEL_w) - 1]), SEL_w, np.array([np.max(SEL_w) + 1])))

w_star = w_star[SEL]
F_0 = PHOENIX_output[1][SEL]



"""
Plot the spectrum and store the figure
"""


fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

ax.plot(w_star * 1e8, F_0, color = 'blue', linewidth = 2)


ax.set_xlabel(r'$\lambda\,[\AA]$')
ax.set_ylabel(r'$F_{\nu}\,[\rm{erg\,s^{-1}\,cm^{-2}\,cm^{-1}}]$')

ax.minorticks_on()
ax.tick_params(which = 'both', direction = 'in', right = True, top = True)

ax.set_xlim(np.min(w_star * 1e8), np.max(w_star * 1e8))


plt.tight_layout()

plt.savefig('../../' + paramsFilename + '_stellarSpectrumPlot.pdf', dpi = 150)