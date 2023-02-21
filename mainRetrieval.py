# coding=utf-8
"""
Main code to run the radiative transfer calculation in inverse modelling.
Created on 21. February 2023 by Andrea Gebek.
"""

import prometheusScripts.geometryHandler as geom
import prometheusScripts.celestialBodies as bodies
import prometheusScripts.gasProperties as gasprop
import prometheusScripts.constants as const
import numpy as np
import matplotlib.pyplot as plt



"""
Free parameters
"""

midTransitOrbphase = 0.5 * np.pi # Orbital phase of the moon at planet mid-transit (in rad)
a_moon = 3. # Orbital radius of the moon (in planetary radii)
N_Na = 1e34 # Number of sodium atoms in the exosphere of the moon

"""
Set up the forward model
"""

wavelengthGrid = gasprop.WavelengthGrid(5890e-8, 5900e-8, 3e-8, 5e-9, 5e-10) # Wavelength grid specifications in cm: lower border, upper border, width of high-res. region (centered on the Na D line peaks), 
#resolution for the low-res. region, resolution for the high-res. region. The wavelength grid looks like this: |  |  |||  |  |||  |  |

W49b = bodies.AvailablePlanets().findPlanet('WASP-49b')
spatialGrid = geom.Grid(W49b.a, 5. * W49b.R, 20, W49b.hostStar.R, 35, 35, 0.25, 25) # Spatial grid specifications: Midpoint of grid in x-direction (corresponding to the orbital radius of
#WASP-49b), extent of the grid in x-direction, number of bins for the x-axis, upper boundary for the radial grid (corresponding to the host star radius), number of bins for the radial axis,
#number of bins for the angular axis, boundary for the orbital phase grid (in rad, i.e. np.pi would correspont to a full orbit, from -pi to +pi), number of bins for the orbital phase axis
W49b.hostStar.addCLVparameters(0.34, 0.28) # Add center-to-limb variations with parameters from Wyttenbach+2017
# Could also add RM-effect here...


moon = bodies.Moon(midTransitOrbphase, const.R_Io, a_moon * W49b.R, W49b) # Setting up to moon with an Io radius, and an unknown orbital radius / orbital phase

moonExosphere = gasprop.MoonExosphere(N_Na, 3.34, moon) # Setting up an exosphere with a power law, with an exponent of -3.34 and an unknown amount of Na atoms
moonExosphere.addConstituent('NaI', 1e6) # Velocity dispersion of 10 km/s
moonExosphere.constituents[0].addLookupFunctionToConstituent(wavelengthGrid)
atmos = gasprop.Atmosphere([moonExosphere], True)

main = gasprop.Transit(atmos, wavelengthGrid, spatialGrid)
main.addWavelength()


"""
Execute forward model
"""

R = main.sumOverChords()

wavelength = wavelengthGrid.constructWavelengthGrid([moonExosphere])


"""
Extract lightcurve
"""

filterBandwidth = 2e-9 # In cm. This gives the width of the wavelength interval over which we average the transit depth signal
NaD2_center = 5891.583253e-8 # In cm
NaD1_center = 5897.558147e-8 # In cm

SEL_filter = ((wavelength >= NaD2_center - filterBandwidth / 2.) * (wavelength <= NaD2_center + filterBandwidth / 2.)) + ((wavelength >= NaD1_center - filterBandwidth / 2.) * (wavelength <= NaD1_center + filterBandwidth / 2.))


LIGHTCURVE = np.mean(R[:, SEL_filter], axis = 1) # Lightcurve, not normalized
ORBPHASE = spatialGrid.constructOrbphaseAxis() # Orbital phase, in rad


"""
Plotting for quick testing
"""

"""
fig = plt.figure(figsize = (10., 10.))
ax = fig.add_subplot(111)

ax.plot(ORBPHASE, LIGHTCURVE)
ax.set_xlabel('Orbital phase [rad]')
ax.set_ylabel('Transit depth')

plt.tight_layout()
plt.savefig('lightcurve.png', dpi = 150.)
"""
