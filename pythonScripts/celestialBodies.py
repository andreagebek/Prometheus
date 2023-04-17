"""
This file has properties and methods for the
celestial bodies (star, planet, moon) involved
in the calculation.
Created on 16. September 2022 by Andrea Gebek.
"""

import pythonScripts.constants as const
import pythonScripts.geometryHandler as geom
import numpy as np
import shutil
import urllib.request as request
from contextlib import closing
import astropy.io.fits as fits
import os
from scipy.interpolate import interp1d



class Star:
    def __init__(self, R, M, T_eff, log_g, Z, alpha): # Stellar radius (cm), Stellar mass (g), Stellar effective temperature (K), Stellar surface gravity (log10(cm/s^2)), Metallicity [Fe/H], Alpha-enhancement [alpha/Fe]
        self.R = R
        self.M = M
        self.T_eff = T_eff
        self.log_g = log_g
        self.Z = Z
        self.alpha = alpha
        self.CLV_u1 = 0.
        self.CLV_u2 = 0.
        self.vsiniStarrot = 0.
        self.phiStarrot = 0.

    def addCLVparameters(self, CLV_u1, CLV_u2):
        self.CLV_u1 = CLV_u1
        self.CLV_u2 = CLV_u2
    
    def addRMparameters(self, vsiniStarrot, phiStarrot):
        self.vsiniStarrot = vsiniStarrot
        self.phiStarrot = phiStarrot

    def getSurfaceVelocity(self, phi, rho): # vsiniStarrot is the maximum velocity of the stellar rotation at the edge of the stellar disk, 
        # phiStarrot is the angle at which the maximum velocity (towards the observer) is obtained (measured in the same coordinate system as phi)

        v_los = self.vsiniStarrot * rho / self.R * np.cos(phi - self.phiStarrot)

        return v_los

    def round_to_grid(grid, value):

        diff = np.subtract(value, grid)
        arg = np.argmin(np.abs(diff))

        return grid[arg]

    def getSpectrum(self):
        """Querying a PHOENIX photosphere model, either from disk or from the
            PHOENIX website if the spectrum hasn't been downloaded before, and
            returning the names of the files in the data subfolder. These may then
            be read by the wrapper function read_spectrum() below.

            Parameters
            ----------
            T : int, float
                The model photosphere temperature. Acceptable values are:
                2300 - 7000 in steps of 100, and 7200 - 12000 in steps of 200.
            logg : int, float
                The model log(g) value. Acceptable values are:
                0.0 - 6.0 in steps of 0.5.
            Z : int, float
                The model metallicity [Fe/H] value. Acceptable values are:
                -4.0, -3.0, -2.0 and -1.5 to +1.0 in steps of 0.5.
            a : int, float
                The model alpha element enhancement [alpha/M]. Acceptable values are:
                -0.2 to 1.2 in steps of 0.2, but only for Fe/H of -3.0 to 0.0.
            Returns
            -------
            wlname, specname: str, str
                The names of the files containing the wavelength and flux axes of
                the requested spectrum.
            """

        #These contain the acceptable values.
        T_grid = np.concatenate((np.arange(2300,7100,100),np.arange(7200,12200,200)))
        log_g_grid = np.arange(0,6.5,0.5)
        Z_grid = np.concatenate((np.arange(-4,-1,1),np.arange(-1.5,1.5,0.5)))
        alpha_grid = np.arange(0,1.6,0.2)-0.2
        
        T_a = Star.round_to_grid(T_grid, self.T_eff)
        log_g_a = Star.round_to_grid(log_g_grid, self.log_g)
        Z_a = Star.round_to_grid(Z_grid, self.Z)
        alpha_a = Star.round_to_grid(alpha_grid, self.alpha)


        #This is where phoenix spectra are located.
        root = 'ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/'

        #We assemble a combination of strings to parse the user input into the URL,
        z_string = '{:.1f}'.format(float(Z_a))
        if Z_a > 0:
            z_string = '+' + z_string
        elif Z_a == 0:
            z_string = '-' + z_string
        else:
            z_string = z_string
        a_string=''
        if alpha_a > 0:
            a_string ='.Alpha=+'+'{:.2f}'.format(float(alpha_a))
        if alpha_a < 0:
            a_string ='.Alpha='+'{:.2f}'.format(float(alpha_a))
        t_string = str(int(T_a))
        if T_a < 10000:
            t_string = '0'+t_string
        g_string = '-'+'{:.2f}'.format(float(log_g_a))


        #These are URLS for the input files.
        waveurl = root+'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'
        specurl = root+'PHOENIX-ACES-AGSS-COND-2011/Z'+z_string+a_string+'/lte'+t_string+g_string+z_string+a_string+'.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'

        #These are the output filenames, they will also be returned so that the wrapper
        #of this function can take them in.
        wavename = 'WAVE.fits'

        specname = 'lte'+t_string+g_string+z_string+a_string+'.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'

        #Download PHOENIX spectra:

        with closing(request.urlopen(waveurl)) as r:
            with open(wavename, 'wb') as f:
                shutil.copyfileobj(r, f)

        with closing(request.urlopen(specurl)) as r:
            with open(specname, 'wb') as f:
                shutil.copyfileobj(r, f)

        F = fits.getdata(specname)
        w = fits.getdata(wavename)

        os.remove(wavename)
        os.remove(specname)

        return (w * 1e-8, F / np.pi) # Conversion to cgs-units. Note that Jens divides F by
        # a seemingly random factor of pi, but this should not bother the transit calculations here. 

    def calculateCLV(self, rho):
        # Quadratic limb darkening law

        arg = 1. - np.sqrt(1. - rho**2 / self.R**2)

        return 1. - self.CLV_u1 * arg - self.CLV_u2 * arg**2

    def calculateRM(self, phi, rho, wavelength):

        v_los = self.getSurfaceVelocity(phi, rho)
        shift = const.calculateDopplerShift(v_los)

        F_shifted = 10.**self.Fstar_function(wavelength / shift) # This contains the shifted flux incident on the exoplanet, instead of shifting the spectra
        # at all positions on the stellar surface just read in the spectrum at different wavelengths depending on the position in the spatial grid (phi and rho)

        return F_shifted

    def addFstarFunction(self, wavelength):

        PHOENIX_output = self.getSpectrum()
    
        w_star = PHOENIX_output[0]
        w_max = np.max(wavelength) * const.calculateDopplerShift(-self.vsiniStarrot)
        w_min = np.min(wavelength) * const.calculateDopplerShift(self.vsiniStarrot)

        SEL = (w_star >= w_min) * (w_star <= w_max)
        minArg = max(min(np.argwhere(SEL)).item() - 1, 0)
        maxArg = max(np.argwhere(SEL)).item() + 2 # Add one value at the beginning and end of wavelength array for safety

        w_starSEL = w_star[minArg:maxArg]
        F_0 = PHOENIX_output[1][minArg:maxArg]

        Fstar_function = interp1d(w_starSEL, np.log10(F_0), kind = 'linear')

        self.Fstar_function = Fstar_function


    def getFstarIntegrated(self, wavelength, grid):
        # Return F_star(wavelength, phi, rho) integrated over the stellar disk and over the upper part of the stellar disk (if rho_border is not the stellar radius)
        if self.vsiniStarrot == 0.: # No RM-effect

            FstarIntegrated = np.pi * self.R**2 * (1. - self.CLV_u1 / 3. - self.CLV_u2 / 6.) * np.ones_like(wavelength)

            upperTerm = 0.5 * (-self.CLV_u2 * self.R**2 - self.CLV_u1 * self.R**2 + self.R**2)
            term1 = -4. * self.R**2 * self.CLV_u1 * (1. - grid.rho_border**2 / self.R**2)**1.5
            term2 = self.R**2 * self.CLV_u2 * (6 * grid.rho_border**2 / self.R**2 + 8. * (1. - grid.rho_border**2 / self.R**2)**1.5 - 3. * (self.R**2 - grid.rho_border**2)**2 / self.R**4)
            lowerTerm = 1. / 12. * (term1 - term2 - 6. * self.CLV_u1 * grid.rho_border**2 + 6. * grid.rho_border**2)

            FstarUpper = 2. * np.pi * (upperTerm - lowerTerm) * np.ones_like(wavelength)
        

        else: # With RM-effect

            phiArray = grid.constructPhiAxis()
            delta_phi = grid.getDeltaPhi()
            rhoArray = grid.constructRhoAxis()
            delta_rho = grid.getDeltaRho()

            FstarIntegrated = np.zeros_like(wavelength)

            for phi in phiArray:
                for rho in rhoArray:
            
                    Fstar = self.calculateRM(phi, rho, wavelength)

                    Fstar *= self.calculateCLV(rho)
                    
                    FstarIntegrated += Fstar * delta_phi * delta_rho * rho

                    FstarUpper = np.zeros_like(wavelength) # To calculate this with the RM-effect we would need the spatial grid to run over the entire stellar disk

        return FstarIntegrated, FstarUpper

    def getFstar(self, phi, rho, wavelength):

        # Calculate F_star(wavelength, phi, rho)

        if self.vsiniStarrot == 0.: # No RM-effect
            Fstar = np.ones_like(wavelength) * self.calculateCLV(rho)

        else: # RM_effect = True

            Fstar = self.calculateRM(phi, rho, wavelength)

            Fstar *= self.calculateCLV(rho)

        return Fstar




class Planet:
    def __init__(self, name, R, M, a, hostStar): #  Reference radius (cm), Planetary mass(g), Orbital distance (cm)
        self.name = name
        self.R = R
        self.M = M
        self.a = a
        self.hostStar = hostStar

    def getPosition(self, orbphase):
        x_p = self.a * np.cos(orbphase)
        y_p = self.a * np.sin(orbphase)
        return x_p, y_p

    def getLOSvelocity(self, orbphase):
        v_los = -np.sin(orbphase) * np.sqrt(const.G * self.hostStar.M / self.a)
        return v_los

    def getDistanceFromPlanet(self, x, phi, rho, orbphase):

        y, z = geom.Grid.getCartesianFromCylinder(phi, rho)        
        x_p, y_p = self.getPosition(orbphase)
        r_fromPlanet = np.sqrt((x - x_p)**2 + (y - y_p)**2 + z**2)

        return r_fromPlanet

    def getTorusCoords(self, x, phi, rho, orbphase): # Get coordinates for a circumplanetary torus coordinate system

        y, z = geom.Grid.getCartesianFromCylinder(phi, rho)        
        x_p, y_p = self.getPosition(orbphase)
        a = np.sqrt((x - x_p)**2 + (y - y_p)**2)

        return a, z


class Moon:
    def __init__(self, midTransitOrbphase, R, a, hostPlanet):
        self.midTransitOrbphase = midTransitOrbphase
        self.R = R
        self.a = a
        self.hostPlanet = hostPlanet

    def getOrbphase(self, orbphase):

        orbphase_moon = self.midTransitOrbphase + orbphase * np.sqrt((self.hostPlanet.a**3 * self.hostPlanet.M) / (self.a**3 * self.hostPlanet.hostStar.M))

        return orbphase_moon

    def getPosition(self, orbphase):

        orbphase_moon = self.getOrbphase(orbphase)

        x_p, y_p = self.hostPlanet.getPosition(orbphase)

        x_moon = x_p + self.a * np.cos(orbphase_moon)
        y_moon = y_p + self.a * np.sin(orbphase_moon)

        return x_moon, y_moon

    def getLOSvelocity(self, orbphase):

        v_los_planet = self.hostPlanet.getLOSvelocity(orbphase)
        orbphase_moon = self.getOrbphase(orbphase)
        v_los = v_los_planet - np.sin(orbphase_moon) * np.sqrt(const.G * self.hostPlanet.M / self.a)

        return v_los

    def getDistanceFromMoon(self, x, phi, rho, orbphase):

        y, z = geom.Grid.getCartesianFromCylinder(phi, rho)
        x_moon, y_moon = self.getPosition(orbphase)        
        r_fromMoon = np.sqrt((x - x_moon)**2 + (y - y_moon)**2 + z**2)

        return r_fromMoon


class AvailablePlanets:
    def __init__(self):
        WASP49a = Star(1.038 * const.R_sun, 1.003 * const.M_sun,  5600., 4.5, -0.08, 0.) # Wyttenbach et al. 2017, Metallicity from Sousa et al. 2018, Alpha-enhancement unknown
        WASP49b = Planet('WASP-49b', 1.198 * const.R_J, 0.399 * const.M_J, 0.03873 * const.AU, hostStar = WASP49a) # Wyttenbach et al. 2017
        HD189733a = Star(0.756 * const.R_sun, 0.823 * const.M_sun, 5201., 4.64, -0.02, 0.) # Wyttenbach et al. 2015, T_eff, log_g, and Metallicity from Chavero et al. 2019, Alpha-enhancement unknown
        HD189733b = Planet('HD189733b', 1.138 * const.R_J, 1.138 * const.M_J, 0.0312 * const.AU, hostStar = HD189733a) # Wyttenbach et al. 2015
        Cancri55a = Star(0.98 * const.R_sun, 1.015 * const.M_sun, 5172., 4.43, 0.35, 0.) # Crida et al. 2018, a_p, T_eff, log_g, and Metallcity from Bourrier et al. 2018, Alpha-enhancement unknown) 
        Cancri55e = Planet('55-Cancri-e', 0.1737 * const.R_J, 0.02703 * const.M_J, 0.01544 * const.AU, hostStar = Cancri55a) # Crida et al. 2018
        WASP39a = Star(0.895 * const.R_sun, 0.93 * const.M_sun,  5400., 4.503, -0.12, 0.) # Faedi et al. 2018, Alpha-enhancement unknown
        WASP39b = Planet('WASP-39b', 1.27 * const.R_J, 0.28 * const.M_J, 0.0486 * const.AU, hostStar = WASP39a) # Faedi et al. 2018

        self.planetList = [WASP49b, HD189733b, Cancri55e, WASP39b]

    def listPlanetNames(self):
        planetNames = []
        for planet in self.planetList:
            planetNames.append(planet.name)
        return planetNames

    def findPlanet(self, namePlanet):
        for planet in self.planetList:
            if planet.name == namePlanet:
                return planet

        print('System', namePlanet, 'was not found.')
        return None
