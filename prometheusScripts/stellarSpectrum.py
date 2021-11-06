"""
Get stellar input spectrum from the synthetic PHOENIX library.
Created on 9. August 2021 by Andrea Gebek (Initial code
from Jens Hoeijmakers' StarRotator code on GitHub,
https://github.com/Hoeijmakers/StarRotator/blob/master/lib/stellar_spectrum.py)
"""

import numpy as np
import requests
import shutil
import urllib.request as request
from contextlib import closing
import sys
import astropy.io.fits as fits
import os

def round_to_grid(grid, value):

    diff = np.subtract(value, grid)
    arg = np.argmin(np.abs(diff))

    return grid[arg]

def getSpectrum(T,logg,Z,a):
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


    #This is where phoenix spectra are located.
    root = 'ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/'

    #We assemble a combination of strings to parse the user input into the URL,
    z_string = '{:.1f}'.format(float(Z))
    if Z > 0:
        z_string = '+' + z_string
    elif Z == 0:
        z_string = '-' + z_string
    else:
        z_string = z_string
    a_string=''
    if a > 0:
        a_string ='.Alpha=+'+'{:.2f}'.format(float(a))
    if a < 0:
        a_string ='.Alpha='+'{:.2f}'.format(float(a))
    t_string = str(int(T))
    if T < 10000:
        t_string = '0'+t_string
    g_string = '-'+'{:.2f}'.format(float(logg))


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

    return(wavename,specname)

def readSpectrum(T,logg,metallicity=0.0,alpha=0.0):
    """Wrapper for the function get_spectrum() above, that checks that the input
        T, log(g), metallicity and alpha are in accordance with what is provided by
        PHOENIX (as of November 1st, 2019), and passes them on to get_spectrum().
        Parameters
        ----------
        T : int,float
            The model photosphere temperature. Acceptable values are:
            2300 - 7000 in steps of 100, and 7200 - 12000 in steps of 200.
        logg : int,float
            The model log(g) value. Acceptable values are:
            0.0 - 6.0 in steps of 0.5.
        metallicity : int,float (optional, default = 0)
            The model metallicity [Fe/H] value. Acceptable values are:
            -4.0, -3.0, -2.0 and -1.5 to +1.0 in steps of 0.5.
            If no location is given, the ``location`` attribute of the Time
            object is used.
        alpha : int,float (optional, default = 0)
            The model alpha element enhancement [alpha/M]. Acceptable values are:
            -0.2 to 1.2 in steps of 0.2, but only for Fe/H of -3.0 to 0.0.
        Returns
        -------
        w,f : np.array(),np.array()
            The wavelength (cm) and flux (erg/s/cm^2/cm) axes of the requested
            spectrum.
        """

    #These contain the acceptable values.
    T_a = np.concatenate((np.arange(2300,7100,100),np.arange(7200,12200,200)))
    logg_a = np.arange(0,6.5,0.5)
    FeH_a = np.concatenate((np.arange(-4,-1,1),np.arange(-1.5,1.5,0.5)))
    alpha_a = np.arange(0,1.6,0.2)-0.2

    T = round_to_grid(T_a, T)
    logg = round_to_grid(logg_a, logg)
    metallicity = round_to_grid(FeH_a, metallicity)
    alpha = round_to_grid(alpha_a, alpha)

    #Retrieve the spectra:
    wavename,specname = getSpectrum(T,logg,metallicity,alpha)
    f = fits.getdata(specname)
    w = fits.getdata(wavename)

    os.remove(wavename)
    os.remove(specname)

    return(w * 1e-8, f / np.pi) # Conversion to cgs-units. Note that Jens divides f by
    # a seemingly random factor of pi, but this should not bother the transit calculations here. 