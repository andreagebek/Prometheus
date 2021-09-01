"""
This file stores number density profiles.
Created on 15. July 2021 by Andrea Gebek.
"""

from constants import *
from scipy.special import erf



def barometric(r, params):
    """Barometric law (to benchmark the numerical solution against analytical
    results)
    """
    T = params['T']
    P_0 = params['P_0']
    mu = params['mu']
    M_p = params['M_p']
    R_0 = params['R_0']

    n_0 = P_0 / (k_B * T)
    H = k_B * T * R_0**2 / (G * mu * M_p)
    n = n_0 * np.exp((R_0 - r) / H)
    return n 


def hydrostatic(r, params):
    """Hydrostatic equilibrium with constant temperature and
    mean molecular weight
    """
    T = params['T']
    P_0 = params['P_0']
    mu = params['mu']
    M_p = params['M_p']
    R_0 = params['R_0']

    n_0 = P_0 / (k_B * T)
    Jeans_0 = G * mu * M_p / (k_B * T * R_0)
    Jeans = G * mu * M_p / (k_B * T * r) * np.heaviside(r - R_0, 1) 
    # When the number density is evaluated inside the planet the Jeans parameter gets very large, leading to overflowing n
    n = n_0 * np.exp(Jeans - Jeans_0) 

    return n    



def escaping(r, params):
    """Power law with index q, for a hydrodynamically
    escaping atmosphere
    """
    q_esc = params['q_esc']
    R_0 = params['R_0']

    if 'T' in params.keys():
        T = params['T']
        P_0 = params['P_0']
        n_0_esc = P_0 / (k_B * T)
    
    else:
        n_0_esc = (q_esc - 3.) / (4. * np.pi * R_0**3)

    n = n_0_esc * (R_0 / r)**q_esc

    return n


def exomoon(r, params):
    """Scaled to observed number density profile at Io (Burger 2001)
    """
    R_moon = params['R_moon']

    n_0_exomoon = 0.34 / (4 * np.pi * R_moon**3)
    n = n_0_exomoon * (r / R_moon)**(-3.34)

    return n

def torus(a, z, params):
    """ Toroidal atmosphere (Johnson & Huggins 2006)
    """
    a_torus = params['a_torus']
    M_p = params['M_p']
    v_ej = params['v_ej']

    v_orbit = np.sqrt(G * M_p / a_torus)
    H_torus = a_torus * v_ej / v_orbit

    n_a = np.exp(- ((a - a_torus) / (4 * H_torus))**2)
    n_z = np.exp(- (z / H_torus)**2)

    term1 = 8. * H_torus**2 * np.exp(-a_torus**2 / (16 * H_torus**2))
    term2 = 2. * np.sqrt(np.pi) * a_torus * H_torus * (erf(a_torus / (4 * H_torus)) + 1.)

    n_0_torus = 1. / (2. * np.pi**1.5 * H_torus * (term1 + term2))

    n = n_0_torus * np.multiply(n_a, n_z) 
    
    return n


