"""
Fireball scaling relations for magnetar giant flares.

Based on the scaling relations connecting prompt gamma-ray luminosity
to fireball parameters (peak kinetic luminosity, temperature, Lorentz factor).
"""

import numpy as np


def L0_from_Lgamma(Lgamma_47, xi_gamma=1.0):
    """
    Peak kinetic luminosity from prompt gamma luminosity.

    L0 ~ 3e47 * xi_{gamma,-0.5}^{-1} * L_{gamma,47} erg/s

    Parameters
    ----------
    Lgamma_47 : float or array
        Prompt gamma luminosity in units of 1e47 erg/s.
    xi_gamma : float
        Gamma-ray efficiency parameter xi_{gamma,-0.5} (default 1.0).

    Returns
    -------
    float or array
        Peak kinetic luminosity [erg/s].
    """
    return 3.0e47 * (xi_gamma**-1.0) * Lgamma_47


def T0_from_Lgamma(Lgamma_47, xi_gamma=1.0, r0_6=1.0):
    """
    Characteristic temperature from prompt gamma luminosity.

    T0 ~ 275 * xi_{gamma,-0.5}^{-1/4} * L_{gamma,47}^{1/4} * r_{0,6}^{-1/2} keV

    Parameters
    ----------
    Lgamma_47 : float or array
        Prompt gamma luminosity in units of 1e47 erg/s.
    xi_gamma : float
        Gamma-ray efficiency parameter.
    r0_6 : float
        Initial radius in units of 1e6 cm.

    Returns
    -------
    float or array
        Characteristic temperature [keV].
    """
    return 275.0 * (xi_gamma**-0.25) * (Lgamma_47**0.25) * (r0_6**-0.5)


def eta_from_Lgamma(Lgamma_47, xi_gamma=1.0, r0_6=1.0):
    """
    Bulk Lorentz factor (eta_*) from prompt gamma luminosity.

    eta_* ~ 140 * xi_{gamma,-0.5}^{-1/4} * L_{gamma,47}^{1/4} * r_{0,6}^{-1/4}

    Parameters
    ----------
    Lgamma_47 : float or array
        Prompt gamma luminosity in units of 1e47 erg/s.
    xi_gamma : float
        Gamma-ray efficiency parameter.
    r0_6 : float
        Initial radius in units of 1e6 cm.

    Returns
    -------
    float or array
        Bulk Lorentz factor eta_*.
    """
    return 140.0 * (xi_gamma**-0.25) * (Lgamma_47**0.25) * (r0_6**-0.25)


def Ekiso(Eg, xi_e):
    """
    Isotropic kinetic energy from gamma-ray energy and efficiency.

    E_k,iso = E_gamma * (1 - xi_e) / xi_e

    Parameters
    ----------
    Eg : float
        Gamma-ray energy [erg].
    xi_e : float
        Radiative efficiency (0 < xi_e < 1).

    Returns
    -------
    float
        Isotropic kinetic energy [erg].
    """
    return Eg * (1.0 - xi_e) / xi_e


def fireball_params(Lgamma_47, xi_gamma=1.0, r0_6=1.0):
    """
    Compute all fireball parameters at once.

    Parameters
    ----------
    Lgamma_47 : float or array
        Prompt gamma luminosity in units of 1e47 erg/s.
    xi_gamma : float
        Gamma-ray efficiency parameter.
    r0_6 : float
        Initial radius in units of 1e6 cm.

    Returns
    -------
    dict
        Dictionary with keys 'L0', 'T0_keV', 'eta_star'.
    """
    return {
        'L0': L0_from_Lgamma(Lgamma_47, xi_gamma),
        'T0_keV': T0_from_Lgamma(Lgamma_47, xi_gamma, r0_6),
        'eta_star': eta_from_Lgamma(Lgamma_47, xi_gamma, r0_6),
    }
