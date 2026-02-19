"""
Cosmological rate and volume calculations.

Provides functions for computing comoving volume elements and
observable volumes as a function of redshift, used for constraining
transient event rates.
"""

import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u


def comoving_volume_element(z_array, H0=70.0, Om0=0.3):
    """
    Compute the differential comoving volume element dV/dz per steradian.

    Parameters
    ----------
    z_array : array-like
        Redshift values.
    H0 : float
        Hubble constant [km/s/Mpc].
    Om0 : float
        Matter density parameter.

    Returns
    -------
    dVdz : ndarray
        Differential comoving volume [Mpc^3/sr] at each redshift.
    cosmo : FlatLambdaCDM
        The cosmology object used.
    """
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    z_array = np.asarray(z_array)

    dVdz = cosmo.differential_comoving_volume(z_array).to(u.Mpc**3 / u.sr).value
    return dVdz, cosmo


def observable_volume(z_array, H0=70.0, Om0=0.3):
    """
    Compute the cumulative observable volume as a function of redshift,
    accounting for time dilation (dividing by 1+z).

    Parameters
    ----------
    z_array : array-like
        Redshift values (should be finely sampled).
    H0 : float
        Hubble constant [km/s/Mpc].
    Om0 : float
        Matter density parameter.

    Returns
    -------
    dict
        'z': redshift array,
        'dVdz': differential volume [Mpc^3/sr],
        'dVdz_corrected': dV/dz / (1+z) [Mpc^3/sr],
        'V_cumulative': cumulative volume [Mpc^3/sr],
        'age_Gyr': age of universe at each z [Gyr].
    """
    z_array = np.asarray(z_array)
    dVdz, cosmo = comoving_volume_element(z_array, H0, Om0)

    # Time-dilation corrected integrand
    dVdz_corrected = dVdz / (1.0 + z_array)

    # Cumulative integral
    V_cumulative = np.zeros_like(z_array)
    for i in range(1, len(z_array)):
        V_cumulative[i] = np.trapz(dVdz_corrected[:i+1], z_array[:i+1])

    # Age of universe
    age_Gyr = cosmo.age(z_array).to(u.Gyr).value

    return {
        'z': z_array,
        'dVdz': dVdz,
        'dVdz_corrected': dVdz_corrected,
        'V_cumulative': V_cumulative,
        'age_Gyr': age_Gyr,
    }
