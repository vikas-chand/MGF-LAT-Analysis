"""
Flux, energy, and spectral conversions for Fermi-LAT analysis.

All functions accept energies in MeV and fluxes in ph/cm2/s unless noted.
"""

import numpy as np
from ..constants import MeV2erg, Mpc_cm


# ─── Power-law integration helpers ──────────────────────────────────────────

def _pow_int(E1, E2, gamma):
    """
    Integral of E^gamma from E1 to E2.

    Parameters
    ----------
    E1, E2 : float
        Integration bounds [MeV].
    gamma : float
        Spectral index (e.g., -2.0).

    Returns
    -------
    float
        Value of the integral.
    """
    if abs(gamma + 1.0) < 1e-10:
        return np.log(E2 / E1)
    gp1 = gamma + 1.0
    return (E2**gp1 - E1**gp1) / gp1


def K_from_photon_flux(Fph, gamma, emin, emax):
    """
    Convert integral photon flux to spectral normalization K.

    For dN/dE = K * E^gamma, the photon flux is:
        Fph = integral(K * E^gamma, emin, emax)

    Parameters
    ----------
    Fph : float
        Integral photon flux [ph/cm2/s].
    gamma : float
        Photon index (negative, e.g., -2.0).
    emin, emax : float
        Energy bounds [MeV].

    Returns
    -------
    float
        Normalization K [ph/cm2/s/MeV].
    """
    return Fph / _pow_int(emin, emax, gamma)


def energy_flux_from_K(K, gamma, emin, emax):
    """
    Convert spectral normalization K to energy flux.

    Parameters
    ----------
    K : float
        Spectral normalization [ph/cm2/s/MeV].
    gamma : float
        Photon index.
    emin, emax : float
        Energy bounds [MeV].

    Returns
    -------
    float
        Energy flux [erg/cm2/s].
    """
    # Energy flux = integral(K * E^gamma * E * MeV2erg, emin, emax)
    return K * _pow_int(emin, emax, gamma + 1.0) * MeV2erg


def photon_flux_from_K(K, gamma, emin, emax):
    """
    Convert spectral normalization K to integral photon flux.

    Parameters
    ----------
    K : float
        Spectral normalization [ph/cm2/s/MeV].
    gamma : float
        Photon index.
    emin, emax : float
        Energy bounds [MeV].

    Returns
    -------
    float
        Photon flux [ph/cm2/s].
    """
    return K * _pow_int(emin, emax, gamma)


def convert_flux_ul(Fph_ul, gamma=-2.0, band_in=(100, 100000),
                    band_out=(100, 10000)):
    """
    Convert a photon-flux upper limit from one energy band to another,
    assuming a power-law spectrum with index gamma.

    Parameters
    ----------
    Fph_ul : float
        Upper limit photon flux [ph/cm2/s] in band_in.
    gamma : float
        Assumed photon index.
    band_in : tuple
        (emin, emax) in MeV for the input band.
    band_out : tuple
        (emin, emax) in MeV for the output band.

    Returns
    -------
    float
        Photon flux upper limit [ph/cm2/s] in band_out.
    """
    K = K_from_photon_flux(Fph_ul, gamma, band_in[0], band_in[1])
    return photon_flux_from_K(K, gamma, band_out[0], band_out[1])


def e2dnde_from_band_flux(Fband, emin, emax, gamma=-2.0):
    """
    Convert band photon flux to E^2 dN/dE at the log-center energy.

    Parameters
    ----------
    Fband : float
        Integral photon flux in the band [ph/cm2/s].
    emin, emax : float
        Band energy bounds [MeV].
    gamma : float
        Assumed photon index.

    Returns
    -------
    tuple
        (Eref, e2dnde) where Eref is log-center energy [MeV]
        and e2dnde is E^2 dN/dE [MeV/cm2/s].
    """
    Eref = np.sqrt(emin * emax)
    K = K_from_photon_flux(Fband, gamma, emin, emax)
    dnde_at_Eref = K * Eref**gamma
    e2dnde = Eref**2 * dnde_at_Eref
    return Eref, e2dnde


def isotropic_energy(d_Mpc, energy_flux_erg, dt_s):
    """
    Compute isotropic-equivalent energy from energy flux.

    Parameters
    ----------
    d_Mpc : float
        Luminosity distance [Mpc].
    energy_flux_erg : float
        Energy flux [erg/cm2/s].
    dt_s : float
        Effective observation time [s].

    Returns
    -------
    float
        Isotropic-equivalent energy [erg].
    """
    d_cm = d_Mpc * Mpc_cm
    return 4.0 * np.pi * d_cm**2 * energy_flux_erg * dt_s


def photon_flux_to_energy_flux(Fph, gamma=-2.0, emin=100.0, emax=100000.0):
    """
    Convert photon flux to energy flux.

    Parameters
    ----------
    Fph : float
        Photon flux [ph/cm2/s].
    gamma : float
        Photon index.
    emin, emax : float
        Energy bounds [MeV].

    Returns
    -------
    float
        Energy flux [erg/cm2/s].
    """
    K = K_from_photon_flux(Fph, gamma, emin, emax)
    return energy_flux_from_K(K, gamma, emin, emax)
