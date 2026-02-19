"""
Unbinned likelihood analysis utilities for Fermi-LAT.

Provides functions to set up, configure, and fit unbinned likelihood
models using the Fermitools UnbinnedAnalysis framework.

Requires fermitools to be installed.
"""

import numpy as np


def setup_unbinned(evfile, scfile, expmap, ltcube, srcmdl, irfs):
    """
    Initialize an UnbinnedAnalysis likelihood object.

    Parameters
    ----------
    evfile : str
        Filtered + GTI event file.
    scfile : str
        Spacecraft (FT2) file.
    expmap : str
        Exposure map file.
    ltcube : str
        Livetime cube file.
    srcmdl : str
        Source model XML file.
    irfs : str
        Instrument response functions.

    Returns
    -------
    UnbinnedAnalysis.UnbinnedAnalysis
        Configured likelihood object.
    """
    from UnbinnedAnalysis import UnbinnedAnalysis, UnbinnedObs

    obs = UnbinnedObs(
        str(evfile), str(scfile),
        expMap=str(expmap),
        expCube=str(ltcube),
        irfs=irfs,
    )
    like = UnbinnedAnalysis(obs, str(srcmdl), optimizer='NEWMINUIT')
    return like


def fit(like, optimizer='NEWMINUIT', tol=1e-5, verbosity=0):
    """
    Run the likelihood fit.

    Parameters
    ----------
    like : UnbinnedAnalysis
        Likelihood object.
    optimizer : str
        Optimizer name (default 'NEWMINUIT').
    tol : float
        Convergence tolerance.
    verbosity : int
        Verbosity level (0=quiet).

    Returns
    -------
    float
        Final -log(likelihood) value.
    """
    like.tol = tol
    logL = like.fit(verbosity=verbosity, optimizer=optimizer)
    return logL


def freeze_all_params(like):
    """
    Freeze all parameters of all sources in the model.

    Parameters
    ----------
    like : UnbinnedAnalysis
        Likelihood object.
    """
    for src_name in like.sourceNames():
        src = like[src_name]
        for param_name in like.model[src_name].funcs['Spectrum'].paramNames:
            idx = like.par_index(src_name, param_name)
            like.freeze(idx)


def free_diffuse_norms(like, galdiff_name=None, iso_name=None):
    """
    Free only the normalization parameters of diffuse sources.

    Identifies diffuse sources by name pattern if names are not provided:
    - Galactic: contains 'gll_iem' or 'GAL'
    - Isotropic: contains 'iso' or 'ISO'

    Parameters
    ----------
    like : UnbinnedAnalysis
        Likelihood object.
    galdiff_name : str, optional
        Exact name of Galactic diffuse source. Auto-detected if None.
    iso_name : str, optional
        Exact name of isotropic diffuse source. Auto-detected if None.
    """
    for src_name in like.sourceNames():
        name_lower = src_name.lower()

        is_galdiff = (galdiff_name and src_name == galdiff_name) or \
                     ('gll_iem' in name_lower or 'gal' in name_lower)
        is_iso = (iso_name and src_name == iso_name) or \
                 ('iso' in name_lower)

        if is_galdiff:
            # Free Prefactor (normalization)
            try:
                idx = like.par_index(src_name, 'Prefactor')
                like.thaw(idx)
            except Exception:
                # Try 'Value' or 'Normalization'
                for pname in ['Value', 'Normalization', 'norm']:
                    try:
                        idx = like.par_index(src_name, pname)
                        like.thaw(idx)
                        break
                    except Exception:
                        continue

        if is_iso:
            # Free Normalization
            try:
                idx = like.par_index(src_name, 'Normalization')
                like.thaw(idx)
            except Exception:
                for pname in ['Prefactor', 'Value', 'norm']:
                    try:
                        idx = like.par_index(src_name, pname)
                        like.thaw(idx)
                        break
                    except Exception:
                        continue


def prepare_transient_pl2(like, src_name, index=-2.0,
                          emin=100.0, emax=100000.0):
    """
    Configure a PowerLaw2 transient source for profile likelihood scanning.

    Sets the spectral index to a fixed value and configures energy bounds.

    Parameters
    ----------
    like : UnbinnedAnalysis
        Likelihood object.
    src_name : str
        Name of the transient source in the model.
    index : float
        Fixed photon index (negative, e.g., -2.0).
    emin, emax : float
        Energy bounds [MeV].
    """
    # Set and freeze index
    idx_par = like.par_index(src_name, 'Index')
    like[idx_par] = index
    like.freeze(idx_par)

    # Set energy bounds
    try:
        ll_par = like.par_index(src_name, 'LowerLimit')
        like[ll_par] = emin
        like.freeze(ll_par)

        ul_par = like.par_index(src_name, 'UpperLimit')
        like[ul_par] = emax
        like.freeze(ul_par)
    except Exception:
        pass  # Not all models have LowerLimit/UpperLimit

    # Ensure flux (Integral) is free
    flux_par = like.par_index(src_name, 'Integral')
    like.thaw(flux_par)


def set_flux(like, src_name, flux):
    """
    Set the integral flux of a source and freeze it.

    Parameters
    ----------
    like : UnbinnedAnalysis
        Likelihood object.
    src_name : str
        Source name.
    flux : float
        Integral flux value [ph/cm2/s].
    """
    idx = like.par_index(src_name, 'Integral')

    # Ensure bounds accommodate the value
    par = like.model[src_name].funcs['Spectrum'].getParam('Integral')
    lo, hi = par.getBounds()
    if flux < lo:
        par.setBounds(flux * 0.1, hi)
    if flux > hi:
        par.setBounds(lo, flux * 10.0)

    like[idx] = flux
    like.freeze(idx)
