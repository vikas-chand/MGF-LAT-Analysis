#!/usr/bin/env python3
"""
Main analysis client for fermi_lat_transients.

Usage
-----
    python -m fermi_lat_transients.client inputs.yaml

Or via the console entry point (after pip install):
    lat-transients-run inputs.yaml

This script orchestrates the full analysis pipeline:
1. Load configuration and source catalog
2. For each source: run Fermitools pipeline, likelihood fit, profile scan
3. Stack profiles across sources
4. Compute combined upper limits and SEDs
"""

import sys
import logging
import numpy as np
from pathlib import Path

from .config import load_config
from .data.targets import load_catalog
from .data.io import save_profile
from .pipeline.gtapps import run_pipeline
from .pipeline.model_builder import build_source_model, add_transient_source
from .likelihood.unbinned import (
    setup_unbinned, fit, freeze_all_params,
    free_diffuse_norms, prepare_transient_pl2,
)
from .likelihood.profile import profile_scan, upper_limit_from_profile
from .likelihood.sed import sed_profile_scan, compute_sed_points
from .stacking.stack import stack_profiles, stacked_upper_limit
from .stacking.sed_stack import stack_sed_profiles, stacked_sed_points
from .plotting.profile_plots import plot_profile, plot_stacked_profile
from .plotting.sed_plots import plot_sed, plot_stacked_sed

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
)
log = logging.getLogger(__name__)


def analyze_single_source(target, cfg, work_dir=None):
    """
    Run the full analysis pipeline for a single source.

    Parameters
    ----------
    target : dict
        Target dictionary from data.targets.
    cfg : Config
        Configuration object.
    work_dir : Path, optional
        Output directory for this source.

    Returns
    -------
    dict
        Analysis results including profile scan and optional SED.
    """
    name = target['name']
    if work_dir is None:
        work_dir = Path(cfg.output_dir) / name
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    log.info(f"=== Analyzing {name} ===")
    log.info(f"  RA={target['ra']:.4f}, Dec={target['dec']:.4f}")
    log.info(f"  Time: {target['tmin']:.1f} - {target['tmax']:.1f} MET")

    # Step 1: Fermitools pipeline
    log.info("  Running Fermitools pipeline...")
    files = run_pipeline(target, cfg, work_dir=work_dir)

    # Step 2: Build source model
    log.info("  Building source model...")
    model_xml = str(work_dir / f"{name}_model.xml")
    build_source_model(
        catalog_fits=cfg.catalog_fits,
        ra=target['ra'], dec=target['dec'], roi=cfg.roi,
        galdiff=cfg.galdiff, isodiff=cfg.isodiff,
        output_xml=model_xml,
    )

    # Add transient source
    model_with_src = str(work_dir / f"{name}_model_with_src.xml")
    add_transient_source(
        model_xml, name=cfg.transient_src_name,
        ra=target['ra'], dec=target['dec'],
        spectrum=cfg.transient_spectrum,
        emin=cfg.emin, emax=cfg.emax,
        index=cfg.spectral_index,
        output_xml=model_with_src,
    )

    # Step 3: Diffuse responses
    log.info("  Computing diffuse responses...")
    from .pipeline.gtapps import gtdiffrsp
    gtdiffrsp(files['evfile_gti'], cfg.ft2, model_with_src, cfg.irfs)

    # Step 4: Likelihood setup and fit
    log.info("  Setting up likelihood...")
    like = setup_unbinned(
        files['evfile_gti'], cfg.ft2,
        files['expmap'], files['ltcube'],
        model_with_src, cfg.irfs,
    )

    freeze_all_params(like)
    free_diffuse_norms(like)
    prepare_transient_pl2(
        like, cfg.transient_src_name,
        index=cfg.spectral_index,
        emin=cfg.emin, emax=cfg.emax,
    )
    fit(like)

    # Step 5: Profile likelihood scan
    log.info("  Running profile likelihood scan...")
    flux_grid = np.concatenate((
        [0.0],
        np.logspace(
            np.log10(cfg.flux_min),
            np.log10(cfg.flux_max),
            cfg.flux_npts,
        ),
    ))
    prof = profile_scan(like, cfg.transient_src_name, flux_grid=flux_grid)

    save_profile(
        str(work_dir / f"{name}_profile.npz"),
        prof['flux_grid'], prof['logL'],
        logL_null=prof['logL_null'],
        test_statistic=prof['ts'],
    )

    ul = upper_limit_from_profile(
        prof['flux_grid'], prof['ts'], delta_ts=cfg.delta_ts_95,
    )

    log.info(f"  TS_max = {prof['ts_max']:.2f}")
    log.info(f"  95%% CL UL = {ul['ul95']:.3e} ph/cm2/s")

    plot_profile(
        prof['flux_grid'], prof['ts'], ul_result=ul,
        title=name,
        save_path=str(work_dir / f"{name}_profile.png"),
    )

    result = {
        'name': name,
        'target': target,
        'files': files,
        'profile': prof,
        'upper_limit': ul,
    }

    # Optional: SED
    if cfg.get('sed_bins'):
        log.info("  Running SED profile scan...")
        sed = sed_profile_scan(
            like, cfg.transient_src_name,
            energy_bins=cfg.sed_bins,
            flux_grid=flux_grid,
            index=cfg.spectral_index,
        )
        sed_pts = compute_sed_points(sed)
        result['sed'] = sed
        result['sed_points'] = sed_pts

        plot_sed(
            sed_pts, title=f"{name} SED",
            save_path=str(work_dir / f"{name}_sed.png"),
        )

    return result


def run_analysis(config_path):
    """
    Run the complete analysis pipeline from a YAML config file.

    Parameters
    ----------
    config_path : str
        Path to inputs.yaml configuration file.
    """
    cfg = load_config(config_path)
    cfg.validate()

    log.info(f"Configuration loaded from {config_path}")

    targets = load_catalog(
        cfg.catalog_file,
        fmt=cfg.catalog_format,
        column_name=cfg.column_name,
        column_ra=cfg.column_ra,
        column_dec=cfg.column_dec,
        column_trigger=cfg.column_trigger_time,
        dt_pre=cfg.dt_pre,
        dt_post=cfg.dt_post,
    )
    log.info(f"Loaded {len(targets)} targets from {cfg.catalog_file}")

    # Analyze each source
    all_results = []
    for target in targets:
        try:
            result = analyze_single_source(target, cfg)
            all_results.append(result)
        except Exception as e:
            log.error(f"Failed on {target['name']}: {e}")
            continue

    if not all_results:
        log.warning("No sources analyzed successfully")
        return

    # Stacking
    log.info(f"\n=== Stacking {len(all_results)} sources ===")
    profile_list = [r['profile'] for r in all_results]
    stacked = stack_profiles(profile_list)
    stacked_ul = stacked_upper_limit(stacked, delta_ts=cfg.delta_ts_95)

    log.info(f"  Stacked TS_max = {stacked_ul['ts_max']:.2f}")
    log.info(f"  Stacked 95%% CL UL = {stacked_ul['ul95']:.3e} ph/cm2/s")

    out_dir = Path(cfg.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    plot_stacked_profile(
        stacked, ul_result=stacked_ul,
        title='Stacked Profile Likelihood',
        save_path=str(out_dir / 'stacked_profile.png'),
    )

    # SED stacking
    if cfg.get('sed_bins'):
        sed_list = [r['sed'] for r in all_results if 'sed' in r]
        if sed_list:
            stacked_sed_result = stack_sed_profiles(sed_list)
            stacked_sed_pts = stacked_sed_points(stacked_sed_result)

            individual_sed_pts = [
                r['sed_points'] for r in all_results if 'sed_points' in r
            ]
            source_names = [
                r['name'] for r in all_results if 'sed_points' in r
            ]

            plot_stacked_sed(
                individual_sed_pts, stacked_sed_pts,
                title='Stacked SED',
                source_names=source_names,
                save_path=str(out_dir / 'stacked_sed.png'),
            )

    log.info("\n=== Analysis complete ===")


def main():
    """Console entry point."""
    if len(sys.argv) < 2:
        print("Usage: lat-transients-run <inputs.yaml>")
        sys.exit(1)
    run_analysis(sys.argv[1])


if __name__ == '__main__':
    main()
