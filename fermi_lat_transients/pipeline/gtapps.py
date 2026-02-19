"""
Thin wrappers around Fermi Science Tools (gt_apps).

Each function accepts explicit parameters (no hardcoded values)
and calls the corresponding Fermi tool.

These wrappers require fermitools to be installed in the environment.
"""

import os
from pathlib import Path


def gtselect(evfile, outfile, ra, dec, roi, tmin, tmax,
             emin=100.0, emax=100000.0, zmax=100.0,
             evclass=8, evtype=3):
    """
    Run gtselect to filter LAT events.

    Parameters
    ----------
    evfile : str
        Input event file (FT1) or @filelist.
    outfile : str
        Output filtered event file.
    ra, dec : float
        Source coordinates [degrees].
    roi : float
        Region of interest radius [degrees].
    tmin, tmax : float
        Time range [MET seconds].
    emin, emax : float
        Energy range [MeV].
    zmax : float
        Maximum zenith angle [degrees].
    evclass : int
        Event class (8=Transient, 128=Source).
    evtype : int
        Event type (3=FRONT+BACK).

    Returns
    -------
    str
        Path to the output file.
    """
    import gt_apps as my_apps

    my_apps.filter['evclass'] = evclass
    my_apps.filter['evtype'] = evtype
    my_apps.filter['ra'] = ra
    my_apps.filter['dec'] = dec
    my_apps.filter['rad'] = roi
    my_apps.filter['emin'] = emin
    my_apps.filter['emax'] = emax
    my_apps.filter['zmax'] = zmax
    my_apps.filter['tmin'] = tmin
    my_apps.filter['tmax'] = tmax
    my_apps.filter['infile'] = str(evfile)
    my_apps.filter['outfile'] = str(outfile)
    my_apps.filter.run()

    return str(outfile)


def gtmktime(scfile, evfile, outfile, ra=None, dec=None, roi=None,
             zmax=100.0, filter_expr=None):
    """
    Run gtmktime to apply good time interval selections.

    Parameters
    ----------
    scfile : str
        Spacecraft (FT2) file.
    evfile : str
        Input event file (from gtselect).
    outfile : str
        Output event file with GTI applied.
    ra, dec : float, optional
        Source coordinates for limb angle cut.
    roi : float, optional
        ROI radius for ANGSEP limb cut.
    zmax : float
        Maximum zenith angle for limb cut [degrees].
    filter_expr : str, optional
        Custom filter expression. If None, uses standard:
        (DATA_QUAL>0)&&(LAT_CONFIG==1)&&(LIVETIME>0)

    Returns
    -------
    str
        Path to the output file.
    """
    import gt_apps as my_apps

    if filter_expr is None:
        filter_expr = "(DATA_QUAL>0)&&(LAT_CONFIG==1)&&(LIVETIME>0)"

        # Add ANGSEP limb cut if coordinates provided
        if ra is not None and dec is not None and roi is not None:
            limb_angle = zmax - roi
            filter_expr += (
                f"&&(ANGSEP(RA_ZENITH,DEC_ZENITH,{ra},{dec})"
                f">({zmax}-{roi}))"
            )

    my_apps.maketime['scfile'] = str(scfile)
    my_apps.maketime['filter'] = filter_expr
    my_apps.maketime['roicut'] = 'no'
    my_apps.maketime['evfile'] = str(evfile)
    my_apps.maketime['outfile'] = str(outfile)
    my_apps.maketime.run()

    return str(outfile)


def gtltcube(scfile, evfile, outfile, zmax=100.0):
    """
    Run gtltcube to compute the livetime cube.

    Parameters
    ----------
    scfile : str
        Spacecraft (FT2) file.
    evfile : str
        Filtered event file.
    outfile : str
        Output livetime cube file.
    zmax : float
        Maximum zenith angle [degrees].

    Returns
    -------
    str
        Path to the output file.
    """
    import gt_apps as my_apps

    my_apps.expCube['evfile'] = str(evfile)
    my_apps.expCube['scfile'] = str(scfile)
    my_apps.expCube['outfile'] = str(outfile)
    my_apps.expCube['zmax'] = zmax
    my_apps.expCube['dcostheta'] = 0.025
    my_apps.expCube['binsz'] = 1
    my_apps.expCube.run()

    return str(outfile)


def gtexpmap(evfile, scfile, ltcube, outfile, irfs,
             roi=12.0, nlong=120, nlat=120, ebin=20):
    """
    Run gtexpmap to compute the exposure map.

    Parameters
    ----------
    evfile : str
        Filtered event file.
    scfile : str
        Spacecraft (FT2) file.
    ltcube : str
        Livetime cube file.
    outfile : str
        Output exposure map file.
    irfs : str
        Instrument response functions (e.g., 'P8R3_TRANSIENT020E_V3').
    roi : float
        ROI radius + margin [degrees].
    nlong, nlat : int
        Number of pixels.
    ebin : int
        Number of energy bins.

    Returns
    -------
    str
        Path to the output file.
    """
    import gt_apps as my_apps

    my_apps.expMap['evfile'] = str(evfile)
    my_apps.expMap['scfile'] = str(scfile)
    my_apps.expMap['expcube'] = str(ltcube)
    my_apps.expMap['outfile'] = str(outfile)
    my_apps.expMap['irfs'] = irfs
    my_apps.expMap['srcrad'] = roi + 10.0  # margin beyond ROI
    my_apps.expMap['nlong'] = nlong
    my_apps.expMap['nlat'] = nlat
    my_apps.expMap['nenergies'] = ebin
    my_apps.expMap.run()

    return str(outfile)


def gtdiffrsp(evfile, scfile, srcmdl, irfs):
    """
    Run gtdiffrsp to compute diffuse source responses.

    Parameters
    ----------
    evfile : str
        Filtered event file.
    scfile : str
        Spacecraft (FT2) file.
    srcmdl : str
        Source model XML file.
    irfs : str
        Instrument response functions.

    Returns
    -------
    str
        Path to the event file (modified in-place).
    """
    import gt_apps as my_apps

    my_apps.diffResps['evfile'] = str(evfile)
    my_apps.diffResps['scfile'] = str(scfile)
    my_apps.diffResps['srcmdl'] = str(srcmdl)
    my_apps.diffResps['irfs'] = irfs
    my_apps.diffResps.run()

    return str(evfile)


def run_pipeline(target, cfg, work_dir=None):
    """
    Run the full Fermitools pipeline for a single target.

    Steps: gtselect -> gtmktime -> gtltcube -> gtexpmap

    Parameters
    ----------
    target : dict
        Target dictionary from data.targets.build_target().
    cfg : Config
        Configuration object.
    work_dir : str or Path, optional
        Working directory for output files. Defaults to
        cfg.output_dir / target['name'].

    Returns
    -------
    dict
        Paths to all generated files:
        'evfile_filt', 'evfile_gti', 'ltcube', 'expmap'.
    """
    from pathlib import Path

    if work_dir is None:
        work_dir = Path(cfg.output_dir) / target['name']
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    prefix = target['name']

    # Step 1: gtselect
    evfile_filt = str(work_dir / f"{prefix}_filt.fits")
    gtselect(
        evfile=cfg.ft1,
        outfile=evfile_filt,
        ra=target['ra'], dec=target['dec'],
        roi=cfg.roi,
        tmin=target['tmin'], tmax=target['tmax'],
        emin=cfg.emin, emax=cfg.emax,
        zmax=cfg.zmax,
        evclass=cfg.evclass, evtype=cfg.evtype,
    )

    # Step 2: gtmktime
    evfile_gti = str(work_dir / f"{prefix}_gti.fits")
    gtmktime(
        scfile=cfg.ft2,
        evfile=evfile_filt,
        outfile=evfile_gti,
        ra=target['ra'], dec=target['dec'],
        roi=cfg.roi, zmax=cfg.zmax,
    )

    # Step 3: gtltcube (use pre-computed if available)
    if cfg.get('ltcube') is not None:
        ltcube_file = cfg.ltcube
    else:
        ltcube_file = str(work_dir / f"{prefix}_ltcube.fits")
        gtltcube(
            scfile=cfg.ft2,
            evfile=evfile_gti,
            outfile=ltcube_file,
            zmax=cfg.zmax,
        )

    # Step 4: gtexpmap
    expmap_file = str(work_dir / f"{prefix}_expmap.fits")
    gtexpmap(
        evfile=evfile_gti,
        scfile=cfg.ft2,
        ltcube=ltcube_file,
        outfile=expmap_file,
        irfs=cfg.irfs,
        roi=cfg.roi,
    )

    return {
        'evfile_filt': evfile_filt,
        'evfile_gti': evfile_gti,
        'ltcube': ltcube_file,
        'expmap': expmap_file,
    }
