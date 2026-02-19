"""Likelihood analysis: unbinned fitting, profile scans, and SED construction."""

from .unbinned import (
    setup_unbinned,
    fit,
    freeze_all_params,
    free_diffuse_norms,
    prepare_transient_pl2,
    set_flux,
)
from .profile import (
    profile_scan,
    upper_limit_from_profile,
)
from .sed import (
    sed_profile_scan,
    compute_sed_points,
)
