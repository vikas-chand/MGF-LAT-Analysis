"""Data handling: target catalogs, file I/O, and GTI utilities."""

from .targets import utc_to_met, build_target, load_catalog
from .io import find_ft1_ft2, ensure_nonempty_gti, save_profile, load_profile
