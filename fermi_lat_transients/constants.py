"""
Physical constants and unit conversions for Fermi-LAT transient analysis.

All values in CGS units unless otherwise noted.
"""

from astropy.time import Time

# ─── Fundamental constants (CGS) ────────────────────────────────────────────
c_cgs = 2.99792458e10            # speed of light [cm/s]
a_rad = 7.5657e-15               # radiation constant [erg cm^-3 K^-4]
kB_erg = 1.380649e-16            # Boltzmann constant [erg/K]
sigma_T = 6.6524587321e-25       # Thomson cross-section [cm^2]
m_p = 1.67262192369e-24          # proton mass [g]
m_e = 9.1093837015e-28           # electron mass [g]

# ─── Distance / cosmology ───────────────────────────────────────────────────
Mpc_cm = 3.085677581491367e24    # 1 Mpc in cm
kpc_cm = 3.085677581491367e21    # 1 kpc in cm
pc_cm = 3.085677581491367e18     # 1 pc in cm

# ─── Energy unit conversions ────────────────────────────────────────────────
MeV2erg = 1.602176634e-6         # 1 MeV in erg
GeV2erg = 1.602176634e-3         # 1 GeV in erg
keV2erg = 1.602176634e-9         # 1 keV in erg
erg2MeV = 1.0 / MeV2erg         # 1 erg in MeV

# ─── Fermi mission reference ────────────────────────────────────────────────
FERMI_EPOCH = Time('2001-01-01T00:00:00', scale='utc')
"""Fermi Mission Elapsed Time (MET) epoch: 2001-01-01 00:00:00 UTC."""

# ─── Statistical thresholds ─────────────────────────────────────────────────
DELTA_TS_95 = 2.71               # Delta-TS for 95% CL one-sided upper limit
DELTA_TS_99 = 6.63               # Delta-TS for 99% CL one-sided upper limit
TS_DETECTION = 25.0              # TS threshold for source detection (~5 sigma)
