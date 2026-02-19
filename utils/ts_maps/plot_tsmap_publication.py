#!/usr/bin/env python3
"""
Publication-quality TS map plotter for Fermi-LAT MGF analysis.
Works without astropy — reads FITS headers/data manually and uses
a linearized WCS appropriate for ~16° fields.

Author: Generated for Vikas's MGF LAT analysis
"""
import math, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patheffects
import xml.etree.ElementTree as ET

# ── Journal-quality matplotlib defaults ──────────────────────────────
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['DejaVu Serif', 'Times New Roman', 'Times',
                   'Computer Modern Roman'],
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'text.usetex': False,
    'mathtext.fontset': 'dejavuserif',
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.minor.width': 0.5,
    'ytick.minor.width': 0.5,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
})


# ══════════════════════════════════════════════════════════════════════
#  FITS reader (no astropy)
# ══════════════════════════════════════════════════════════════════════
def read_fits_image(filepath):
    """Read a 2D image from a simple FITS primary HDU."""
    with open(filepath, 'rb') as f:
        hdr_bytes = b''
        while True:
            block = f.read(2880)
            if not block:
                raise IOError("EOF reading FITS header")
            hdr_bytes += block
            for i in range(0, len(block), 80):
                if block[i:i+3] == b'END':
                    break
            else:
                continue
            break

        hdr = {}
        for i in range(0, len(hdr_bytes), 80):
            card = hdr_bytes[i:i+80].decode('ascii', errors='replace')
            if card.startswith('END'):
                break
            if '=' in card[:9]:
                key = card[:8].strip()
                val = card[10:].split('/')[0].strip()
                if val.startswith("'"):
                    hdr[key] = val.strip("'").strip()
                elif val in ('T', 'F'):
                    hdr[key] = val == 'T'
                else:
                    try:
                        hdr[key] = (float(val)
                                    if ('.' in val or 'E' in val.upper())
                                    else int(val))
                    except ValueError:
                        hdr[key] = val

        bitpix = hdr['BITPIX']
        nx = int(hdr['NAXIS1'])
        ny = int(hdr['NAXIS2'])
        dtype = {-32: '>f4', -64: '>f8', 16: '>i2', 32: '>i4'}[bitpix]
        data = np.frombuffer(f.read(nx * ny * abs(bitpix) // 8),
                             dtype=dtype).reshape(ny, nx).astype(np.float64)
    return data, hdr


# ══════════════════════════════════════════════════════════════════════
#  Linearized WCS for AIT small-field TS maps
# ══════════════════════════════════════════════════════════════════════
class LinearWCS:
    """Linearized WCS valid for the ~16° Fermi-LAT TS map fields.
    Uses the standard plate-carrée local approximation:
        x_iwc = -(RA - CRVAL1)*cos(CRVAL2)   y_iwc = Dec - CRVAL2
    with CDELT mapping IWC → pixel.
    """
    def __init__(self, hdr):
        self.crpix1 = float(hdr['CRPIX1'])      # 1-based ref pixel
        self.crpix2 = float(hdr['CRPIX2'])
        self.crval1 = float(hdr['CRVAL1'])       # RA  of ref pixel (deg)
        self.crval2 = float(hdr['CRVAL2'])       # Dec of ref pixel (deg)
        self.cdelt1 = float(hdr['CDELT1'])       # deg/pix (negative for RA)
        self.cdelt2 = float(hdr['CDELT2'])       # deg/pix

    def world_to_pixel(self, ra_deg, dec_deg):
        """World (RA, Dec) → 0-based pixel (x, y)."""
        dra = ra_deg - self.crval1
        # Wrap RA difference to [-180, 180]
        if dra > 180:
            dra -= 360
        elif dra < -180:
            dra += 360
        ddec = dec_deg - self.crval2
        # IWC  (the cos(CRVAL2) factor converts RA offset to on-sky offset)
        x_iwc = dra * math.cos(math.radians(self.crval2))
        y_iwc = ddec
        # Pixel (0-based).  CDELT1 is negative → RA increases to the left
        xpix = (self.crpix1 - 1.0) + x_iwc / self.cdelt1
        ypix = (self.crpix2 - 1.0) + y_iwc / self.cdelt2
        return xpix, ypix

    def pixel_to_world(self, xpix, ypix):
        """0-based pixel → World (RA, Dec)."""
        x_iwc = (xpix - (self.crpix1 - 1.0)) * self.cdelt1
        y_iwc = (ypix - (self.crpix2 - 1.0)) * self.cdelt2
        dec_deg = self.crval2 + y_iwc
        cos_ref = math.cos(math.radians(self.crval2))
        ra_deg  = self.crval1 + x_iwc / cos_ref if cos_ref > 1e-6 else self.crval1
        return ra_deg, dec_deg


# ══════════════════════════════════════════════════════════════════════
#  Helpers
# ══════════════════════════════════════════════════════════════════════
def ang_sep(ra1, dec1, ra2, dec2):
    """Great-circle separation (degrees)."""
    d = math.pi / 180
    val = (math.sin(dec1*d)*math.sin(dec2*d) +
           math.cos(dec1*d)*math.cos(dec2*d)*math.cos((ra1-ra2)*d))
    return math.degrees(math.acos(max(-1, min(1, val))))


def parse_xml_sources(xml_path):
    """Parse Fermi-LAT XML model → [(name, ra, dec), ...]."""
    tree = ET.parse(xml_path)
    out = []
    for src in tree.getroot().findall('.//source'):
        if src.get('type') != 'PointSource':
            continue
        sp = src.find('.//spatialModel')
        if sp is None:
            continue
        ra_p  = sp.find("./parameter[@name='RA']")
        dec_p = sp.find("./parameter[@name='DEC']")
        if ra_p is not None and dec_p is not None:
            out.append((src.get('name',''),
                        float(ra_p.get('value')),
                        float(dec_p.get('value'))))
    return out


def _smooth(img, sigma=1.2):
    """Simple Gaussian smooth using numpy only (no scipy needed)."""
    ksize = int(4 * sigma + 0.5) * 2 + 1
    x = np.arange(ksize) - ksize // 2
    kern1d = np.exp(-0.5 * (x / sigma) ** 2)
    kern1d /= kern1d.sum()
    # separable 2D convolution
    out = np.apply_along_axis(lambda r: np.convolve(r, kern1d, mode='same'), 1, img)
    out = np.apply_along_axis(lambda c: np.convolve(c, kern1d, mode='same'), 0, out)
    return out


def nice_ticks(vmin, vmax, target_n=5):
    """Return nicely-spaced tick values."""
    rng = vmax - vmin
    if rng == 0:
        return [vmin]
    raw = rng / target_n
    mag = 10 ** math.floor(math.log10(max(abs(raw), 1e-10)))
    for s in [1, 2, 5, 10]:
        step = mag * s
        if step >= raw * 0.7:
            break
    start = math.ceil(vmin / step) * step
    ticks = []
    v = start
    while v <= vmax + step * 0.01:
        ticks.append(round(v, 10))
        v += step
    return ticks


# ══════════════════════════════════════════════════════════════════════
#  Main plotter
# ══════════════════════════════════════════════════════════════════════
def plot_tsmap_publication(
    nosrc_fits, withgrb_fits, xml_path,
    src_name, src_ra, src_dec,
    roi_radius=12.0, grb_trigger=None,
    output_dir='.', vmin=0, vmax=None,
    cmap='inferno', label_radius=7.0,
    show_labels=True, max_labels=10,
):
    # ── Data ──────────────────────────────────────────────────────
    data_nosrc,   hdr_nosrc   = read_fits_image(nosrc_fits)
    data_withgrb, hdr_withgrb = read_fits_image(withgrb_fits)

    img_nosrc   = np.sqrt(np.maximum(data_nosrc,   0))
    img_withgrb = np.sqrt(np.maximum(data_withgrb, 0))

    wcs_nosrc   = LinearWCS(hdr_nosrc)
    wcs_withgrb = LinearWCS(hdr_withgrb)

    if vmax is None:
        vmax = max(np.nanmax(img_nosrc), np.nanmax(img_withgrb))
        vmax = math.ceil(vmax * 10) / 10          # round up to 0.1

    # ── Catalog sources ───────────────────────────────────────────
    all_srcs = parse_xml_sources(xml_path)
    nearby = [(n, r, d) for n, r, d in all_srcs
              if ang_sep(r, d, src_ra, src_dec) < roi_radius + 2]
    labeled = sorted(
        [(n, r, d, ang_sep(r, d, src_ra, src_dec))
         for n, r, d in nearby
         if ang_sep(r, d, src_ra, src_dec) < label_radius],
        key=lambda x: x[3])[:max_labels]
    labeled_set = {t[0] for t in labeled}

    # ── Figure ────────────────────────────────────────────────────
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.8))
    fig.subplots_adjust(wspace=0.35)

    panels = [
        (ax1, img_nosrc,   wcs_nosrc,
         f"$\\sqrt{{\\mathrm{{TS}}}}$ (no {src_name})"),
        (ax2, img_withgrb, wcs_withgrb,
         f"$\\sqrt{{\\mathrm{{TS}}}}$ ({src_name} in model)"),
    ]

    for ax, img, wcs, title in panels:
        ny, nx = img.shape

        # ── Smoothed image ──
        im = ax.imshow(img, origin='lower', vmin=vmin, vmax=vmax,
                        cmap=cmap, aspect='equal',
                        interpolation='gaussian')
        ax.set_xlim(-0.5, nx - 0.5)
        ax.set_ylim(-0.5, ny - 0.5)
        ax.set_facecolor('black')

        # ── Colorbar ──
        cbar = fig.colorbar(im, ax=ax, pad=0.02,
                            fraction=0.046, shrink=0.92)
        cbar.set_label(r'$\sqrt{\mathrm{TS}}$', fontsize=11)
        cbar.ax.tick_params(labelsize=9, width=0.5)
        cbar.outline.set_linewidth(0.5)

        # ── Title ──
        ax.set_title(title, fontsize=11, pad=8)

        # ── Axis ticks (RA / Dec) ──
        _set_wcs_ticks(ax, wcs, nx, ny)

        # ── Catalog markers ──
        pe_mark = [patheffects.withStroke(linewidth=2.2, foreground='black')]
        pe_text = [patheffects.withStroke(linewidth=1.5, foreground='black',
                                          alpha=0.85)]

        for name, ra, dec in nearby:
            xp, yp = wcs.world_to_pixel(ra, dec)
            if -3 < xp < nx + 3 and -3 < yp < ny + 3:
                ax.plot(xp, yp, '+', color='#00dddd',
                        ms=5, mew=0.7, alpha=0.65,
                        zorder=15, path_effects=pe_mark)

                if show_labels and name in labeled_set:
                    short = (name.replace('4FGL ', '')
                                  .replace('4FGL_', ''))
                    ax.annotate(
                        short, (xp, yp),
                        xytext=(4, 4), textcoords='offset points',
                        fontsize=5, color='#88eeee', alpha=0.8,
                        path_effects=pe_text, zorder=16)

        # ── GRB position marker ──
        xg, yg = wcs.world_to_pixel(src_ra, src_dec)
        pe_grb = [patheffects.withStroke(linewidth=3.5, foreground='black')]
        ax.plot(xg, yg, 'x', color='#ff4040', ms=10, mew=2.2,
                zorder=30, path_effects=pe_grb)

        # ── ROI circle ──
        _draw_sky_circle(ax, wcs, src_ra, src_dec, roi_radius,
                         ls='--', lw=1.5, color='white', alpha=0.75)

    # ── Shared axis labels ──
    ax1.set_xlabel('R.A. (J2000)', fontsize=11)
    ax1.set_ylabel('Decl. (J2000)', fontsize=11)
    ax2.set_xlabel('R.A. (J2000)', fontsize=11)
    ax2.set_ylabel('Decl. (J2000)', fontsize=11)

    # ── Save ──
    os.makedirs(output_dir, exist_ok=True)
    png = os.path.join(output_dir, f'{src_name}_TSmap.png')
    pdf = os.path.join(output_dir, f'{src_name}_TSmap.pdf')
    fig.savefig(png, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(pdf, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: {png}")
    print(f"  Saved: {pdf}")
    return png


# ── axis tick helpers ─────────────────────────────────────────────────
def _set_wcs_ticks(ax, wcs, nx, ny):
    """Place RA / Dec tick labels derived from the WCS."""
    # Corner world coords
    corners = [wcs.pixel_to_world(x, y)
               for x in (0, nx-1) for y in (0, ny-1)]
    ras  = [c[0] for c in corners]
    decs = [c[1] for c in corners]

    ra_lo, ra_hi   = min(ras), max(ras)
    dec_lo, dec_hi = min(decs), max(decs)

    # RA increases right-to-left (pixel 0 is high RA), ensure order
    ra_cen = wcs.crval1

    # ── RA ticks along the middle row ──
    ra_ticks = nice_ticks(ra_lo, ra_hi, 5)
    xpos, xlbl = [], []
    dec_mid = wcs.crval2
    for rt in ra_ticks:
        xp, _ = wcs.world_to_pixel(rt, dec_mid)
        if -1 < xp < nx + 1:
            xpos.append(xp)
            xlbl.append(f'{rt:.0f}°')
    ax.set_xticks(xpos)
    ax.set_xticklabels(xlbl)

    # ── Dec ticks along the middle column ──
    dec_ticks = nice_ticks(dec_lo, dec_hi, 5)
    ypos, ylbl = [], []
    for dt in dec_ticks:
        _, yp = wcs.world_to_pixel(ra_cen, dt)
        if -1 < yp < ny + 1:
            ypos.append(yp)
            sign = '+' if dt >= 0 else ''
            ylbl.append(f'{sign}{dt:.0f}°')
    ax.set_yticks(ypos)
    ax.set_yticklabels(ylbl)


def _draw_sky_circle(ax, wcs, ra0, dec0, r_deg, n=120, **kw):
    """Draw a great-circle of angular radius *r_deg* centred on (ra0, dec0)."""
    d = math.radians
    xs, ys = [], []
    for phi in np.linspace(0, 2*math.pi, n):
        dec_p = math.asin(
            math.sin(d(dec0))*math.cos(d(r_deg)) +
            math.cos(d(dec0))*math.sin(d(r_deg))*math.cos(phi))
        dra = math.atan2(
            math.sin(phi)*math.sin(d(r_deg))*math.cos(d(dec0)),
            math.cos(d(r_deg)) - math.sin(d(dec0))*math.sin(dec_p))
        ra_p  = ra0 + math.degrees(dra)
        dec_pd = math.degrees(dec_p)
        xp, yp = wcs.world_to_pixel(ra_p, dec_pd)
        xs.append(xp); ys.append(yp)
    ax.plot(xs, ys, zorder=20, **kw)


# ══════════════════════════════════════════════════════════════════════
#  Event table
# ══════════════════════════════════════════════════════════════════════
BASE = '/sessions/happy-optimistic-cerf/mnt/MGFs_LAT/MGF_LAT_analysis_trigger'
OUT  = f'{BASE}/TS_Maps'

EVENTS = {
    'GRB200415A': dict(
        src_ra=11.888058, src_dec=-25.288800,
        trigger='bn200415367', roi=3.0,
        nosrc=f'{BASE}/GRB_200415A/GRB200415A_Analysis/tsmap_nogrb_psc35.fits',
        withgrb=f'{BASE}/GRB_200415A/GRB200415A_Analysis/tsmap_with_grb.fits',
        xml=f'{BASE}/GRB_200415A/GRB200415A_Analysis/model_without_grb_psc35.xml',
    ),
    'GRB081213A': dict(
        src_ra=11.888058, src_dec=-25.288800,
        trigger='bn081213173', roi=3.0,
        nosrc=f'{BASE}/GRB_081213A/GRB081213A_Analysis/tsmap_nogrb_psc35.fits',
        withgrb=f'{BASE}/GRB_081213A/GRB081213A_Analysis/tsmap_with_grb.fits',
        xml=f'{BASE}/GRB_081213A/GRB081213A_Analysis/model_without_grb_psc35.xml',
    ),
    'GRB180128A': dict(
        src_ra=11.888058, src_dec=-25.288800,
        trigger='bn180128215', roi=3.0,
        nosrc=f'{BASE}/GRB_180128A/GRB180128A_Analysis/tsmap_nogrb_psc35.fits',
        withgrb=f'{BASE}/GRB_180128A/GRB180128A_Analysis/tsmap_with_grb.fits',
        xml=f'{BASE}/GRB_180128A/GRB180128A_Analysis/model_without_grb_psc35.xml',
    ),
    'GRB200423A': dict(
        src_ra=308.718050, src_dec=60.153678,
        trigger='bn200423579', roi=3.0,
        nosrc=f'{BASE}/GRB_200423A/GRB200423A_Analysis/tsmap_nogrb_psc35.fits',
        withgrb=f'{BASE}/GRB_200423A/GRB200423A_Analysis/tsmap_with_grb.fits',
        xml=f'{BASE}/GRB_200423A/GRB200423A_Analysis/model_without_grb_psc35.xml',
    ),
    'GRB231024A': dict(
        src_ra=11.888058, src_dec=-25.288800,
        trigger='bn231024556', roi=3.0,
        nosrc=f'{BASE}/GRB_231024A/GRB231024A_Analysis/tsmap_nogrb_psc35.fits',
        withgrb=f'{BASE}/GRB_231024A/GRB231024A_Analysis/tsmap_with_grb.fits',
        xml=f'{BASE}/GRB_231024A/GRB231024A_Analysis/model_without_grb_psc35.xml',
    ),
    'GRB231115A': dict(
        src_ra=148.968458, src_dec=69.679703,
        trigger='bn231115650', roi=3.0,
        nosrc=f'{BASE}/GRB_231115A/GRB231115A_Analysis/tsmap_nogrb_psc35.fits',
        withgrb=f'{BASE}/GRB_231115A/GRB231115A_Analysis/tsmap_with_grb.fits',
        xml=f'{BASE}/GRB_231115A/GRB231115A_Analysis/model_without_grb_psc35.xml',
    ),
}


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    targets = sys.argv[1:] if len(sys.argv) > 1 else list(EVENTS.keys())
    for name in targets:
        if name not in EVENTS:
            print(f"Unknown event: {name}"); continue
        ev = EVENTS[name]
        print(f"\n{'='*50}\n  Plotting {name}\n{'='*50}")
        plot_tsmap_publication(
            nosrc_fits=ev['nosrc'], withgrb_fits=ev['withgrb'],
            xml_path=ev['xml'], src_name=name,
            src_ra=ev['src_ra'], src_dec=ev['src_dec'],
            roi_radius=ev['roi'], grb_trigger=ev['trigger'],
            output_dir=OUT, cmap='inferno',
            label_radius=7.0, max_labels=10,
        )
