"""
Demo mode: synthetic SPHEREx-like data, so the app can be developed and shown without internet.
Turn it on with  python backend/app.py --demo  (or DEMO=1).

Everything here is MADE UP. The frontends label demo data as synthetic.

What it imitates:
  - three sky sweeps (May 2025, Nov 2025, May 2026), six bands each, two images per band
  - each sweep is rotated differently, like real SPHEREx images, so alignment code gets exercised
  - a fixed star field, plus one fast-moving object (10.4"/yr north, like Barnard's Star) at the centre
  - an asteroid that crosses the field during the second sweep, in band 3 only
  - the FITS layout of real cutouts: PRIMARY, IMAGE, VARIANCE, ZODI
"""
import math
import random
import struct
import sys
from array import array

PIXEL_DEG = 6.15 / 3600                     # SPHEREx pixel size
SWEEPS = [(60820.0, 0.0), (61005.0, 30.0), (61190.0, 60.0)]   # (MJD, field rotation in degrees)
BANDS = [(0.75, 1.09), (1.10, 1.62), (1.63, 2.41), (2.42, 3.82), (3.83, 4.41), (4.42, 5.00)]   # microns
PREFIX = "https://irsa.ipac.caltech.edu/ibe/data/spherex/demo/"
MOVER_PM = (-800.0, 10360.0)                # mas/yr in RA*cos(Dec) and Dec
TARGETS = {                                 # for name lookups in demo mode (J2000, mas/yr)
    "barnard's star": (269.452077, 4.693365, -801.551, 10362.394),
    "epsilon indi b": (331.043833, -56.782708, 3966.7, -2536.2),
    "wise 0855-0714": (133.7951321, -7.24515, -8118.9, 679.3),
    "rw aurigae": (76.956068, 30.401326, 4.9, -24.6),
    "r draconis": (248.16768, 66.754984, -20.07, 7.06),
}


def _year(mjd):
    return 2000 + (mjd - 51544.5) / 365.25


# ----------------------------------------------------------------------------- what spherex.py calls

def search_images(ra, dec, collection):
    rows = []
    for s, (mjd0, _) in enumerate(SWEEPS):
        for b, (lo, hi) in enumerate(BANDS, start=1):
            for k in range(2):
                t = mjd0 + b * 1.4 + k * 0.3
                rows.append({"obs_id": f"demo_{s + 1}_D{b}_{k + 1}", "t_min": t, "t_max": t + 0.0013,
                             "em_min": lo * 1e-6, "em_max": hi * 1e-6, "band": f"SPHEREx-D{b}",
                             "access_url": f"{PREFIX}{s}/{b}/{k}.fits"})
    return {"collection": collection, "count": len(rows), "rows": rows, "demo": True}


def gaia_stars(ra, dec, radius):
    stars = [{"ra": sra, "dec": sdec, "pmra": 0.0, "pmdec": 0.0, "g": 19 - 2.5 * math.log10(f)}
             for sra, sdec, f in _star_field(ra, dec)]
    mra, mdec = _mover(ra, dec, 2016.0)
    stars.append({"ra": mra, "dec": mdec, "pmra": MOVER_PM[0], "pmdec": MOVER_PM[1], "g": 8.2})
    return {"epoch": 2016.0, "rows": stars, "demo": True}


def resolve_name(name):
    t = TARGETS.get(name.lower().replace("−", "-"))
    if not t:
        raise LookupError(f'Demo mode only knows the five test targets, not "{name}".')
    return {"name": name, "ra": t[0], "dec": t[1], "pmra": t[2], "pmdec": t[3], "epoch": 2000.0, "demo": True}


def make_cutout(url, ra, dec, size):
    """A FITS cutout centred on (ra, dec), size degrees across, for one of the demo images."""
    sweep, band, k = (int(p) for p in url[len(PREFIX):].replace(".fits", "").split("/"))
    mjd0, rot = SWEEPS[sweep]
    mjd = mjd0 + band * 1.4 + k * 0.3 + 0.00065
    n = max(8, round(size / PIXEL_DEG))
    th = math.radians(rot)
    cd = [-PIXEL_DEG * math.cos(th), PIXEL_DEG * math.sin(th), PIXEL_DEG * math.sin(th), PIXEL_DEG * math.cos(th)]
    crpix = (n + 1) / 2
    rng = random.Random(hash((round(ra, 3), round(dec, 3), sweep, band, k)) & 0xFFFFFFFF)
    img = [rng.gauss(0, 0.05) for _ in range(n * n)]

    def add(sra, sdec, flux):
        x, y = _sky_to_pixel(sra, sdec, ra, dec, cd, crpix)
        x, y = x - 1, y - 1                                   # 1-based FITS -> 0-based index
        for j in range(max(0, int(y) - 4), min(n, int(y) + 5)):
            for i in range(max(0, int(x) - 4), min(n, int(x) + 5)):
                img[j * n + i] += flux * math.exp(-((i - x) ** 2 + (j - y) ** 2) / (2 * 0.75 ** 2))

    for sra, sdec, f in _star_field(ra, dec):
        add(sra, sdec, f * (0.7 + 0.1 * band))
    mra, mdec = _mover(ra, dec, _year(mjd))
    add(mra, mdec, 40 * (1 + 0.15 * band))
    if sweep == 1 and band == 3:                              # asteroid crossing during sweep 2
        add(ra + 0.02 + 0.004 * k, dec + 0.015 - 0.002 * k, 25)
    zodi = 0.3 + 0.1 * sweep
    img = [v + zodi for v in img]

    wcs = [("CTYPE1", "RA---TAN"), ("CTYPE2", "DEC--TAN"), ("CRVAL1", ra), ("CRVAL2", dec),
           ("CRPIX1", crpix), ("CRPIX2", crpix), ("CD1_1", cd[0]), ("CD1_2", cd[1]),
           ("CD2_1", cd[2]), ("CD2_2", cd[3]), ("CUNIT1", "deg"), ("CUNIT2", "deg")]
    meta = [("MJD-AVG", mjd), ("DETECTOR", band), ("BUNIT", "MJy/sr"), ("ORIGIN", "DEMO - SYNTHETIC DATA")]
    return (_hdu([("SIMPLE", True), ("BITPIX", 8), ("NAXIS", 0), ("EXTEND", True), ("COMMENT", "Synthetic demo data")])
            + _image_hdu("IMAGE", img, n, wcs + meta)
            + _image_hdu("VARIANCE", [0.0025] * (n * n), n, wcs)
            + _image_hdu("ZODI", [zodi] * (n * n), n, wcs))


# ----------------------------------------------------------------------------- sky helpers

def _star_field(ra, dec):
    """The same fixed stars every time this patch is asked for."""
    rng = random.Random(round(ra * 100) * 100003 + round(dec * 100))
    c = math.cos(math.radians(dec))
    return [(ra + rng.uniform(-0.17, 0.17) / c, dec + rng.uniform(-0.17, 0.17), rng.uniform(0.6, 1.0) ** -6)
            for _ in range(160)]


def _mover(ra, dec, year):
    """Position of the demo's fast-moving object: at the centre in late 2025, drifting like Barnard's Star."""
    dy = year - 2025.9
    return (ra + MOVER_PM[0] * dy / 3.6e6 / math.cos(math.radians(dec)), dec + MOVER_PM[1] * dy / 3.6e6)


def _sky_to_pixel(ra, dec, ra0, dec0, cd, crpix):
    """Gnomonic (TAN) projection, then the CD matrix inverse: sky degrees -> 1-based FITS pixel."""
    a, d, a0, d0 = map(math.radians, (ra, dec, ra0, dec0))
    cosc = math.sin(d0) * math.sin(d) + math.cos(d0) * math.cos(d) * math.cos(a - a0)
    xi = math.degrees(math.cos(d) * math.sin(a - a0) / cosc)
    eta = math.degrees((math.cos(d0) * math.sin(d) - math.sin(d0) * math.cos(d) * math.cos(a - a0)) / cosc)
    det = cd[0] * cd[3] - cd[1] * cd[2]
    return ((cd[3] * xi - cd[1] * eta) / det + crpix, (-cd[2] * xi + cd[0] * eta) / det + crpix)


# ----------------------------------------------------------------------------- tiny FITS writer

def _card(key, value):
    if key == "COMMENT":
        text = f"COMMENT {value}"
    elif isinstance(value, bool):
        text = f"{key:<8}= {'T' if value else 'F':>20}"
    elif isinstance(value, int):
        text = f"{key:<8}= {value:>20}"
    elif isinstance(value, float):
        text = f"{key:<8}= {value:>20.12E}"
    else:
        text = f"{key:<8}= '{value:<8}'"
    return text[:80].ljust(80)


def _hdu(cards, data=b""):
    header = "".join(_card(k, v) for k, v in cards) + "END".ljust(80)
    header += " " * (-len(header) % 2880)
    return header.encode("ascii") + data + b"\0" * (-len(data) % 2880)


def _image_hdu(name, values, n, extra):
    data = array("f", values)
    if sys.byteorder == "little":
        data.byteswap()                          # FITS stores numbers big-endian
    cards = [("XTENSION", "IMAGE"), ("BITPIX", -32), ("NAXIS", 2), ("NAXIS1", n), ("NAXIS2", n),
             ("PCOUNT", 0), ("GCOUNT", 1), ("EXTNAME", name)] + extra
    return _hdu(cards, data.tobytes())
