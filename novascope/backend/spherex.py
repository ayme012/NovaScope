"""
Talking to the data archives. Each public function here is one thing the frontend can ask for:

    search_images(ra, dec, collection)  -> list of SPHEREx images covering a position   (NASA IRSA)
    get_cutout(url, ra, dec, size)      -> FITS bytes of a small square cut from one image (NASA IRSA)
    gaia_stars(ra, dec, radius)         -> Gaia DR3 stars around a position                (ESA Gaia archive)
    resolve_name(name)                  -> coordinates and motion of a named object        (CDS Sesame)

Only Python's standard library is used, so there is nothing to pip install.
In demo mode (config.DEMO) the same functions return synthetic data from demo.py.
"""
import csv
import hashlib
import io
import os
import threading
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

import config
import demo


# ----------------------------------------------------------------------------- helpers

def fetch(url, timeout=120):
    """Download a URL and return its bytes."""
    req = urllib.request.Request(url, headers={"User-Agent": config.USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def _tag(element):
    """XML tag without its namespace: '{http://...}FIELD' -> 'FIELD'."""
    return element.tag.rsplit("}", 1)[-1]


def _num(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_votable(data):
    """Read a VOTable (the XML table format archives answer with). Returns a list of dicts."""
    root = ET.fromstring(data)
    for el in root.iter():
        if _tag(el) == "INFO" and el.get("name") == "QUERY_STATUS" and el.get("value") == "ERROR":
            raise RuntimeError((el.text or "The archive reported an error").strip())
    names = [f.get("name") for f in root.iter() if _tag(f) == "FIELD"]
    rows = []
    for tr in root.iter():
        if _tag(tr) == "TR":
            cells = [(td.text or "").strip() for td in tr if _tag(td) == "TD"]
            rows.append(dict(zip(names, cells)))
    return rows


# ----------------------------------------------------------------------------- SPHEREx

def search_images(ra, dec, collection=None):
    """
    Every SPHEREx image that covers (ra, dec), oldest first.
    Each image: obs_id, t_min / t_max (MJD), em_min / em_max (metres), band name, access_url.
    """
    collection = collection or config.COLLECTIONS[0]
    if collection not in config.COLLECTIONS:
        raise ValueError(f"Unknown collection {collection}. Allowed: {', '.join(config.COLLECTIONS)}")
    if config.DEMO:
        return demo.search_images(ra, dec, collection)

    params = {"COLLECTION": collection, "POS": f"CIRCLE {ra} {dec} {config.SEARCH_RADIUS_DEG}"}
    url = config.IRSA_SIA + "?" + urllib.parse.urlencode(params, quote_via=urllib.parse.quote)
    images = []
    for r in parse_votable(fetch(url, config.TIMEOUT_SEARCH)):
        if not r.get("access_url"):
            continue
        images.append({
            "obs_id": r.get("obs_id", ""),
            "t_min": _num(r.get("t_min")), "t_max": _num(r.get("t_max")),
            "em_min": _num(r.get("em_min")), "em_max": _num(r.get("em_max")),
            "band": r.get("energy_bandpassname", ""),
            "access_url": r["access_url"],
        })
    images.sort(key=lambda i: i["t_min"] or 0)
    return {"collection": collection, "count": len(images), "rows": images}


def get_cutout(url, ra, dec, size):
    """
    FITS bytes of a size x size degree square around (ra, dec), cut from one SPHEREx image.
    Cutouts are cached in config.CACHE_DIR, so asking twice only downloads once.
    """
    if not url.startswith(config.CUTOUT_PREFIX) or "?" in url:
        raise ValueError("Only SPHEREx image links from IRSA are allowed")
    size = max(config.CUTOUT_MIN_DEG, min(config.CUTOUT_MAX_DEG, size))
    if config.DEMO:
        return demo.make_cutout(url, ra, dec, size)

    full = f"{url}?center={ra:.7f},{dec:.7f}&size={size:.4f}"
    os.makedirs(config.CACHE_DIR, exist_ok=True)
    path = os.path.join(config.CACHE_DIR, hashlib.sha1(full.encode()).hexdigest() + ".fits")
    if os.path.exists(path):
        with open(path, "rb") as f:
            return f.read()
    data = fetch(full, config.TIMEOUT_CUTOUT)
    if not data.startswith(b"SIMPLE"):
        raise RuntimeError("IRSA did not return a FITS file for this cutout")
    tmp = f"{path}.{threading.get_ident()}.part"          # write then rename, so a half-written file is never used
    with open(tmp, "wb") as f:
        f.write(data)
    os.replace(tmp, path)
    return data


# ----------------------------------------------------------------------------- Gaia

def gaia_stars(ra, dec, radius):
    """Gaia DR3 stars within radius degrees: ra, dec (epoch 2016.0), pmra, pmdec (mas/yr), g (magnitude)."""
    radius = min(config.GAIA_MAX_RADIUS_DEG, radius)
    if config.DEMO:
        return demo.gaia_stars(ra, dec, radius)
    adql = (f"SELECT TOP {config.GAIA_MAX_ROWS} source_id, ra, dec, pmra, pmdec, phot_g_mean_mag "
            "FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS', ra, dec), "
            f"CIRCLE('ICRS', {ra}, {dec}, {radius})) AND phot_g_mean_mag < {config.GAIA_MAG_LIMIT} "
            "ORDER BY phot_g_mean_mag")
    url = config.GAIA_TAP + "?" + urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": adql})
    text = fetch(url, config.TIMEOUT_GAIA).decode("utf-8", "replace")
    stars = [{"ra": _num(r["ra"]), "dec": _num(r["dec"]), "pmra": _num(r["pmra"]) or 0.0,
              "pmdec": _num(r["pmdec"]) or 0.0, "g": _num(r["phot_g_mean_mag"])}
             for r in csv.DictReader(io.StringIO(text))]
    return {"epoch": 2016.0, "rows": stars}


# ----------------------------------------------------------------------------- names

def resolve_name(name):
    """Coordinates (J2000, epoch 2000) and proper motion of a named object, e.g. "Barnard's Star"."""
    name = name.strip()
    if config.DEMO:
        return demo.resolve_name(name)
    root = ET.fromstring(fetch(config.SESAME + urllib.parse.quote(name), config.TIMEOUT_RESOLVE))
    for resolver in root.iter():
        if _tag(resolver) != "Resolver":
            continue
        v = {_tag(c): (c.text or "").strip() for c in resolver}
        if "jradeg" in v and "jdedeg" in v:
            return {"name": name, "ra": float(v["jradeg"]), "dec": float(v["jdedeg"]),
                    "pmra": _num(v.get("pmRA")) or 0.0, "pmdec": _num(v.get("pmDE")) or 0.0, "epoch": 2000.0}
    raise LookupError(f'No object called "{name}" was found. Try coordinates instead.')
