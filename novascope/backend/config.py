"""
Settings for the backend. Change values here, not in the other files.

Every setting can also be overridden with an environment variable of the same name,
for example:  set PORT=9000   (Windows)   or   PORT=9000 python backend/app.py
"""
import os

# --- Server --------------------------------------------------------------------
HOST = os.environ.get("HOST", "127.0.0.1")       # 0.0.0.0 to accept other computers (needed when hosting online)
PORT = int(os.environ.get("PORT", 8000))
OPEN_BROWSER = os.environ.get("OPEN_BROWSER", "1") == "1"   # open the page automatically when started locally

# Demo mode: no internet needed. Search results and images are synthetic (made by demo.py).
DEMO = os.environ.get("DEMO", "0") == "1"

# --- Folders -------------------------------------------------------------------
BACKEND_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(BACKEND_DIR)
FRONTEND_DIR = os.path.join(PROJECT_DIR, "frontend")      # files served to the browser
CACHE_DIR = os.path.join(PROJECT_DIR, "cache")            # downloaded cutouts are kept here

# --- NASA / ESA / CDS services ---------------------------------------------------
IRSA_SIA = "https://irsa.ipac.caltech.edu/SIA"            # SPHEREx image search (IVOA SIA2)
CUTOUT_PREFIX = "https://irsa.ipac.caltech.edu/ibe/data/spherex/"   # only image links that start like this are allowed
GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"  # Gaia DR3 star catalogue
SESAME = "https://cds.unistra.fr/cgi-bin/nph-sesame/-ox/SNV?"       # object name -> coordinates

# SPHEREx collections the frontend may ask for. The first one is the default.
COLLECTIONS = ["spherex_dr1", "spherex_dr1_deep", "spherex_qr3", "spherex_qr3_deep"]

# --- Limits ----------------------------------------------------------------------
SEARCH_RADIUS_DEG = 0.0003       # images must contain this circle around the position (about 1 arcsec)
CUTOUT_MIN_DEG = 0.01            # smallest cutout side allowed
CUTOUT_MAX_DEG = 0.5             # largest cutout side allowed (bigger = slower downloads)
GAIA_MAX_RADIUS_DEG = 0.6
GAIA_MAX_ROWS = 4000
GAIA_MAG_LIMIT = 20              # faintest Gaia G magnitude returned

TIMEOUT_SEARCH = 120             # seconds
TIMEOUT_CUTOUT = 180
TIMEOUT_GAIA = 90
TIMEOUT_RESOLVE = 30

USER_AGENT = "starlight-time-machine/1.0 (NASA Space Apps project)"
