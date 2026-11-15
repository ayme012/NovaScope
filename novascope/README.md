# NovaScope: Navigate the stars, Discover the Universe

NASA Space Apps 2026, "Planet X and SPHEREx".
A star map you can run 200,000 years back or forward, plus a SPHEREx time-lapse:
enter coordinates and watch that patch of sky change across SPHEREx's scans.

## Run it

You need Python 3.8 or newer. **Nothing to `pip install`**: the backend uses only Python's
standard library, so Windows Smart App Control won't block it.

From this folder, in a terminal (VS Code: Ctrl + `):

```
py backend/app.py            # Windows (or: python backend/app.py)
                             # Windows shortcut: double-click start-novascope.bat
python3 backend/app.py       # Mac / Linux
```

Your browser opens **http://localhost:8000**. Stop with **Ctrl + C**.

| Mode | Command | What you get |
|---|---|---|
| Real data | `py backend/app.py` | SPHEREx images from NASA's IRSA archive, Gaia DR3 stars, name lookups |
| Demo | `py backend/app.py --demo` | Synthetic images, no internet needed. Pages label them "Demo data (synthetic)" |
| No backend | double-click `frontend/index.html` | Star map works; SPHEREx tab plays simulated images and real M101 frames |

Other options: `--port 9000`, `--host 0.0.0.0` (let other computers connect), `--no-browser`.

## Folders

```
novascope/
├── backend/                 Python web server (standard library only)
│   ├── app.py               starts the server, defines the /api routes
│   ├── config.py            ALL settings: port, archive URLs, limits, timeouts
│   ├── spherex.py           talks to the archives: IRSA, Gaia, CDS
│   └── demo.py              synthetic data for --demo
├── frontend/                what the browser loads
│   ├── index.html           star map page layout
│   ├── style.css            star map styling
│   ├── app.js               star map logic (sections marked // ---------- name ----------)
│   ├── data/stars.js        star catalogue, constellation lines, Milky Way (generated)
│   ├── data/spherex.js      real M101 frames, IRSA figures, 3I/ATLAS track (generated)
│   └── timelapse/           full SPHEREx time-lapse tool (index.html, timelapse.css, timelapse.js)
├── tools/                   scripts that generated frontend/data (optional, need pip packages)
├── cache/                   downloaded cutouts (created when needed, safe to delete)
├── Dockerfile               for hosting online
└── requirements.txt         empty on purpose: no packages needed
```

## Where to change things

| You want to… | Edit | What to change |
|---|---|---|
| Use another port | `backend/config.py` | `PORT` |
| Search a different data release by default | `backend/config.py` | order of `COLLECTIONS` (first = default) |
| Allow bigger or smaller cutouts | `backend/config.py` | `CUTOUT_MIN_DEG`, `CUTOUT_MAX_DEG` |
| Get more or fewer Gaia stars | `backend/config.py` | `GAIA_MAX_ROWS`, `GAIA_MAG_LIMIT` |
| Wait longer for slow archives | `backend/config.py` | `TIMEOUT_*` |
| Add a new API endpoint | `backend/app.py` | write `def api_x(q): ... return {...}`, add it to `ROUTES` |
| Change how IRSA is searched | `backend/spherex.py` | `search_images()` |
| Download more images in the star map's SPHEREx tab | `frontend/app.js` | `const LIVE = { collection, maxImages, parallel }` |
| Change the test-target buttons | `frontend/app.js` and `frontend/timelapse/timelapse.js` | `TL_PRESETS` / `PRESETS` |
| Change colours or fonts | `frontend/style.css` | the `:root` block at the top |
| Change the guided tours | `frontend/app.js` | `JOURNEYS` |
| Add or edit a star's fun facts | `frontend/app.js` | `STAR_FACTS` (key = star name, or `"HIP 12345"`) |
| Change the "light left this star when…" events | `frontend/app.js` | `LIGHT_EVENTS` |
| Speed of the 3I/ATLAS animation | `frontend/app.js` | `ATLAS_DAY` (seconds per real day), `ATLAS_GAP` (seconds to cross the Sep–Nov gap) |
| Timing of the M101 zoom animation | `frontend/app.js` | `M101_T` (seconds for each step), `m101Spot()` (where it sits and how big) |
| The enlarged player on the map (click the video square) | `frontend/app.js` / `style.css` | `cutBigHalf()` (how big), `cutSpot()` (where), `openCutPop` / `.cut-pop` (control bar) |
| Zoom speed and limits | `frontend/app.js` | wheel listener (`k`, how fast each scroll step zooms), `FOV_MIN` / `FOV_MAX` |

After editing a Python file, stop the server (Ctrl + C) and start it again.
After editing a frontend file, just reload the page.

## API

All requests are `GET`. Answers are JSON, except `/api/cutout`, which answers a FITS file.
Errors answer `{"error": "..."}` with HTTP 400 (bad request) or 502 (archive problem).

| Endpoint | Parameters | Answer | Archive behind it |
|---|---|---|---|
| `/api/ping` | none | `{"ok": true, "demo": false}` | none |
| `/api/search` | `ra`, `dec` (degrees), `collection` (optional) | `{"collection", "count", "rows": [{"obs_id", "t_min", "t_max", "em_min", "em_max", "band", "access_url"}]}`, oldest first; `t_*` in MJD, `em_*` in metres | IRSA SIA2 `https://irsa.ipac.caltech.edu/SIA` |
| `/api/cutout` | `url` (an `access_url` from search), `ra`, `dec`, `size` (degrees, default 0.1) | FITS file: IMAGE, FLAGS, VARIANCE, ZODI, ... | IRSA cutout: `<url>?center=ra,dec&size=deg` |
| `/api/gaia` | `ra`, `dec`, `r` (radius, degrees) | `{"epoch": 2016.0, "rows": [{"ra", "dec", "pmra", "pmdec", "g"}]}` | ESA Gaia archive, `gaiadr3.gaia_source` |
| `/api/resolve` | `name` | `{"name", "ra", "dec", "pmra", "pmdec", "epoch": 2000.0}` | CDS Sesame |

Example: `http://localhost:8000/api/search?ra=269.44629&dec=4.76792`

Cutout links are checked: only `https://irsa.ipac.caltech.edu/ibe/data/spherex/...` is accepted,
so the server can't be used to download from other sites.

## Put it online

Docker:

```
docker build -t starlight .
docker run -p 8000:8000 starlight
```

Render.com (free): push this folder to GitHub, create a Web Service from it, set the start command
to `python backend/app.py --host 0.0.0.0 --no-browser`. Render sets `PORT` automatically.

## Data credits

- Stars: HYG Database v4.1 by Astronexus (CC BY-SA 4.0).
- Constellation lines and Milky Way outline: d3-celestial by Olaf Frohn (BSD).
- SPHEREx images and tutorial figures: NASA/IPAC Infrared Science Archive (IRSA); IRSA tutorials © Caltech (BSD 3-Clause).
- Gaia DR3: ESA/Gaia/DPAC. Name lookups: CDS Sesame, Strasbourg.
- 3I/ATLAS track: SPHEREx observation list (ObsCore, Quick Release 2) from IRSA.

## Known limits

- The star map includes stars brighter than magnitude 7 today, moving in straight lines.
- The demo backend's images are synthetic. Real-data mode was tested against a simulated copy of IRSA's
  answers; the first run against the real archive is the true test.
- Some 3I/ATLAS image links point to Quick Release 2, which IRSA has replaced with DR1 and QR3.
