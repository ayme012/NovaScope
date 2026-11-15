"""
Step 3c: pull the positions of interstellar comet 3I/ATLAS out of a SPHEREx observation list.

Input : sources/obscore_3I_Atlas.vot
        IRSA ObsCore search result: 259 SPHEREx images (Quick Release 2) that contain 3I/ATLAS,
        each with the comet's position at that moment (comet_ra, comet_dec).
Output: atlas.json   {mjd: [...], ra: [...], dec: [...], det: [...]}  sorted by time
        (det = SPHEREx detector / band number 1-6)
Needs : astropy (only for reading the VOTable)
"""
import json
import warnings
import numpy as np
from astropy.io.votable import parse_single_table

warnings.simplefilter("ignore")
t = parse_single_table("sources/obscore_3I_Atlas.vot").to_table()
order = np.argsort(t["t_min"])
atlas = {
    "mjd": [round(float(t["t_min"][i]), 5) for i in order],
    "ra": [round(float(t["comet_ra"][i]), 5) for i in order],
    "dec": [round(float(t["comet_dec"][i]), 5) for i in order],
    "det": [int(str(t["spherex_detector"][i])[1]) for i in order],   # "D3" -> 3
}
json.dump(atlas, open("atlas.json", "w"), separators=(",", ":"))
print(len(atlas["mjd"]), "positions,", atlas["mjd"][0], "to", atlas["mjd"][-1], "(MJD)")
