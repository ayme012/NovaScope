"""
Step 1: turn the Milky Way outline into a cloud of sky points.

Input : sources/mw.json   (from the d3-celestial npm package, data/mw.json, BSD licence)
Output: milky_way_points.npy   columns: RA-ish longitude (deg, -180..180), Dec (deg), brightness level 1-5

The outline is five nested polygons (brightness levels). We scatter 40,000 evenly spaced
points over the sphere and count how many polygons contain each point.
"""
import json
import numpy as np
from matplotlib.path import Path

m = json.load(open("sources/mw.json"))
N = 40000
i = np.arange(N) + 0.5
lat = np.degrees(np.arcsin(1 - 2 * i / N))                      # Fibonacci sphere: even spacing
lon = np.degrees(np.pi * (1 + 5 ** 0.5) * i) % 360 - 180

level = np.zeros(N, int)
for feature in m["features"]:
    rings = [np.array(r, float) for poly in feature["geometry"]["coordinates"] for r in poly]
    inside = np.zeros(N, bool)
    for r in rings:
        L = np.degrees(np.unwrap(np.radians(r[:, 0])))
        B = r[:, 1]
        if abs(L[-1] - L[0]) > 180:          # ring wraps all the way round the sky: close it through the south pole
            L = np.r_[L, L[-1], L[0]]
            B = np.r_[B, -90, -90]
        path = Path(np.c_[L, B])
        hit = np.zeros(N, bool)
        for shift in (-360, 0, 360):         # test each point at lon and lon +/- 360 so wrapped rings work
            hit |= path.contains_points(np.c_[lon + shift, lat])
        inside ^= hit                         # even-odd rule: holes cancel out
    level += inside

np.save("milky_way_points.npy", np.c_[lon, lat, level])
print("points per level:", np.bincount(level))
