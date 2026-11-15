"""
Step 3a: recover 8 real SPHEREx images of M101 from a published figure.

Input : sources/m101_cutouts.png
        The 2x5 grid of cutouts in IRSA's "SPHEREx cutouts" tutorial
        (https://caltech-ipac.github.io/irsa-tutorials/spherex-cutouts/), BSD 3-Clause, Caltech.
        File: _public/4a2990c2668f392c5538e45e117da9b4.png on the gh-pages branch of
        https://github.com/Caltech-IPAC/irsa-tutorials
Output: m101_frames.json  (8 frames, 58x58 pixels, brightness 0-255, with wavelength)

The figure draws each image with matplotlib's "viridis" colours on a log scale.
We find each panel, map every colour back to its position on the viridis scale (0-255),
paint over the white wavelength labels, and shrink 3x3 screen pixels back to 1 SPHEREx pixel.
The two partial cutouts in the figure are skipped.
"""
import base64
import json
import matplotlib
import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.spatial import cKDTree

viridis = matplotlib.colormaps["viridis"](np.linspace(0, 1, 256))[:, :3] * 255
tree = cKDTree(viridis)
im = np.asarray(Image.open("sources/m101_cutouts.png").convert("RGB")).astype(float)

# (top, bottom, left, right) of each full panel in the PNG, and the wavelength printed on it (microns)
panels = [((10, 184, 51, 225), 1.506), ((10, 184, 292, 466), 1.543), ((10, 184, 532, 706), 1.58),
          ((10, 184, 773, 947), 1.618), ((220, 394, 292, 466), 1.113), ((220, 394, 532, 706), 1.139),
          ((220, 394, 773, 947), 1.166), ((220, 394, 1014, 1188), 1.194)]
frames = []
for (y0, y1, x0, x1), wl in panels:
    p = im[y0:y1, x0:x1]
    coloured = (p.max(2) - p.min(2)) > 25                     # viridis pixels are never grey
    dist, idx = tree.query(p.reshape(-1, 3))
    val = idx.reshape(p.shape[:2]).astype(float)
    bad = ~coloured | (dist.reshape(p.shape[:2]) > 30)        # label text and anti-aliasing
    if bad.any():
        _, (iy, ix) = ndimage.distance_transform_edt(bad, return_indices=True)
        val = val[iy, ix]                                     # fill from the nearest good pixel
    n = round(val.shape[1] / 3)
    blocks = val[:n * 3, :n * 3].reshape(n, 3, n, 3).transpose(0, 2, 1, 3).reshape(n, n, 9)
    img = np.median(blocks, axis=2)[::-1]                     # flip so row 0 is the bottom, like FITS
    frames.append({"wl": wl, "w": n, "h": n,
                   "px": base64.b64encode(np.clip(img, 0, 255).astype(np.uint8).tobytes()).decode()})
json.dump(frames, open("m101_frames.json", "w"))
print(len(frames), "frames")
