"""
Step 3a-2: cut clean, round images of the whole Pinwheel Galaxy (M101) out of two figures in
IRSA's SPHEREx mosaic tutorial, for the "Into the Pinwheel Galaxy" zoom animation.

  sources/mosaic.png  129ae5fa892b3d71d9a72deea3311404.png  all 102 wavelengths added (inferno, log)
  sources/pah.png     04203b8936167314fc90b300c63e981c.png  3.3 um dust (PAH) map     (inferno, linear)
  (gh-pages branch of https://github.com/Caltech-IPAC/irsa-tutorials, _public/ folder, BSD 3-Clause, Caltech)

The tutorial's mosaic is 267 x 267 pixels at 9 arcsec/px (40 arcmin), centred on M101
(RA 210.8022671, Dec 54.34895). Both figures put that image in the same box: x 65-450, y 10-395.
We keep a disc of radius 98 mosaic pixels (14.7 arcmin), just inside the outer background ring
and clear of the legend. On the mosaic figure, the white aperture rings (radius 20 and 60 px) are
painted over: white is not an inferno colour, so those pixels are easy to find, and they are
refilled from their neighbours. The dust map has nothing drawn on it and is left as it is.

Output: m101_disc.json  {"mosaic": data URI, "pah": data URI, "deg": width of each image in degrees}
Needs: numpy, scipy, pillow, matplotlib (only for the inferno colour table).
"""
import base64
import io
import json

import matplotlib
import numpy as np
from PIL import Image
from scipy import ndimage

BOX = (65, 10, 451, 396)            # left, top, right, bottom of the image area in both figures
MOSAIC_PX, ARCSEC_PX = 267, 9.0     # mosaic size and pixel scale (from the tutorial)
KEEP_R = 98                         # disc radius, in mosaic pixels
RINGS = (20, 60)                    # aperture rings drawn on the mosaic figure, in mosaic pixels

lut = matplotlib.colormaps["inferno"](np.linspace(0, 1, 256))[:, :3] * 255

def disc(path, rings):
    a = np.asarray(Image.open(path).convert("RGB").crop(BOX)).astype(float)
    n = a.shape[0]
    yy, xx = np.mgrid[0:n, 0:n]
    c = (n - 1) / 2
    r = np.hypot(xx - c, yy - c) * MOSAIC_PX / n            # radius in mosaic pixels
    if rings:
        flat = a.reshape(-1, 3)
        off = np.sqrt(((flat[:, None, :] - lut[None]) ** 2).sum(-1).min(1)).reshape(n, n)   # distance from the inferno colours
        near = np.zeros_like(off, bool)
        for rr in rings:
            near |= np.abs(r - rr) < 2.2
        bad = ndimage.binary_dilation((off > 18) & near, iterations=1) & near
        for ch in range(3):                                   # refill from neighbours, working inwards
            v, good = a[:, :, ch].copy(), ~bad
            while not good.all():
                s = ndimage.uniform_filter(np.where(good, v, 0), 3)
                w = ndimage.uniform_filter(good.astype(float), 3)
                f = ~good & (w > 0.01)
                v[f] = s[f] / w[f]
                good |= f
            a[:, :, ch] = v
    alpha = np.clip((KEEP_R - r) / 7, 0, 1)                   # soft edge
    img = Image.fromarray(np.dstack([a.clip(0, 255).astype(np.uint8), (alpha * 255).astype(np.uint8)]))
    R = int(np.ceil(KEEP_R * n / MOSAIC_PX)) + 1
    ci = int(round(c))
    img = img.crop((ci - R, ci - R, ci + R + 1, ci + R + 1))
    deg = img.width * (MOSAIC_PX / n) * ARCSEC_PX / 3600
    return img, deg

def webp(img):
    b = io.BytesIO()
    img.save(b, "WEBP", quality=90, alpha_quality=90, method=6)
    return "data:image/webp;base64," + base64.b64encode(b.getvalue()).decode()

mos, deg = disc("sources/mosaic.png", RINGS)
pah, _ = disc("sources/pah.png", ())
json.dump({"mosaic": webp(mos), "pah": webp(pah), "deg": round(deg, 5)}, open("m101_disc.json", "w"))
print(f"M101 discs: {mos.width} px, {deg:.4f} degrees across")
