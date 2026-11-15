"use strict";
// ================= constants =================
const PIX = 6.15;                       // SPHEREx pixel size, arcsec
const BANDS = [
  { n: 1, lo: 0.75, hi: 1.10, color: "#7fa8ff" }, { n: 2, lo: 1.10, hi: 1.62, color: "#69d0d8" },
  { n: 3, lo: 1.63, hi: 2.41, color: "#93d27c" }, { n: 4, lo: 2.42, hi: 3.82, color: "#e3c95c" },
  { n: 5, lo: 3.83, hi: 4.41, color: "#ee9a55" }, { n: 6, lo: 4.42, hi: 5.00, color: "#e8697c" }];
// Positions are J2000 at epoch 2000; proper motions in mas/yr (RA motion includes cos Dec).
const PRESETS = [
  { name: "Barnard's Star", ra: 269.452077, dec: 4.693365, pmra: -801.551, pmdec: 10362.394, note: "The fastest-moving star in our sky: about 10.4 arcseconds per year." },
  { name: "Epsilon Indi B", ra: 331.043833, dec: -56.782708, pmra: 3966.7, pmdec: -2536.2, note: "Two brown dwarfs orbiting each other, 12 light-years away." },
  { name: "WISE 0855−0714", ra: 133.7951321, dec: -7.24515, pmra: -8118.9, pmdec: 679.3, note: "The coldest brown dwarf known, about as cold as the North Pole. Faint in short bands, brightest near 4.5 µm (band 6)." },
  { name: "RW Aurigae", ra: 76.956068, dec: 30.401326, pmra: 4.9, pmdec: -24.6, note: "A young star that fades when clumps of dust pass in front of it." },
  { name: "R Draconis", ra: 248.16768, dec: 66.754984, pmra: -20.07, pmdec: 7.06, note: "A pulsating red giant (Mira variable) that brightens and fades over about 245 days." }];
const $ = id => document.getElementById(id);
const d2r = Math.PI / 180;

// ================= state =================
const S = {
  target: null, rows: [], frames: [], items: [], idx: 0, playing: false, run: 0,
  band: 0, framesMode: "visit", show: "image", gaia: [], grid: null, cache: new Map(), chartMode: "time"
};

// ================= coordinates =================
function parseCoords(txt) {
  const t = txt.trim();
  if (!/^[-+\d\s.:hdms°'′"″,]+$/i.test(t) || !/\d/.test(t)) return null;
  const sexa = /[hms:°'′"″]/i.test(t) || t.replace(/,/g, " ").trim().split(/\s+/).length === 6;
  const toks = t.replace(/[hdms°'′"″:,]/gi, " ").trim().split(/\s+/);
  if (!sexa && toks.length === 2) {
    const ra = +toks[0], dec = +toks[1];
    return isFinite(ra) && isFinite(dec) && Math.abs(dec) <= 90 ? { ra: ((ra % 360) + 360) % 360, dec } : null;
  }
  if (toks.length === 6) {
    const [h, m, s, dd, dm, ds] = toks.map(Number);
    if ([h, m, s, dd, dm, ds].some(v => !isFinite(v))) return null;
    const sign = toks[3].startsWith("-") ? -1 : 1;
    const ra = (h + m / 60 + s / 3600) * 15, dec = sign * (Math.abs(dd) + dm / 60 + ds / 3600);
    return Math.abs(dec) <= 90 ? { ra, dec } : null;
  }
  return null;
}
const fmtRA = ra => { const h = ra / 15, H = Math.floor(h), M = Math.floor((h - H) * 60), Sx = ((h - H) * 60 - M) * 60;
  return `${String(H).padStart(2, "0")}h${String(M).padStart(2, "0")}m${Sx.toFixed(2).padStart(5, "0")}s`; };
const fmtDec = dec => { const a = Math.abs(dec), D = Math.floor(a), M = Math.floor((a - D) * 60), Sx = ((a - D) * 60 - M) * 60;
  return `${dec < 0 ? "−" : "+"}${String(D).padStart(2, "0")}°${String(M).padStart(2, "0")}′${Sx.toFixed(1).padStart(4, "0")}″`; };
const mjdToYear = mjd => 2000 + (mjd - 51544.5) / 365.25;
const mjdToDate = mjd => new Date((mjd - 40587) * 864e5);
function targetAt(year) {           // target position at a given year, following its proper motion
  const t = S.target, dy = year - (t.epoch || 2000);
  const dec = t.dec + t.pmdec * dy / 3.6e6;
  return { ra: t.ra + t.pmra * dy / 3.6e6 / Math.cos(t.dec * d2r), dec };
}

// ================= FITS & WCS =================
function parseCardValue(s) {
  s = s.trim();
  if (s.startsWith("'")) { const e = s.indexOf("'", 1); return s.slice(1, e < 0 ? undefined : e).trim(); }
  const v = s.split("/")[0].trim();
  if (v === "T") return true; if (v === "F") return false;
  const n = parseFloat(v.replace("D", "E")); return isNaN(n) ? v : n;
}
function parseFITS(buf) {
  const u8 = new Uint8Array(buf), hdus = []; let off = 0;
  while (off + 2880 <= buf.byteLength) {
    const hdr = {}; let end = false;
    while (!end) {
      if (off + 2880 > buf.byteLength) return hdus;
      for (let k = 0; k < 36; k++) {
        const card = String.fromCharCode.apply(null, u8.subarray(off + k * 80, off + k * 80 + 80));
        const key = card.slice(0, 8).trim();
        if (key === "END") { end = true; break; }
        if (card[8] === "=" && card[9] === " ") hdr[key] = parseCardValue(card.slice(10));
      }
      off += 2880;
    }
    const naxis = hdr.NAXIS || 0; let n = naxis ? 1 : 0;
    for (let a = 1; a <= naxis; a++) n *= hdr["NAXIS" + a];
    const size = Math.abs(hdr.BITPIX || 8) / 8 * (hdr.GCOUNT || 1) * ((hdr.PCOUNT || 0) + n);
    hdus.push({ hdr, off, size });
    off += Math.ceil(size / 2880) * 2880;
  }
  return hdus;
}
function readImage(buf, h) {
  const dv = new DataView(buf), w = h.hdr.NAXIS1, ht = h.hdr.NAXIS2, n = w * ht, out = new Float32Array(n);
  const bp = h.hdr.BITPIX, bs = h.hdr.BSCALE ?? 1, bz = h.hdr.BZERO ?? 0, o = h.off;
  for (let i = 0; i < n; i++) {
    let v;
    if (bp === -32) v = dv.getFloat32(o + i * 4, false); else if (bp === -64) v = dv.getFloat64(o + i * 8, false);
    else if (bp === 16) v = dv.getInt16(o + i * 2, false); else if (bp === 32) v = dv.getInt32(o + i * 4, false); else v = dv.getUint8(o + i);
    out[i] = v * bs + bz;
  }
  return { w, h: ht, data: out };
}
function makeWCS(h) {
  if (!h.CTYPE1 || !String(h.CTYPE1).includes("TAN") || h.CRVAL1 === undefined) return null;
  let cd;
  if (h.CD1_1 !== undefined) cd = [h.CD1_1, h.CD1_2 || 0, h.CD2_1 || 0, h.CD2_2];
  else { const c1 = h.CDELT1 ?? 1, c2 = h.CDELT2 ?? 1; cd = [(h.PC1_1 ?? 1) * c1, (h.PC1_2 || 0) * c1, (h.PC2_1 || 0) * c2, (h.PC2_2 ?? 1) * c2]; }
  const det = cd[0] * cd[3] - cd[1] * cd[2], inv = [cd[3] / det, -cd[1] / det, -cd[2] / det, cd[0] / det];
  const terms = p => { const t = [], ord = h[p + "_ORDER"] || 0;
    for (let i = 0; i <= ord; i++) for (let j = 0; j <= ord - i; j++) { const c = h[`${p}_${i}_${j}`]; if (c) t.push([i, j, c]); } return t; };
  const A = terms("A"), B = terms("B");
  const poly = (t, u, v) => { let s = 0; for (const [i, j, c] of t) s += c * Math.pow(u, i) * Math.pow(v, j); return s; };
  const ra0 = h.CRVAL1 * d2r, de0 = h.CRVAL2 * d2r, sd0 = Math.sin(de0), cd0 = Math.cos(de0);
  return {
    pix2sky(x, y) {           // 1-based FITS pixel -> radians
      const u = x - h.CRPIX1, v = y - h.CRPIX2, uu = u + poly(A, u, v), vv = v + poly(B, u, v);
      const xi = (cd[0] * uu + cd[1] * vv) * d2r, eta = (cd[2] * uu + cd[3] * vv) * d2r, den = cd0 - eta * sd0;
      return [ra0 + Math.atan2(xi, den), Math.atan((sd0 + eta * cd0) / Math.hypot(xi, den))];
    },
    sky2pix(ra, dec) {        // radians -> 1-based FITS pixel
      const cosd = Math.cos(dec), sind = Math.sin(dec), dra = ra - ra0, cdra = Math.cos(dra);
      const cc = sd0 * sind + cd0 * cosd * cdra;
      const xi = cosd * Math.sin(dra) / cc / d2r, eta = (cd0 * sind - sd0 * cosd * cdra) / cc / d2r;
      const U0 = inv[0] * xi + inv[1] * eta, V0 = inv[2] * xi + inv[3] * eta; let u = U0, v = V0;
      if (A.length || B.length) for (let k = 0; k < 8; k++) { const nu = U0 - poly(A, u, v), nv = V0 - poly(B, u, v); u = nu; v = nv; }
      return [u + h.CRPIX1, v + h.CRPIX2];
    }
  };
}
function headerTime(h) {
  for (const k of ["MJD-AVG", "MJD-MID", "MJD_AVG", "MJDMID", "MJD-OBS"]) if (typeof h[k] === "number") return h[k];
  if (typeof h["MJD-BEG"] === "number") return typeof h["MJD-END"] === "number" ? (h["MJD-BEG"] + h["MJD-END"]) / 2 : h["MJD-BEG"];
  for (const k of ["DATE-AVG", "DATE-MID", "DATE-OBS", "DATE-BEG"]) if (typeof h[k] === "string" && /^\d{4}-/.test(h[k])) return Date.parse(h[k] + (h[k].endsWith("Z") ? "" : "Z")) / 864e5 + 40587;
  return null;
}

// ================= output grid (north up, east left) =================
function makeGrid(ra0, dec0, fov) {
  const N = Math.max(16, Math.min(200, Math.round(fov * 3600 / PIX))), s = fov / N, c = (N - 1) / 2;
  const RA = new Float64Array(N * N), DE = new Float64Array(N * N), a0 = ra0 * d2r, b0 = dec0 * d2r;
  for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
    const xi = -(i - c) * s * d2r, eta = (j - c) * s * d2r, den = Math.cos(b0) - eta * Math.sin(b0);
    RA[j * N + i] = a0 + Math.atan2(xi, den); DE[j * N + i] = Math.atan((Math.sin(b0) + eta * Math.cos(b0)) / Math.hypot(xi, den));
  }
  return { N, s, c, ra0, dec0, fov, RA, DE,
    toGrid(ra, dec) {        // degrees -> grid pixel (0-based, j up)
      const a = ra * d2r, b = dec * d2r, cc = Math.sin(b0) * Math.sin(b) + Math.cos(b0) * Math.cos(b) * Math.cos(a - a0);
      const xi = Math.cos(b) * Math.sin(a - a0) / cc / d2r, eta = (Math.cos(b0) * Math.sin(b) - Math.sin(b0) * Math.cos(b) * Math.cos(a - a0)) / cc / d2r;
      return [c - xi / s, c + eta / s];
    } };
}
function reproject(img, wcs, G) {
  const out = new Float32Array(G.N * G.N).fill(NaN), w = img.w, h = img.h, d = img.data;
  for (let k = 0; k < out.length; k++) {
    let [x, y] = wcs.sky2pix(G.RA[k], G.DE[k]); x -= 1; y -= 1;
    const x0 = Math.floor(x), y0 = Math.floor(y);
    if (x0 < 0 || y0 < 0 || x0 >= w - 1 || y0 >= h - 1) continue;
    const fx = x - x0, fy = y - y0, o = y0 * w + x0;
    out[k] = d[o] * (1 - fx) * (1 - fy) + d[o + 1] * fx * (1 - fy) + d[o + w] * (1 - fx) * fy + d[o + w + 1] * fx * fy;
  }
  return out;
}

// ================= statistics =================
function finiteSorted(a) { const t = []; for (const v of a) if (Number.isFinite(v)) t.push(v); return t.sort((x, y) => x - y); }
const pct = (sorted, p) => sorted.length ? sorted[Math.min(sorted.length - 1, Math.max(0, Math.floor(p / 100 * sorted.length)))] : 0;
function robust(a) { const s = finiteSorted(a), med = pct(s, 50); const dev = finiteSorted(s.map(v => Math.abs(v - med))); return { med, sig: (pct(dev, 50) * 1.4826) || 1 }; }
function photometry(arr, G, x, y) {
  let ap = 0, n = 0; const ann = [];
  for (let j = Math.floor(y - 8); j <= y + 8; j++) for (let i = Math.floor(x - 8); i <= x + 8; i++) {
    if (i < 0 || j < 0 || i >= G.N || j >= G.N) continue;
    const v = arr[j * G.N + i]; if (!Number.isFinite(v)) continue;
    const r = Math.hypot(i - x, j - y);
    if (r <= 1.8) { ap += v; n++; } else if (r >= 4 && r <= 7.5) ann.push(v);
  }
  if (n < 5 || ann.length < 10) return null;
  ann.sort((p, q) => p - q);
  return ap - ann[ann.length >> 1] * n;
}

// ================= bands =================
function bandOf(row) {
  let mid = (row.em_min + row.em_max) / 2;
  if (Number.isFinite(mid) && mid > 0) {
    if (mid < 1e-3) mid *= 1e6;                       // metres -> microns
    let best = 1, bd = 1e9;
    for (const b of BANDS) { const d = mid < b.lo ? b.lo - mid : mid > b.hi ? mid - b.hi : 0; if (d < bd) { bd = d; best = b.n; } }
    return best;
  }
  const m = /D\s*([1-6])\b/i.exec(row.band || "") || /([1-6])\s*$/.exec(row.band || "");
  return m ? +m[1] : 0;
}

// ================= API =================
async function api(path, params) {
  const r = await fetch(path + "?" + new URLSearchParams(params));
  const type = r.headers.get("content-type") || "";
  if (type.includes("json")) { const j = await r.json(); if (!r.ok || j.error) throw new Error(j.error || "Request failed"); return j; }
  if (!r.ok) throw new Error("Request failed (HTTP " + r.status + ")");
  return r.arrayBuffer();
}

// ================= main flow =================
function setStatus(title, text, progress) {
  $("status").hidden = !title; $("statusTitle").textContent = title || ""; $("statusText").textContent = text || "";
  $("bar").hidden = progress === undefined; if (progress !== undefined) $("barFill").style.width = Math.round(progress * 100) + "%";
}
async function start(target) {
  const run = ++S.run; stop();
  S.target = target; S.frames = []; S.items = []; S.rows = []; S.gaia = []; S.idx = 0; S.cache.clear();
  document.querySelectorAll("#presets .chip").forEach(c => c.classList.toggle("on", c.textContent === target.name));
  $("target").value = target.name && !target.typed ? target.name : `${fmtRA(target.ra)} ${fmtDec(target.dec)}`;
  renderTargetCard(); updateURL(); renderAll();
  const fov = +$("fov").value, center = targetAt(2025.9);
  S.grid = makeGrid(center.ra, center.dec, fov);
  S.sim = false; $("simTag").hidden = !DEMO_BACKEND; $("simTag").textContent = "Demo data (synthetic)";
  if (OFFLINE) return startDemo(target);
  $("goBtn").disabled = true;
  try {
    setStatus("Searching IRSA…", "Asking NASA's infrared archive for every SPHEREx image of this spot.", 0);
    const res = await api("/api/search", { ra: center.ra.toFixed(6), dec: center.dec.toFixed(6), collection: $("release").value });
    if (run !== S.run) return;
    S.rows = res.rows.filter(r => r.t_min).map(r => ({ ...r, bandN: bandOf(r) })).sort((a, b) => a.t_min - b.t_min);
    renderBands();
    if (!S.rows.length) { setStatus("No SPHEREx images here yet", "This release has no images covering this position. Try another data release, or check back as new weeks are added."); return; }
    api("/api/gaia", { ra: center.ra.toFixed(6), dec: center.dec.toFixed(6), r: (fov * 0.75).toFixed(3) })
      .then(g => { if (run === S.run) { S.gaia = g.rows; draw(); } }).catch(() => {});
    renderTargetCard();
    await loadFrames(run);
  } catch (e) {
    if (run === S.run) setStatus("Something went wrong", e.message === "Failed to fetch" ? "The page can't reach its server. Start it with: python server.py, then open http://localhost:8000" : e.message);
  } finally { if (run === S.run) $("goBtn").disabled = false; }
}
function pickRows() {
  const max = +$("maxn").value;
  let rows = S.rows.filter(r => S.band === 0 || r.bandN === S.band);
  if (rows.length > max) { const out = []; for (let k = 0; k < max; k++) out.push(rows[Math.round(k * (rows.length - 1) / (max - 1))]); rows = [...new Set(out)]; }
  return rows;
}
async function loadFrames(run) {
  const want = pickRows().filter(r => !S.frames.some(f => f.row.access_url === r.access_url));
  if (!want.length) { rebuild(); setStatus(null); return; }
  let done = 0, failed = 0, lastErr = "";
  const center = { ra: S.grid.ra0, dec: S.grid.dec0 }, size = Math.min(0.5, S.grid.fov * 1.45 + 0.01);   // covers the grid corners at any rotation
  setStatus("Downloading images…", `0 of ${want.length} cutouts`, 0);
  let next = 0, lastDraw = 0;
  async function worker() {
    while (next < want.length) {
      const row = want[next++];
      try {
        const buf = await api("/api/cutout", { url: row.access_url, ra: center.ra.toFixed(7), dec: center.dec.toFixed(7), size: size.toFixed(4) });
        if (run !== S.run) return;
        const f = buildFrame(buf, row); if (f) S.frames.push(f); else failed++;
      } catch (e) { failed++; lastErr = e.message; }
      done++;
      if (run !== S.run) return;
      setStatus("Downloading images…", `${done} of ${want.length} cutouts${failed ? ` · ${failed} failed` : ""}`, done / want.length);
      if (performance.now() - lastDraw > 700) { lastDraw = performance.now(); rebuild(); }
    }
  }
  await Promise.all([worker(), worker(), worker(), worker()]);
  if (run !== S.run) return;
  rebuild();
  if (!S.frames.length) setStatus("Couldn't load any images", lastErr || "The archive didn't return usable cutouts. Try again in a minute.");
  else setStatus(null);
}
function buildFrame(buf, row) {
  const hdus = parseFITS(buf);
  const imgs = hdus.filter(x => x.hdr.NAXIS === 2 && x.hdr.XTENSION !== "BINTABLE" && x.hdr.NAXIS1 > 1);
  const ih = imgs.find(x => x.hdr.EXTNAME === "IMAGE") || imgs[0];
  if (!ih) return null;
  const img = readImage(buf, ih), wcs = makeWCS(ih.hdr);
  if (!wcs) return null;
  const zh = imgs.find(x => x.hdr.EXTNAME === "ZODI" && x.hdr.NAXIS1 === img.w && x.hdr.NAXIS2 === img.h);
  if (zh) { const z = readImage(buf, zh).data; for (let i = 0; i < img.data.length; i++) img.data[i] -= z[i]; }
  const data = reproject(img, wcs, S.grid);
  const { med, sig } = robust(data), norm = new Float32Array(data.length);
  for (let i = 0; i < data.length; i++) norm[i] = (data[i] - med) / sig;
  const mjd = row.t_min && row.t_max ? (row.t_min + row.t_max) / 2 : row.t_min || headerTime(ih.hdr) || headerTime(hdus[0].hdr);
  const year = mjdToYear(mjd), tp = targetAt(year), [tx, ty] = S.grid.toGrid(tp.ra, tp.dec);
  return { row, mjd, year, band: row.bandN, data, norm, flux: photometry(data, S.grid, tx, ty), tx, ty };
}

// ================= frame list & display =================
function rebuild() {
  const frames = S.frames.filter(f => S.band === 0 || f.band === S.band).sort((a, b) => a.mjd - b.mjd);
  let items;
  if (S.framesMode === "each") items = frames.map(f => ({ mjd: f.mjd, year: f.year, data: f.data, norm: f.norm, frames: [f] }));
  else {
    const groups = [];
    for (const f of frames) { const g = groups[groups.length - 1]; if (g && f.mjd - g[g.length - 1].mjd < 20) g.push(f); else groups.push([f]); }
    items = groups.map(g => {
      const n = g[0].norm.length, sum = new Float32Array(n), cnt = new Uint16Array(n);
      for (const f of g) for (let i = 0; i < n; i++) { const v = f.norm[i]; if (Number.isFinite(v)) { sum[i] += v; cnt[i]++; } }
      for (let i = 0; i < n; i++) sum[i] = cnt[i] ? sum[i] / Math.sqrt(cnt[i]) : NaN;     // signal-to-noise of the stack
      const mjd = g.reduce((s, f) => s + f.mjd, 0) / g.length;
      return { mjd, year: mjdToYear(mjd), data: sum, norm: sum, frames: g };
    });
  }
  S.items = items; S.cache.clear(); S.median = null;
  S.idx = Math.min(S.idx, Math.max(0, items.length - 1));
  $("frameSlider").max = Math.max(0, items.length - 1);
  renderAll();
}
function medianOf(items) {
  if (S.median) return S.median;
  const n = items[0].norm.length, out = new Float32Array(n), col = [];
  for (let i = 0; i < n; i++) { col.length = 0; for (const it of items) { const v = it.norm[i]; if (Number.isFinite(v)) col.push(v); }
    col.sort((a, b) => a - b); out[i] = col.length ? col[col.length >> 1] : NaN; }
  return (S.median = out);
}
function displayArray(k) {
  const key = S.show + ":" + k; if (S.cache.has(key)) return S.cache.get(key);
  const it = S.items[k]; let arr = it.data, kind = "gray";
  if (S.show !== "image" && S.items.length > 1) {
    const ref = S.show === "median" ? medianOf(S.items) : S.items[0].norm;
    arr = new Float32Array(it.norm.length); for (let i = 0; i < arr.length; i++) arr[i] = it.norm[i] - ref[i]; kind = "diff";
  }
  const sorted = finiteSorted(arr), out = { arr, kind, sorted, sig: kind === "diff" ? robust(arr).sig : 1 };
  S.cache.set(key, out); return out;
}

// ================= drawing =================
const view = $("view"), vctx = view.getContext("2d"), off = document.createElement("canvas"), octx = off.getContext("2d");
function sizeView() { const r = view.getBoundingClientRect(), dpr = Math.min(2, devicePixelRatio || 1); view.width = Math.round(r.width * dpr); view.height = Math.round(r.height * dpr); vctx.setTransform(dpr, 0, 0, dpr, 0, 0); draw(); }
new ResizeObserver(sizeView).observe($("viewer"));
function draw() {
  const W = view.clientWidth, H = view.clientHeight;
  vctx.fillStyle = "#000"; vctx.fillRect(0, 0, W, H);
  if (!S.items.length || !S.grid) { drawFrameInfo(); return; }
  const G = S.grid, N = G.N, it = S.items[S.idx], D = displayArray(S.idx), c = +$("contrast").value;
  if (off.width !== N) { off.width = N; off.height = N; }
  const img = octx.createImageData(N, N), px = img.data;
  if (D.kind === "gray") {
    const lo = pct(D.sorted, 30), hi = pct(D.sorted, 96 + c * 0.039), ash = Math.asinh(8), span = hi - lo || 1;
    for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
      const v = D.arr[j * N + i], o = ((N - 1 - j) * N + i) * 4;
      const t = Number.isFinite(v) ? Math.min(1, Math.asinh(8 * Math.max(0, (v - lo) / span)) / ash) : 0;
      px[o] = px[o + 1] = px[o + 2] = 255 * t; px[o + 3] = 255;
    }
  } else {
    const floor = D.sig * 2, scale = D.sig * (12 - c * 0.1);       // changes under 2 sigma are treated as noise
    for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
      const v = D.arr[j * N + i], o = ((N - 1 - j) * N + i) * 4;
      const t = Number.isFinite(v) ? Math.sign(v) * Math.min(1, Math.max(0, Math.abs(v) - floor) / scale) : 0;
      if (t >= 0) { px[o] = 121 * t; px[o + 1] = 198 * t; px[o + 2] = 234 * t; } else { px[o] = 235 * -t; px[o + 1] = 154 * -t; px[o + 2] = 92 * -t; }
      px[o + 3] = 255;
    }
  }
  octx.putImageData(img, 0, 0);
  vctx.imageSmoothingEnabled = $("smooth").checked;
  vctx.drawImage(off, 0, 0, W, H);
  const k = W / N, toCanvas = (gi, gj) => [(gi + 0.5) * k, (N - 0.5 - gj) * k];
  // Gaia stars at this date
  if ($("ovGaia").checked && S.gaia.length) {
    const dy = it.year - 2016; vctx.lineWidth = 1;
    for (const s of S.gaia) {
      if (s.g > 19.5) continue;
      const dec = s.dec + s.pmdec * dy / 3.6e6, ra = s.ra + s.pmra * dy / 3.6e6 / Math.cos(s.dec * d2r);
      const [gi, gj] = G.toGrid(ra, dec); if (gi < -2 || gj < -2 || gi > N + 1 || gj > N + 1) continue;
      const [x, y] = toCanvas(gi, gj), r = Math.max(2.5, Math.min(9, (20 - s.g) * 0.7)) * Math.max(1, k / 8);
      vctx.strokeStyle = Math.hypot(s.pmra, s.pmdec) > 200 ? "rgba(124,200,255,0.95)" : "rgba(121,198,234,0.55)";
      vctx.beginPath(); vctx.arc(x, y, r, 0, 7); vctx.stroke();
    }
  }
  // predicted target position
  if ($("ovTarget").checked && S.target) {
    const tp = targetAt(it.year), [gi, gj] = G.toGrid(tp.ra, tp.dec), [x, y] = toCanvas(gi, gj), a = Math.max(10, k * 2.2);
    vctx.strokeStyle = "#7cc8ff"; vctx.lineWidth = 1.6; vctx.beginPath();
    vctx.moveTo(x - a * 2, y); vctx.lineTo(x - a, y); vctx.moveTo(x + a, y); vctx.lineTo(x + a * 2, y);
    vctx.moveTo(x, y - a * 2); vctx.lineTo(x, y - a); vctx.moveTo(x, y + a); vctx.lineTo(x, y + a * 2); vctx.stroke();
  }
  // date, compass, scale bar
  vctx.font = "500 14px 'IBM Plex Mono', monospace"; vctx.fillStyle = "rgba(0,0,0,0.55)"; vctx.fillRect(10, 10, 210, 26);
  vctx.fillStyle = "#e6f0fb"; vctx.fillText(mjdToDate(it.mjd).toISOString().slice(0, 10) + "  " + labelFor(it, true), 18, 28);
  const arcmin = 60 / (G.s * 3600) * k, bx = 16, by = H - 18;
  vctx.strokeStyle = "#e6f0fb"; vctx.lineWidth = 2; vctx.beginPath(); vctx.moveTo(bx, by); vctx.lineTo(bx + arcmin, by); vctx.stroke();
  vctx.font = "12px 'IBM Plex Sans', sans-serif"; vctx.fillText("1′", bx + arcmin / 2 - 5, by - 6);
  const cx = W - 34, cy = H - 34; vctx.lineWidth = 1.4; vctx.beginPath(); vctx.moveTo(cx, cy); vctx.lineTo(cx, cy - 22); vctx.moveTo(cx, cy); vctx.lineTo(cx - 22, cy); vctx.stroke();
  vctx.fillText("N", cx - 4, cy - 26); vctx.fillText("E", cx - 36, cy + 4);
  drawFrameInfo();
}
function labelFor(it, short) {
  const bands = [...new Set(it.frames.map(f => f.band))].sort();
  if (it.frames.length === 1) { const b = BANDS[bands[0] - 1]; return !b ? "band unknown" : short ? `band ${b.n}` : `Band ${b.n} · ${b.lo}–${b.hi} µm`; }
  return short ? `${it.frames.length} images` : `${it.frames.length} images stacked · band${bands.length > 1 ? "s " + bands.join(", ") : " " + bands[0]}`;
}
function drawFrameInfo() {
  const n = S.items.length, it = S.items[S.idx];
  $("counter").textContent = n ? `${S.idx + 1} / ${n}` : "0 / 0"; $("frameSlider").value = S.idx;
  if (!it) { $("frameInfo").innerHTML = ""; return; }
  const d = mjdToDate(it.mjd).toISOString().slice(0, 16).replace("T", " ") + " UTC";
  const one = it.frames.length === 1 ? it.frames[0].row : null;
  $("frameInfo").innerHTML = `<span><span class="k">Date</span><span class="v">${d}</span></span>
    <span><span class="k">Data</span><span class="v">${labelFor(it)}</span></span>
    ${one ? `<span><span class="k">Image</span><span class="v">${one.obs_id || "—"}</span> · <a href="${one.access_url}">full file</a></span>` : ""}
    ${n > 1 ? `<span><span class="k">Since first</span><span class="v">${Math.round(it.mjd - S.items[0].mjd)} days</span></span>` : ""}`;
}

// ================= chart =================
const chart = $("chart"), cctx = chart.getContext("2d"); let chartPts = [];
function drawChart() {
  const r = chart.getBoundingClientRect(), dpr = Math.min(2, devicePixelRatio || 1);
  chart.width = Math.round(r.width * dpr); chart.height = Math.round(r.height * dpr); cctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  const W = r.width, H = r.height, L = 44, R = 10, T = 10, B = 26;
  cctx.clearRect(0, 0, W, H); chartPts = [];
  const pts = S.frames.filter(f => f.flux !== null && f.band > 0 && (S.band === 0 || f.band === S.band));
  cctx.font = "11px 'IBM Plex Mono', monospace"; cctx.fillStyle = "#8fa0bd";
  if (pts.length < 2) { cctx.fillText(S.frames.length ? "Not enough measurements yet" : "Load a target to see its brightness", L, H / 2); return; }
  const ys = finiteSorted(pts.map(f => f.flux));
  let y0 = pct(ys, 2), y1 = pct(ys, 98); const pad = (y1 - y0) * 0.12 || 1; y0 -= pad; y1 += pad;
  const time = S.chartMode === "time";
  let x0 = time ? Math.min(...pts.map(f => f.mjd)) : 0.7, x1 = time ? Math.max(...pts.map(f => f.mjd)) : 5.05; if (x1 === x0) x1 = x0 + 1;
  const X = v => L + (v - x0) / (x1 - x0) * (W - L - R), Y = v => T + (1 - (Math.max(y0, Math.min(y1, v)) - y0) / (y1 - y0)) * (H - T - B);
  cctx.strokeStyle = "rgba(124,200,255,0.12)"; cctx.lineWidth = 1;
  for (let k = 0; k <= 4; k++) { const v = y0 + (y1 - y0) * k / 4, y = Y(v); cctx.beginPath(); cctx.moveTo(L, y); cctx.lineTo(W - R, y); cctx.stroke(); cctx.fillText(v.toPrecision(2), 2, y + 4); }
  if (time) {
    const e0 = mjdToDate(x0).toISOString().slice(0, 7), e1 = mjdToDate(x1).toISOString().slice(0, 7);
    cctx.fillText(e0, L, H - 8); cctx.fillText(e1, W - R - cctx.measureText(e1).width, H - 8);
    const cur = S.items[S.idx]; if (cur) { const x = X(cur.mjd); cctx.strokeStyle = "rgba(235,228,211,0.5)"; cctx.beginPath(); cctx.moveTo(x, T); cctx.lineTo(x, H - B); cctx.stroke(); }
    const curFrames = new Set((S.items[S.idx] || { frames: [] }).frames);
    for (const f of pts) {
      const x = X(f.mjd), y = Y(f.flux), on = curFrames.has(f);
      cctx.fillStyle = BANDS[f.band - 1].color; cctx.globalAlpha = on ? 1 : 0.75;
      cctx.beginPath(); cctx.arc(x, y, on ? 4.5 : 2.8, 0, 7); cctx.fill();
      if (on) { cctx.strokeStyle = "#e6f0fb"; cctx.stroke(); }
      chartPts.push([x, y, f]);
    }
    cctx.globalAlpha = 1;
  } else {
    for (const w of [1, 2, 3, 4, 5]) cctx.fillText(w + " µm", X(w) - 12, H - 8);
    for (const b of BANDS) {                 // median and 10-90% spread of each band's measurements
      const v = finiteSorted(pts.filter(f => f.band === b.n).map(f => f.flux)); if (!v.length) continue;
      const x = X((b.lo + b.hi) / 2), half = (X(b.hi) - X(b.lo)) / 2;
      cctx.strokeStyle = b.color; cctx.fillStyle = b.color; cctx.lineWidth = 2;
      cctx.beginPath(); cctx.moveTo(x, Y(pct(v, 10))); cctx.lineTo(x, Y(pct(v, 90))); cctx.stroke();
      cctx.lineWidth = 1; cctx.beginPath(); cctx.moveTo(x - half, Y(pct(v, 50))); cctx.lineTo(x + half, Y(pct(v, 50))); cctx.stroke();
      cctx.beginPath(); cctx.arc(x, Y(pct(v, 50)), 4, 0, 7); cctx.fill();
    }
  }
}
chart.addEventListener("click", e => {
  const r = chart.getBoundingClientRect(), x = e.clientX - r.left, y = e.clientY - r.top;
  let best = null, bd = 14; for (const [px, py, f] of chartPts) { const d = Math.hypot(px - x, py - y); if (d < bd) { bd = d; best = f; } }
  if (best) { const k = S.items.findIndex(it => it.frames.includes(best)); if (k >= 0) { S.idx = k; renderAll(); } }
});
$("legend").innerHTML = BANDS.map(b => `<span><i style="background:${b.color}"></i>Band ${b.n} · ${b.lo}–${b.hi} µm</span>`).join("");

// ================= side panels =================
function renderTargetCard() {
  const t = S.target; if (!t) return;
  $("tName").textContent = t.name || "Custom position";
  $("tNote").textContent = (OFFLINE && t.name !== "Your files" ? "Simulated preview. " : "") + (t.note || (t.typed ? "Typed coordinates. No proper motion is applied." : "Position and motion from CDS Sesame."));
  const now = targetAt(2025.9), pm = Math.hypot(t.pmra, t.pmdec);
  let span = "";
  if (S.rows.length > 1) {
    const yrs = (S.rows[S.rows.length - 1].t_min - S.rows[0].t_min) / 365.25, mv = pm / 1000 * yrs;
    span = `<span class="k">Expected shift</span><span class="v">${mv.toFixed(1)}″ (${(mv / PIX).toFixed(1)} px) over ${yrs.toFixed(1)} yr</span>`;
  }
  $("tFacts").innerHTML = `<span class="k">RA (2025.9)</span><span class="v">${fmtRA(now.ra)}</span>
    <span class="k">Dec (2025.9)</span><span class="v">${fmtDec(now.dec)}</span>
    <span class="k">Motion</span><span class="v">${pm >= 1 ? (pm / 1000).toFixed(pm > 1000 ? 2 : 3) + "″ / year" : "none applied"}</span>
    <span class="k">Images found</span><span class="v">${S.rows.length || "—"}</span>${span}
    <span class="k">Context</span><span class="v"><a href="starmap.html?ra=${now.ra.toFixed(4)}&dec=${now.dec.toFixed(4)}">see it on the star map</a></span>`;
}
function renderBands() {
  const counts = [S.rows.length, ...BANDS.map(b => S.rows.filter(r => r.bandN === b.n).length)];
  $("bands").innerHTML = ["All", 1, 2, 3, 4, 5, 6].map((b, k) =>
    `<button data-v="${k}" class="${S.band === k ? "on" : ""}" ${S.rows.length && !counts[k] ? "disabled" : ""} title="${k ? `Band ${k}: ${BANDS[k - 1].lo}–${BANDS[k - 1].hi} µm` : "All six bands"}">${b}<small>${S.rows.length ? counts[k] : "&nbsp;"}</small></button>`).join("");
  $("bands").querySelectorAll("button").forEach(btn => btn.onclick = async () => {
    S.band = +btn.dataset.v; renderBands(); updateURL(); rebuild();
    if (S.rows.length && !OFFLINE && S.target && S.target.name !== "Your files") { const run = S.run; $("goBtn").disabled = true; await loadFrames(run); if (run === S.run) $("goBtn").disabled = false; }
  });
}
function renderAll() { draw(); drawChart(); }
window.addEventListener("resize", drawChart);

// ================= player =================
let timer = null;
function stop() { S.playing = false; clearInterval(timer); $("playBtn").textContent = "▶"; $("playBtn").setAttribute("aria-label", "Play"); }
function play() {
  if (S.items.length < 2) return; S.playing = true; $("playBtn").textContent = "❚❚"; $("playBtn").setAttribute("aria-label", "Pause");
  clearInterval(timer); timer = setInterval(() => { S.idx = (S.idx + 1) % S.items.length; renderAll(); }, 1000 / +$("fps").value);
}
$("playBtn").onclick = () => S.playing ? stop() : play();
$("prevBtn").onclick = () => { stop(); if (S.items.length) { S.idx = (S.idx - 1 + S.items.length) % S.items.length; renderAll(); } };
$("nextBtn").onclick = () => { stop(); if (S.items.length) { S.idx = (S.idx + 1) % S.items.length; renderAll(); } };
$("frameSlider").oninput = e => { stop(); S.idx = +e.target.value; renderAll(); };
$("fps").onchange = () => { if (S.playing) play(); };
window.addEventListener("keydown", e => {
  if (["INPUT", "SELECT"].includes(e.target.tagName) && e.target.type !== "range") return;
  if (e.code === "Space") { e.preventDefault(); S.playing ? stop() : play(); }
  if (e.key === "ArrowRight") $("nextBtn").click(); if (e.key === "ArrowLeft") $("prevBtn").click();
});

// ================= controls =================
function seg(id, key, after) {
  $(id).querySelectorAll("button").forEach(b => b.onclick = () => {
    S[key] = b.dataset.v; $(id).querySelectorAll("button").forEach(x => x.classList.toggle("on", x === b)); updateURL(); after();
  });
}
seg("framesMode", "framesMode", rebuild);
seg("showMode", "show", () => { S.cache.clear(); renderAll(); });
seg("chartMode", "chartMode", drawChart);
["contrast", "ovGaia", "ovTarget", "smooth"].forEach(id => $(id).addEventListener("input", draw));
["fov", "release", "maxn"].forEach(id => $(id).addEventListener("change", () => { updateURL(); if (S.target && S.target.name !== "Your files") start(S.target); }));
$("presets").innerHTML = PRESETS.map(p => `<button class="chip">${p.name}</button>`).join("");
$("presets").querySelectorAll(".chip").forEach((c, k) => c.onclick = () => start({ ...PRESETS[k] }));
async function go() {
  const txt = $("target").value.trim(); if (!txt) return;
  const preset = PRESETS.find(p => p.name.toLowerCase() === txt.toLowerCase()); if (preset) return start({ ...preset });
  const c = parseCoords(txt); if (c) return start({ name: null, typed: true, ra: c.ra, dec: c.dec, pmra: 0, pmdec: 0 });
  if (OFFLINE) return setStatus("Name lookup needs the server", "In preview mode, type coordinates instead, for example 08h55m10.8s -07d14m42.5s.");
  try {
    setStatus("Looking up " + txt + "…", "Asking the CDS name resolver for its position.");
    const r = await api("/api/resolve", { name: txt });
    start({ name: txt, ra: r.ra, dec: r.dec, pmra: r.pmra, pmdec: r.pmdec, epoch: r.epoch });
  } catch (e) { setStatus("Couldn't find that", e.message === "Failed to fetch" ? "The page can't reach its server. Start it with: python server.py" : e.message); }
}
$("goBtn").onclick = go;
$("target").addEventListener("keydown", e => { if (e.key === "Enter") go(); });


// ================= preview mode (no server: simulated data + your own files) =================
let OFFLINE = false, DEMO_BACKEND = false;
async function checkServer() {          // true = no backend: fall back to the simulated preview
  if (location.protocol === "file:") return true;
  try {
    const r = await fetch("/api/ping");
    if (!(r.ok && (r.headers.get("content-type") || "").includes("json"))) return true;
    DEMO_BACKEND = !!(await r.json()).demo;
    return false;
  } catch { return true; }
}
function seeded(seed) { let s = seed >>> 0 || 1; return () => (s = (s * 1664525 + 1013904223) >>> 0) / 4294967296; }
function gauss(rnd) { return Math.sqrt(-2 * Math.log(rnd() + 1e-12)) * Math.cos(2 * Math.PI * rnd()); }
// Rough brightness of each test target in SPHEREx's six bands (made-up units), plus how it changes with time.
const SIM = {
  "Barnard's Star": { sed: [60, 70, 62, 42, 30, 26], vary: () => 1 },
  "Epsilon Indi B": { sed: [6, 14, 10, 3, 11, 16], vary: () => 1 },
  "WISE 0855−0714": { sed: [0, 0, 0.2, 0.8, 3.5, 6], vary: () => 1 },
  "RW Aurigae": { sed: [28, 34, 40, 46, 42, 38], vary: (mjd, v) => [1, 0.5, 0.85, 0.7][v % 4] },
  "R Draconis": { sed: [70, 115, 105, 68, 48, 40], vary: mjd => 1 + 0.6 * Math.sin(2 * Math.PI * (mjd - 60800) / 245) }
};
function startDemo(target) {
  const G = S.grid, N = G.N, run = S.run, sim = SIM[target.name] || { sed: [20, 22, 20, 16, 12, 10], vary: () => 1 };
  const rnd = seeded(Math.round(target.ra * 1000 + target.dec * 7919));
  const stars = []; for (let k = 0; k < Math.round(70 * (G.fov / 0.1) ** 2); k++) stars.push([rnd() * N, rnd() * N, Math.pow(rnd(), 2.5) * 45 + 0.8, 0.5 + rnd()]);
  const visits = [60820, 61005, 61190];
  S.rows = []; S.frames = [];
  visits.forEach((mjd0, v) => { for (let b = 1; b <= 6; b++) for (let k = 0; k < 3; k++) {
    const t = mjd0 + b * 1.4 + k * 0.35;
    S.rows.push({ obs_id: `simulated-${v + 1}-D${b}-${k + 1}`, t_min: t, t_max: t + 0.002, em_min: BANDS[b - 1].lo * 1e-6, em_max: BANDS[b - 1].hi * 1e-6, bandN: b, access_url: "", v });
  } });
  for (const row of S.rows) {
    const b = row.bandN, year = mjdToYear(row.t_min), noise = 0.35, data = new Float32Array(N * N);
    for (let i = 0; i < data.length; i++) data[i] = gauss(rnd) * noise;
    const add = (x, y, f) => { const r = 4; for (let j = Math.max(0, Math.floor(y - r)); j <= Math.min(N - 1, y + r); j++) for (let i = Math.max(0, Math.floor(x - r)); i <= Math.min(N - 1, x + r); i++)
      data[j * N + i] += f * Math.exp(-((i - x) ** 2 + (j - y) ** 2) / (2 * 0.75 ** 2)); };
    for (const [x, y, f, slope] of stars) add(x, y, f * Math.pow(0.85 + 0.05 * b, slope * 2));
    const tp = targetAt(year), [tx, ty] = G.toGrid(tp.ra, tp.dec);
    add(tx, ty, sim.sed[b - 1] * sim.vary(row.t_min, row.v));
    if (row.v === 1 && b === 3) add(N * 0.25 + (row.t_min % 1) * 6, N * 0.7 - (row.t_min % 1) * 3, 14);   // a passing asteroid
    S.frames.push(frameFromData(row, data));
  }
  S.sim = true; $("simTag").hidden = false; $("simTag").textContent = "Simulated preview";
  renderBands(); renderTargetCard(); rebuild();
  setStatus(null);
}
function frameFromData(row, data) {
  const { med, sig } = robust(data), norm = new Float32Array(data.length);
  for (let i = 0; i < data.length; i++) norm[i] = (data[i] - med) / sig;
  const mjd = row.t_max ? (row.t_min + row.t_max) / 2 : row.t_min, year = mjdToYear(mjd);
  const tp = S.target ? targetAt(year) : { ra: S.grid.ra0, dec: S.grid.dec0 }, [tx, ty] = S.grid.toGrid(tp.ra, tp.dec);
  return { row, mjd, year, band: row.bandN, data, norm, flux: photometry(data, S.grid, tx, ty), tx, ty };
}
async function openFiles(fileList) {
  const files = [...fileList]; if (!files.length) return;
  stop(); S.run++;
  const parsed = [];
  for (const f of files) {
    try {
      const buf = await f.arrayBuffer(), hdus = parseFITS(buf);
      const imgs = hdus.filter(x => x.hdr.NAXIS === 2 && x.hdr.XTENSION !== "BINTABLE" && x.hdr.NAXIS1 > 1);
      const ih = imgs.find(x => x.hdr.EXTNAME === "IMAGE") || imgs[0]; if (!ih) continue;
      const wcs = makeWCS(ih.hdr); if (!wcs) continue;
      const img = readImage(buf, ih);
      const zh = imgs.find(x => x.hdr.EXTNAME === "ZODI" && x.hdr.NAXIS1 === img.w && x.hdr.NAXIS2 === img.h);
      if (zh) { const z = readImage(buf, zh).data; for (let i = 0; i < img.data.length; i++) img.data[i] -= z[i]; }
      const det = +ih.hdr.DETECTOR || +hdus[0].hdr.DETECTOR || 0;
      parsed.push({ name: f.name, img, wcs, mjd: headerTime(ih.hdr) || headerTime(hdus[0].hdr), band: det >= 1 && det <= 6 ? det : 0 });
    } catch (e) { /* skip unreadable files */ }
  }
  if (!parsed.length) { setStatus("Couldn't read those files", "Use SPHEREx spectral image or cutout files (.fits) from IRSA."); return; }
  const first = parsed[0], [cra, cdec] = first.wcs.pix2sky(first.img.w / 2 + 0.5, first.img.h / 2 + 0.5);
  const extent = Math.min(first.img.w, first.img.h) * PIX / 3600 / 1.45;
  S.target = { name: "Your files", typed: true, ra: cra / d2r, dec: cdec / d2r, pmra: 0, pmdec: 0, note: `${parsed.length} file${parsed.length > 1 ? "s" : ""} opened from your computer.` };
  S.grid = makeGrid(S.target.ra, S.target.dec, Math.min(+$("fov").value, Math.max(0.02, extent)));
  S.rows = parsed.map((p, k) => ({ obs_id: p.name, t_min: p.mjd ?? 60800 + k * 30, t_max: null, em_min: p.band ? BANDS[p.band - 1].lo * 1e-6 : NaN, em_max: p.band ? BANDS[p.band - 1].hi * 1e-6 : NaN, bandN: p.band, access_url: "" }));
  S.frames = parsed.map((p, k) => frameFromData(S.rows[k], reproject(p.img, p.wcs, S.grid)));
  S.gaia = []; S.idx = 0; S.band = 0; S.sim = false; $("simTag").hidden = true;   // your own files are real data
  document.querySelectorAll("#presets .chip").forEach(c => c.classList.remove("on"));
  $("target").value = `${fmtRA(S.target.ra)} ${fmtDec(S.target.dec)}`;
  renderBands(); renderTargetCard(); rebuild(); setStatus(null);
  if (parsed.some(p => !p.mjd)) { setStatus("Some files have no date", "Those frames are shown in the order you picked them. Click anywhere on this message to close it."); $("status").onclick = () => setStatus(null); }
}
$("fileIn").addEventListener("change", e => openFiles(e.target.files));
for (const el of [$("viewer"), document.querySelector(".filebtn")]) {
  el.addEventListener("dragover", e => { e.preventDefault(); document.querySelector(".filebtn").classList.add("over"); });
  el.addEventListener("dragleave", () => document.querySelector(".filebtn").classList.remove("over"));
  el.addEventListener("drop", e => { e.preventDefault(); document.querySelector(".filebtn").classList.remove("over"); openFiles(e.dataTransfer.files); });
}

// ================= URL & sharing =================
function shareURL() {
  const p = new URLSearchParams(), t = S.target;
  if (t) { if (t.name && !t.typed) p.set("target", t.name); else { p.set("ra", t.ra.toFixed(6)); p.set("dec", t.dec.toFixed(6)); } }
  p.set("fov", $("fov").value); if (S.band) p.set("band", S.band); if ($("release").value !== "spherex_dr1") p.set("release", $("release").value);
  if (S.framesMode !== "visit") p.set("frames", S.framesMode); if (S.show !== "image") p.set("show", S.show);
  return location.origin + location.pathname + "?" + p;
}
function updateURL() { if (!OFFLINE && location.protocol.startsWith("http")) try { history.replaceState(null, "", shareURL()); } catch { } }
async function copy(text, label) {
  try { await navigator.clipboard.writeText(text); $("shareMsg").textContent = label + " copied"; }
  catch { $("shareMsg").textContent = "Copy didn't work here. Select the text below."; $("embedCode").hidden = false; $("embedCode").textContent = text; }
}
$("copyLink").onclick = () => OFFLINE ? ($("shareMsg").textContent = "Share links work when the tool runs from server.py or is hosted online.") : copy(shareURL(), "Link");
$("copyEmbed").onclick = () => { const code = `<iframe src="${shareURL()}" width="1100" height="820" style="border:0" title="SPHEREx Time-Lapse"></iframe>`; $("embedCode").hidden = false; $("embedCode").textContent = code; copy(code, "Embed code"); };
$("savePng").onclick = () => { if (!S.items.length) return; const a = document.createElement("a"); a.download = `spherex_${(S.target.name || "position").replace(/\W+/g, "_")}_${S.idx + 1}.png`; a.href = view.toDataURL("image/png"); a.click(); };

// ================= boot =================
renderBands(); renderAll();
(async function boot() {
  OFFLINE = await checkServer();
  if (OFFLINE) {
    $("ovGaia").checked = false; $("ovGaia").disabled = true; $("gaiaText").textContent = "Gaia DR3 stars (needs server.py)";
    $("savePng").hidden = true; $("release").disabled = true; $("maxn").disabled = true;
    setStatus("Preview mode", "This page can't reach NASA's archive from here, so the five test targets play simulated images that show how the tool works. Run server.py for real SPHEREx data, or open your own SPHEREx FITS files with the button on the left.");
    return;
  }
  if (DEMO_BACKEND) setStatus("Demo mode", "The backend is running with --demo, so every image here is synthetic. Pick a test target to see how the tool works.");
  const q = new URLSearchParams(location.search);
  if (q.get("fov")) $("fov").value = q.get("fov");
  if (q.get("release")) $("release").value = q.get("release");
  if (q.get("band")) S.band = +q.get("band") || 0;
  for (const [k, id] of [["frames", "framesMode"], ["show", "showMode"]]) if (q.get(k)) { const b = $(id).querySelector(`[data-v="${q.get(k)}"]`); if (b) b.click(); }
  renderBands();
  if (q.get("ra") && q.get("dec")) { const ra = +q.get("ra"), dec = +q.get("dec"); if (isFinite(ra) && isFinite(dec)) return start({ name: null, typed: true, ra, dec, pmra: 0, pmdec: 0 }); }
  if (q.get("target")) { $("target").value = q.get("target"); go(); }
})();
