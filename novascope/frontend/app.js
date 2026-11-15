"use strict";
// ---------- data ----------
const D = window.STARLIGHT_DATA;
// "inferno" colour scale (matplotlib, BSD-style licence): dark purple → magenta → orange → yellow, as in IRSA's M101 figures
const INFERNO = (h => { const a = new Uint8Array(768); for (let i = 0; i < 768; i++) a[i] = parseInt(h.substr(i * 2, 2), 16); return a; })("00000401000501010601010802010a02020c02020e03021004031204031405041706041907051b08051d09061f0a07220b07240c08260d08290e092b10092d110a30120a32140b34150b37160b39180c3c190c3e1b0c411c0c431e0c451f0c48210c4a230c4c240c4f260c51280b53290b552b0b572d0b592f0a5b310a5c320a5e340a5f3609613809623909633b09643d09653e0966400a67420a68440a68450a69470b6a490b6a4a0c6b4c0c6b4d0d6c4f0d6c510e6c520e6d540f6d550f6d57106e59106e5a116e5c126e5d126e5f136e61136e62146e64156e65156e67166e69166e6a176e6c186e6d186e6f196e71196e721a6e741a6e751b6e771c6d781c6d7a1d6d7c1d6d7d1e6d7f1e6c801f6c82206c84206b85216b87216b88226a8a226a8c23698d23698f24699025689225689326679526679727669827669a28659b29649d29649f2a63a02a63a22b62a32c61a52c60a62d60a82e5fa92e5eab2f5ead305dae305cb0315bb1325ab3325ab43359b63458b73557b93556ba3655bc3754bd3853bf3952c03a51c13a50c33b4fc43c4ec63d4dc73e4cc83f4bca404acb4149cc4248ce4347cf4446d04545d24644d34743d44842d54a41d74b3fd84c3ed94d3dda4e3cdb503bdd513ade5238df5337e05536e15635e25734e35933e45a31e55c30e65d2fe75e2ee8602de9612bea632aeb6429eb6628ec6726ed6925ee6a24ef6c23ef6e21f06f20f1711ff1731df2741cf3761bf37819f47918f57b17f57d15f67e14f68013f78212f78410f8850ff8870ef8890cf98b0bf98c0af98e09fa9008fa9207fa9407fb9606fb9706fb9906fb9b06fb9d07fc9f07fca108fca309fca50afca60cfca80dfcaa0ffcac11fcae12fcb014fcb216fcb418fbb61afbb81dfbba1ffbbc21fbbe23fac026fac228fac42afac62df9c72ff9c932f9cb35f8cd37f8cf3af7d13df7d340f6d543f6d746f5d949f5db4cf4dd4ff4df53f4e156f3e35af3e55df2e661f2e865f2ea69f1ec6df1ed71f1ef75f1f179f2f27df2f482f3f586f3f68af4f88ef5f992f6fa96f8fb9af9fc9dfafda1fcffa4");
function inferno(t, px, o) { const k = Math.max(0, Math.min(255, Math.round(t * 255))) * 3; px[o] = INFERNO[k]; px[o + 1] = INFERNO[k + 1]; px[o + 2] = INFERNO[k + 2]; }
function bytes(b64) { const s = atob(b64), u = new Uint8Array(s.length); for (let i = 0; i < s.length; i++) u[i] = s.charCodeAt(i); return u; }
const S = new Float32Array(bytes(D.stars).buffer);       // x y z vx vy vz absmag ci  (pc, pc/yr, J2000 equatorial)
const CON = bytes(D.con);
const HIP = new Int32Array(bytes(D.hip).buffer);
const HD = new Int32Array(bytes(D.hd).buffer);
const MWQ = new Int16Array(bytes(D.mw).buffer), MWL = bytes(D.mwl);
const N = D.n, PC2LY = 3.26156, PCYR2KMS = 977792;
const NOW_YEAR = new Date().getFullYear();
const SPAN = 200000;
const GREEK = {Alp:"α",Bet:"β",Gam:"γ",Del:"δ",Eps:"ε",Zet:"ζ",Eta:"η",The:"θ",Iot:"ι",Kap:"κ",Lam:"λ",Mu:"μ",Nu:"ν",Xi:"ξ",Omi:"ο",Pi:"π",Rho:"ρ",Sig:"σ",Tau:"τ",Ups:"υ",Phi:"φ",Chi:"χ",Psi:"ψ",Ome:"ω"};
const proper = new Array(N).fill(null), bayer = new Array(N).fill(null);
D.names.forEach(([i, n]) => proper[i] = n);
D.bayer.forEach(([i, b]) => { const m = /^([A-Za-z]+)(-?\d*)$/.exec(b); bayer[i] = m && GREEK[m[1]] ? GREEK[m[1]] + (m[2] ? m[2].replace("-", "") : "") : b; });
const conAbbr = i => CON[i] < 255 ? D.cons[CON[i]] : "";
const conName = a => D.cname[a] || a;
function starName(i) {
  if (proper[i]) return proper[i];
  if (bayer[i] && conAbbr(i)) return bayer[i] + " " + conAbbr(i);
  if (HIP[i]) return "HIP " + HIP[i];
  return "HD " + HD[i];
}

// ---------- epoch math ----------
const P = new Float32Array(N * 3), MAG = new Float32Array(N), DIST = new Float32Array(N);
function posAt(i, dt, out) {
  const o = i * 8;
  out[0] = S[o] + S[o + 3] * dt; out[1] = S[o + 1] + S[o + 4] * dt; out[2] = S[o + 2] + S[o + 5] * dt;
  const d = Math.hypot(out[0], out[1], out[2]);
  out[3] = d; out[4] = S[o + 6] + 5 * Math.log10(d) - 5;
  return out;
}
let epochDt = NaN;
function computeEpoch(year) {
  const dt = year - 2000;
  if (dt === epochDt) return;
  epochDt = dt;
  const t = [0, 0, 0, 0, 0];
  for (let i = 0; i < N; i++) {
    posAt(i, dt, t);
    P[i * 3] = t[0] / t[3]; P[i * 3 + 1] = t[1] / t[3]; P[i * 3 + 2] = t[2] / t[3];
    DIST[i] = t[3]; MAG[i] = t[4];
  }
}
// Precession: north celestial pole circles the ecliptic pole once per 25,772 yr (obliquity held at 23.44°).
const EPS = 23.44 * Math.PI / 180, PREC = 25772;
function poleAt(year, south) {
  const lon = (90 - 360 * (year - 2000) / PREC) * Math.PI / 180, lat = (90 - 23.44) * Math.PI / 180;
  const xe = Math.cos(lat) * Math.cos(lon), ye = Math.cos(lat) * Math.sin(lon), ze = Math.sin(lat);
  const v = [xe, ye * Math.cos(EPS) - ze * Math.sin(EPS), ye * Math.sin(EPS) + ze * Math.cos(EPS)];
  return south ? v.map(c => -c) : v;
}
function poleStar(year) {
  const p = poleAt(year, false);
  let best = -1, bestAng = 99;
  for (let i = 0; i < N; i++) {
    if (MAG[i] > 4.0) continue;
    const c = P[i * 3] * p[0] + P[i * 3 + 1] * p[1] + P[i * 3 + 2] * p[2];
    const a = Math.acos(Math.min(1, c)) * 180 / Math.PI;
    if (a < bestAng) { bestAng = a; best = i; }
  }
  return bestAng < 5 ? { i: best, ang: bestAng } : null;
}

// ---------- state ----------
const state = {
  year: NOW_YEAR, ra: 200 * Math.PI / 180, dec: 36 * Math.PI / 180, fov: 95, sel: -1, tab: "pane-star",
  layers: { lines: true, conNames: true, names: true, milky: true, trails: false, poles: true, grid: false, spherex: true, atlas: true }
};
let playing = 0, dirty = true, journey = null, cutoutMark = null, cutRect = null;

// ---------- canvas & projection ----------
const cv = document.getElementById("sky"), ctx = cv.getContext("2d");
let W = 0, H = 0, DPR = 1;
function resize() {
  DPR = Math.min(2, window.devicePixelRatio || 1);
  W = window.innerWidth; H = window.innerHeight;
  cv.width = Math.round(W * DPR); cv.height = Math.round(H * DPR);
  dirty = true;
}
window.addEventListener("resize", resize);
let F = [1, 0, 0], U = [0, 0, 1], R = [0, -1, 0], scale = 1, cx = 0, cy = 0;
function setBasis() {
  const cr = Math.cos(state.ra), sr = Math.sin(state.ra), cd = Math.cos(state.dec), sd = Math.sin(state.dec);
  F = [cd * cr, cd * sr, sd];
  U = [-sd * cr, -sd * sr, cd];
  R = [F[1] * U[2] - F[2] * U[1], F[2] * U[0] - F[0] * U[2], F[0] * U[1] - F[1] * U[0]];
  scale = (Math.max(W, H) / 2) / (2 * Math.tan(state.fov * Math.PI / 720));
  cx = W / 2; cy = H / 2;
}
const PR = [0, 0]; let FRONT = false;   // FRONT: point is in front of the viewer (PR is valid even if off-screen)
function project(x, y, z) {
  const c = x * F[0] + y * F[1] + z * F[2];
  FRONT = c >= -0.9;
  if (!FRONT) return false;
  const k = 2 / (1 + c) * scale;
  PR[0] = cx + (x * R[0] + y * R[1] + z * R[2]) * k;
  PR[1] = cy - (x * U[0] + y * U[1] + z * U[2]) * k;
  return PR[0] > -60 && PR[0] < W + 60 && PR[1] > -60 && PR[1] < H + 60;
}
const TUTORIAL_RA = 304.693508808, TUTORIAL_DEC = 42.4436872991, M101_RA = 210.80227, M101_DEC = 54.34895, CLOUD_RA = 60.8109955, CLOUD_DEC = 51.3404038;
const radec = (ra, dec) => [Math.cos(dec) * Math.cos(ra), Math.cos(dec) * Math.sin(ra), Math.sin(dec)];

// ---------- star colours & sprites ----------
function bvToRGB(bv) {
  bv = Math.max(-0.35, Math.min(2.0, bv));
  const T = 4600 * (1 / (0.92 * bv + 1.7) + 1 / (0.92 * bv + 0.62)), t = T / 100;
  let r, g, b;
  r = t <= 66 ? 255 : 329.7 * Math.pow(t - 60, -0.1332);
  g = t <= 66 ? 99.47 * Math.log(t) - 161.12 : 288.12 * Math.pow(t - 60, -0.0755);
  b = t >= 66 ? 255 : t <= 19 ? 0 : 138.52 * Math.log(t - 10) - 305.04;
  const cl = v => Math.max(0, Math.min(255, v));
  const mix = 0.3;   // soften toward white, as the eye sees them
  return [cl(r) * (1 - mix) + 255 * mix, cl(g) * (1 - mix) + 255 * mix, cl(b) * (1 - mix) + 255 * mix].map(Math.round);
}
const NB = 16, bucketOf = new Uint8Array(N), bucketRGB = [], sprites = [];
for (let k = 0; k < NB; k++) bucketRGB.push(bvToRGB(-0.3 + (k + 0.5) * 2.3 / NB));
for (let i = 0; i < N; i++) bucketOf[i] = Math.max(0, Math.min(NB - 1, Math.floor((S[i * 8 + 7] + 0.3) / 2.3 * NB)));
const order = Array.from({ length: N }, (_, i) => i).sort((a, b) => bucketOf[a] - bucketOf[b]);
for (const [r, g, b] of bucketRGB) {
  const c = document.createElement("canvas"); c.width = c.height = 64;
  const x = c.getContext("2d"), gr = x.createRadialGradient(32, 32, 0, 32, 32, 32);
  gr.addColorStop(0, "rgba(255,255,255,1)"); gr.addColorStop(0.16, `rgba(${r},${g},${b},0.95)`);
  gr.addColorStop(0.38, `rgba(${r},${g},${b},0.28)`); gr.addColorStop(1, `rgba(${r},${g},${b},0)`);
  x.fillStyle = gr; x.fillRect(0, 0, 64, 64); sprites.push(c);
}
const mwSprite = (() => {
  const c = document.createElement("canvas"); c.width = c.height = 32;
  const x = c.getContext("2d"), gr = x.createRadialGradient(16, 16, 0, 16, 16, 16);
  gr.addColorStop(0, "rgba(170,185,230,1)"); gr.addColorStop(1, "rgba(170,185,230,0)");
  x.fillStyle = gr; x.fillRect(0, 0, 32, 32); return c;
})();
const MWV = new Float32Array(MWL.length * 3);
for (let i = 0; i < MWL.length; i++) {
  const v = radec((MWQ[i * 2] / 100) * Math.PI / 180, (MWQ[i * 2 + 1] / 100) * Math.PI / 180);
  MWV.set(v, i * 3);
}
const BAND_COLORS = ["#7fa8ff", "#69d0d8", "#93d27c", "#e3c95c", "#ee9a55", "#e8697c"];
let ATLAS = { n: 0, v: [], mjd: [], det: [], days: [] };
const SPOTS = [{ v: radec(TUTORIAL_RA * Math.PI / 180, TUTORIAL_DEC * Math.PI / 180), label: "One SPHEREx picture (example)", sec: "secIntro", fov: 10 },
  { v: radec(M101_RA * Math.PI / 180, M101_DEC * Math.PI / 180), label: "Pinwheel Galaxy (M101)", m101: true, sec: "secCutouts" },
  { v: radec(CLOUD_RA * Math.PI / 180, CLOUD_DEC * Math.PI / 180), label: "Dust cloud (example)", sec: "secCube", fov: 12 }];
let spotHits = [];      // where each marker (and its label) is on screen, for clicks
// constellation line stars & label anchors
const conStars = {};
for (const [abbr, segs] of Object.entries(D.lines)) conStars[abbr] = [...new Set(segs)];

// ---------- drawing ----------
const SX = new Float32Array(N), SY = new Float32Array(N), SR = new Float32Array(N), VIS = new Uint8Array(N);
function magLimit() { return Math.min(7.0, 6.2 + 1.3 * Math.log10(Math.max(1, 90 / state.fov))); }

function draw() {
  computeEpoch(state.year);
  setBasis();
  ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
  const bg = ctx.createRadialGradient(W * 0.5, H * 0.45, 0, W * 0.5, H * 0.5, Math.max(W, H) * 0.75);
  bg.addColorStop(0, "#0d1226"); bg.addColorStop(1, "#04050c");
  ctx.globalCompositeOperation = "source-over"; ctx.fillStyle = bg; ctx.fillRect(0, 0, W, H);
  const L = state.layers;

  // Milky Way
  if (L.milky) {
    ctx.globalCompositeOperation = "lighter";
    const base = (Math.PI / 180) * 1.02 * scale;
    for (let i = 0; i < MWL.length; i++) {
      const x = MWV[i * 3], y = MWV[i * 3 + 1], z = MWV[i * 3 + 2];
      if (!project(x, y, z)) continue;
      const c = x * F[0] + y * F[1] + z * F[2], sz = base * 2 / (1 + c) * 2.6;
      ctx.globalAlpha = [0, 0.045, 0.07, 0.09, 0.11, 0.13][MWL[i]];
      ctx.drawImage(mwSprite, PR[0] - sz / 2, PR[1] - sz / 2, sz, sz);
    }
    ctx.globalAlpha = 1;
  }

  // grid
  if (L.grid) {
    ctx.globalCompositeOperation = "source-over"; ctx.strokeStyle = "rgba(95,179,212,0.16)"; ctx.lineWidth = 1;
    const line = pts => { ctx.beginPath(); let pen = false;
      for (const v of pts) { if ((project(v[0], v[1], v[2]), FRONT)) { pen ? ctx.lineTo(PR[0], PR[1]) : ctx.moveTo(PR[0], PR[1]); pen = true; } else pen = false; }
      ctx.stroke(); };
    for (let h = 0; h < 24; h += 2) { const pts = []; for (let d = -88; d <= 88; d += 2) pts.push(radec(h * Math.PI / 12, d * Math.PI / 180)); line(pts); }
    for (let d = -75; d <= 75; d += 15) { const pts = []; for (let a = 0; a <= 360; a += 2) pts.push(radec(a * Math.PI / 180, d * Math.PI / 180)); line(pts); }
  }

  // constellation lines
  ctx.globalCompositeOperation = "source-over";
  if (L.lines) {
    ctx.beginPath();
    for (const segs of Object.values(D.lines)) {
      for (let k = 0; k < segs.length; k += 2) {
        const a = segs[k], b = segs[k + 1];
        const va = project(P[a * 3], P[a * 3 + 1], P[a * 3 + 2]), fa = FRONT; const ax = PR[0], ay = PR[1];
        const vb = project(P[b * 3], P[b * 3 + 1], P[b * 3 + 2]), fb = FRONT; const bx = PR[0], by = PR[1];
        if (!fa || !fb || (!va && !vb)) continue;
        const dx = bx - ax, dy = by - ay, len = Math.hypot(dx, dy);
        if (len < 8 || len > Math.max(W, H) * 0.7) continue;
        const g = Math.min(5, len * 0.18) / len;
        ctx.moveTo(ax + dx * g, ay + dy * g); ctx.lineTo(bx - dx * g, by - dy * g);
      }
    }
    // electric-blue glow: soft layers, widest and faintest first, then a thin bright core
    ctx.save(); ctx.lineCap = "round"; ctx.globalCompositeOperation = "lighter";
    for (const [w, c] of [[12, "rgba(30,110,255,0.06)"], [7, "rgba(50,145,255,0.11)"], [3.6, "rgba(100,190,255,0.22)"], [1.7, "rgba(185,228,255,0.5)"], [0.75, "rgba(240,250,255,0.85)"]]) {
      ctx.lineWidth = w; ctx.strokeStyle = c; ctx.stroke();
    }
    ctx.restore();
  }
  if (L.conNames && state.fov < 150) {
    ctx.font = "italic 600 14px 'Cormorant Garamond', Georgia, serif"; ctx.textAlign = "center"; ctx.fillStyle = "rgba(124,200,255,0.62)";
    for (const [abbr, ids] of Object.entries(conStars)) {
      if (!ids.length) continue;
      let x = 0, y = 0, z = 0; for (const i of ids) { x += P[i * 3]; y += P[i * 3 + 1]; z += P[i * 3 + 2]; }
      const n = Math.hypot(x, y, z); if (!project(x / n, y / n, z / n)) continue;
      ctx.fillText(conName(abbr).toUpperCase().split("").join(" "), PR[0], PR[1]);
    }
    ctx.textAlign = "left";
  }

  // poles
  const yr = state.year;
  if (L.poles) {
    for (const south of [false, true]) {
      ctx.strokeStyle = "rgba(95,179,212,0.35)"; ctx.setLineDash([2, 5]); ctx.beginPath(); let pen = false;
      for (let k = 0; k <= 96; k++) { const p = poleAt(2000 + k * PREC / 96, south);
        if ((project(p[0], p[1], p[2]), FRONT)) { pen ? ctx.lineTo(PR[0], PR[1]) : ctx.moveTo(PR[0], PR[1]); pen = true; } else pen = false; }
      ctx.stroke(); ctx.setLineDash([]);
      const p = poleAt(yr, south);
      if (project(p[0], p[1], p[2])) {
        ctx.strokeStyle = "rgba(95,179,212,0.9)"; ctx.lineWidth = 1.2; ctx.beginPath();
        ctx.arc(PR[0], PR[1], 7, 0, 7); ctx.moveTo(PR[0] - 12, PR[1]); ctx.lineTo(PR[0] - 4, PR[1]); ctx.moveTo(PR[0] + 4, PR[1]); ctx.lineTo(PR[0] + 12, PR[1]);
        ctx.moveTo(PR[0], PR[1] - 12); ctx.lineTo(PR[0], PR[1] - 4); ctx.moveTo(PR[0], PR[1] + 4); ctx.lineTo(PR[0], PR[1] + 12); ctx.stroke();
        ctx.font = "11.5px 'IBM Plex Sans', sans-serif"; ctx.fillStyle = "rgba(95,179,212,0.95)";
        ctx.fillText(south ? "South celestial pole" : "North celestial pole", PR[0] + 14, PR[1] + 16);
      }
    }
  }

  // trails: where bright stars head over the next 30,000 years
  const tmp = [0, 0, 0, 0, 0], dt = yr - 2000;
  if (L.trails) {
    ctx.lineWidth = 1.2;
    for (let i = 0; i < N; i++) {
      if (MAG[i] > 3.6) continue;
      if (!project(P[i * 3], P[i * 3 + 1], P[i * 3 + 2])) continue;
      const x0 = PR[0], y0 = PR[1];
      posAt(i, dt + 30000, tmp); if (!(project(tmp[0] / tmp[3], tmp[1] / tmp[3], tmp[2] / tmp[3]), FRONT)) continue;
      const g = ctx.createLinearGradient(x0, y0, PR[0], PR[1]);
      g.addColorStop(0, "rgba(95,179,212,0.75)"); g.addColorStop(1, "rgba(95,179,212,0)");
      ctx.strokeStyle = g; ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(PR[0], PR[1]); ctx.stroke();
    }
  }
  // selected star's full path
  if (state.sel >= 0) {
    const i = state.sel; ctx.lineWidth = 1.4;
    let prev = null;
    for (let k = -40; k <= 40; k++) {
      const d = NOW_YEAR - 2000 + k * SPAN / 40; posAt(i, d, tmp);
      const ok = (project(tmp[0] / tmp[3], tmp[1] / tmp[3], tmp[2] / tmp[3]), FRONT); const cur = ok ? [PR[0], PR[1]] : null;
      if (prev && cur && Math.hypot(cur[0] - prev[0], cur[1] - prev[1]) < W * 0.5) {
        ctx.strokeStyle = k <= 0 ? "rgba(235,154,92,0.4)" : "rgba(95,179,212,0.4)";
        ctx.beginPath(); ctx.moveTo(prev[0], prev[1]); ctx.lineTo(cur[0], cur[1]); ctx.stroke();
      }
      prev = cur;
    }
  }

  // stars
  ctx.globalCompositeOperation = "lighter";
  const mlim = magLimit(), zoomBoost = Math.max(0.8, Math.min(1.8, 1 + 0.25 * Math.log2(90 / state.fov)));
  VIS.fill(0);
  let lastB = -1;
  for (let n = 0; n < N; n++) {
    const i = order[n], m = MAG[i];
    if (m > mlim) continue;
    if (!project(P[i * 3], P[i * 3 + 1], P[i * 3 + 2])) continue;
    const x = PR[0], y = PR[1], mm = mlim - m, r = (0.3 + 0.4 * mm) * zoomBoost;
    SX[i] = x; SY[i] = y; SR[i] = r; VIS[i] = 1;
    const b = bucketOf[i];
    if (r < 1.0) {
      if (b !== lastB) { const c = bucketRGB[b]; ctx.fillStyle = `rgb(${c[0]},${c[1]},${c[2]})`; lastB = b; }
      ctx.globalAlpha = Math.min(1, 0.25 + mm * 0.6); const s = r * 2 + 0.5;
      ctx.fillRect(x - s / 2, y - s / 2, s, s);
    } else {
      ctx.globalAlpha = 1; const s = r * 4.2;
      ctx.drawImage(sprites[b], x - s / 2, y - s / 2, s, s);
      if (m < 0.6) { ctx.globalAlpha = 0.22; const g = r * 11; ctx.drawImage(sprites[b], x - g / 2, y - g / 2, g, g); }
    }
  }
  ctx.globalAlpha = 1; ctx.globalCompositeOperation = "source-over";

  // star names
  if (L.names) {
    const nameLim = state.fov > 120 ? 1.2 : state.fov > 60 ? 2.0 : state.fov > 25 ? 3.2 : 5;
    ctx.font = "12px 'IBM Plex Sans', sans-serif"; ctx.fillStyle = "rgba(235,228,211,0.8)";
    for (let i = 0; i < N; i++) {
      if (!VIS[i] || i === state.sel) continue;
      const lab = proper[i] && MAG[i] < nameLim ? proper[i] : (state.fov < 20 && bayer[i] && MAG[i] < 5.5 ? bayer[i] : null);
      if (lab) ctx.fillText(lab, SX[i] + SR[i] + 4, SY[i] - SR[i] - 1);
    }
  }
  // SPHEREx tutorial target (IRSA spherex-intro notebook)
  spotHits = [];
  if (L.spherex) for (const sp of SPOTS) {
    if (!project(sp.v[0], sp.v[1], sp.v[2])) continue;
    const x = PR[0], y = PR[1];
    ctx.font = "11.5px 'IBM Plex Sans', sans-serif"; spotHits.push({ sp, x0: x - 10, x1: x + 14 + ctx.measureText(sp.label).width, y0: y - 10, y1: y + 10 });
    ctx.strokeStyle = "rgba(124,200,255,0.95)"; ctx.lineWidth = 1.3; ctx.beginPath();
    ctx.moveTo(x, y - 7); ctx.lineTo(x + 7, y); ctx.lineTo(x, y + 7); ctx.lineTo(x - 7, y); ctx.closePath(); ctx.stroke();
    const dup = cutoutMark && sp.v[0] * cutoutMark[0] + sp.v[1] * cutoutMark[1] + sp.v[2] * cutoutMark[2] > 0.99998;      // same place as the cutout box: its label already names it
    if (!dup && !(sp.m101 && m101Anim)) { ctx.font = "11.5px 'IBM Plex Sans', sans-serif"; ctx.fillStyle = "rgba(124,200,255,0.95)"; ctx.fillText(sp.label, x + 11, y + 4); }
  }
  if (m101Anim) drawM101(m101Anim);
  // interstellar comet 3I/ATLAS: positions SPHEREx recorded with each image (only meaningful near the present day)
  if (L.atlas && Math.abs(state.year - 2025.8) < 300) drawAtlas();
  // SPHEREx cutout marker
  cutRect = null;
  if (cutoutMark && !(m101Anim && m101Anim.img > 0) && project(cutoutMark[0], cutoutMark[1], cutoutMark[2])) {
    const has = fitsState.aligned.length > 0, small = Math.max(has ? 28 : 8, ((fitsState.fov || 0.1) / 2 * Math.PI / 180) * scale), x = PR[0], y = PR[1];
    const g = has ? ease(cutGrow) : 0, half = small + (Math.max(small, cutBigHalf()) - small) * g;
    if (has) {      // the current frame of the SPHEREx video, playing right where it was taken
      ctx.save(); ctx.shadowColor = "rgba(40,130,255,0.8)"; ctx.shadowBlur = 20; ctx.fillStyle = "#000"; ctx.fillRect(x - half, y - half, half * 2, half * 2); ctx.restore();
      ctx.imageSmoothingEnabled = false; ctx.drawImage(fc, x - half, y - half, half * 2, half * 2); ctx.imageSmoothingEnabled = true;
      cutRect = { x, y, half };
      if (g > 0.6) {      // frame info and data tag inside the enlarged square
        const t1 = document.getElementById("fitsTag").textContent, sim = document.getElementById("simTag"), real = sim.classList.contains("real");
        ctx.globalAlpha = Math.min(1, (g - 0.6) / 0.3); ctx.font = "12px 'IBM Plex Mono', monospace";
        let w = ctx.measureText(t1).width; ctx.fillStyle = "rgba(0,0,0,0.65)"; ctx.fillRect(x - half + 6, y - half + 6, w + 12, 20); ctx.fillStyle = "#ebe4d3"; ctx.fillText(t1, x - half + 12, y - half + 20);
        if (!sim.hidden) { ctx.font = "11.5px 'IBM Plex Sans', sans-serif"; const t2 = sim.textContent; w = ctx.measureText(t2).width;
          ctx.fillStyle = "rgba(0,0,0,0.65)"; ctx.fillRect(x + half - w - 18, y + half - 26, w + 12, 20); ctx.fillStyle = real ? "#a9de96" : "#9097ab"; ctx.fillText(t2, x + half - w - 12, y + half - 12); }
        ctx.globalAlpha = 1;
      }
    }
    ctx.strokeStyle = "rgba(124,200,255,0.95)"; ctx.lineWidth = 1.2; ctx.strokeRect(x - half, y - half, half * 2, half * 2);
    ctx.font = g > 0.5 ? "italic 600 19px 'Cormorant Garamond', Georgia, serif" : "11.5px 'IBM Plex Sans', sans-serif"; ctx.fillStyle = "rgba(124,200,255,0.95)"; ctx.fillText(fitsState.label || "SPHEREx cutout", x - half, y - half - 8);
    if (has && !cutBig && g < 0.05) { ctx.font = "11px 'IBM Plex Sans', sans-serif"; ctx.fillStyle = "rgba(235,228,211,0.75)"; ctx.fillText("Click to enlarge", x - half, y + half + 14); }
  }
  // selection ring
  if (state.sel >= 0 && project(P[state.sel * 3], P[state.sel * 3 + 1], P[state.sel * 3 + 2])) {
    const r = Math.max(9, (SR[state.sel] || 2) * 2.2);
    ctx.strokeStyle = "rgba(124,200,255,0.95)"; ctx.lineWidth = 1.3; ctx.beginPath(); ctx.arc(PR[0], PR[1], r, 0, 7); ctx.stroke();
    ctx.font = "italic 600 16px 'Cormorant Garamond', Georgia, serif"; ctx.fillStyle = "#7cc8ff";
    ctx.fillText(starName(state.sel), PR[0] + r + 5, PR[1] + 5);
  }
}

function drawAtlas() {
  const A = ATLAS, an = atlasAnim && atlasAnim.started ? atlasAnim : null, upto = an ? an.i : A.n - 1; let prev = null;
  ctx.lineWidth = 1.6;
  for (let i = 0; i <= upto; i++) {
    const v = A.v[i], ok = (project(v[0], v[1], v[2]), FRONT), cur = ok ? [PR[0], PR[1]] : null;
    if (prev && cur) {
      const gap = A.mjd[i] - A.mjd[i - 1] > 30;
      ctx.strokeStyle = gap ? "rgba(130,205,230,0.45)" : "rgba(130,205,230,0.9)"; ctx.setLineDash(gap ? [3, 6] : []);
      ctx.beginPath(); ctx.moveTo(prev[0], prev[1]); ctx.lineTo(cur[0], cur[1]); ctx.stroke(); ctx.setLineDash([]);
    }
    prev = cur;
  }
  if (an && !an.done && prev && project(an.head[0], an.head[1], an.head[2])) {      // last bit, up to the comet
    const gap = A.mjd[upto + 1] - A.mjd[upto] > 30;
    ctx.strokeStyle = gap ? "rgba(130,205,230,0.45)" : "rgba(130,205,230,0.9)"; ctx.setLineDash(gap ? [3, 6] : []);
    ctx.beginPath(); ctx.moveTo(prev[0], prev[1]); ctx.lineTo(PR[0], PR[1]); ctx.stroke(); ctx.setLineDash([]);
  }
  for (let i = 0; i <= upto; i++) {
    const v = A.v[i]; if (!project(v[0], v[1], v[2])) continue;
    ctx.fillStyle = BAND_COLORS[A.det[i] - 1]; ctx.beginPath(); ctx.arc(PR[0], PR[1], 2.4, 0, 7); ctx.fill();
  }
  ctx.font = "11.5px 'IBM Plex Sans', sans-serif"; ctx.fillStyle = "rgba(130,205,230,0.95)";
  const many = state.fov <= 30;
  A.days.forEach((d, k) => {
    if (d.i > upto || (!many && !d.end)) return;
    const v = A.v[d.i]; if (!project(v[0], v[1], v[2])) return;
    ctx.fillText(d.label, PR[0] + 6, PR[1] + (k % 2 ? 15 : -7));
  });
  if (an && !an.done) { drawComet(an); return; }
  const last = A.v[A.n - 1];
  if (project(last[0], last[1], last[2])) { ctx.font = "italic 600 16px 'Cormorant Garamond', Georgia, serif"; ctx.fillText("3I/ATLAS", PR[0] - 70, PR[1] - 12); }
}
// ---------- year text & star panel ----------
const fmt = n => Math.round(n).toLocaleString("en-US");
function yearLabel(y) { y = Math.round(y); const n = y >= 1 ? y : 1 - y; return (n < 10000 ? String(n) : fmt(n)) + (y >= 1 ? " CE" : " BCE"); }
function updateHUD() {
  const y = Math.round(state.year), diff = y - NOW_YEAR;
  const yt = document.getElementById("yearText");
  yt.textContent = yearLabel(y); yt.className = "year" + (diff < 0 ? " past" : diff > 0 ? " future" : "");
  document.getElementById("whenText").textContent = diff === 0 ? "Today" : diff < 0 ? fmt(-diff) + " years ago" : "In " + fmt(diff) + " years";
  const ps = poleStar(state.year);
  document.getElementById("poleText").innerHTML = ps ? `North pole star: <b>${starName(ps.i)}</b>, ${ps.ang.toFixed(1)}° from the pole` : "No bright star near the north pole in this year";
  slider.value = Math.round(yearToSlider(state.year));
  renderStar();
  if (state.tab === "pane-star" && state.sel < 0) renderNames(false);
}
function colorWord(ci) { return ci < 0 ? "Blue-white" : ci < 0.3 ? "White" : ci < 0.6 ? "Yellow-white" : ci < 0.9 ? "Yellow" : ci < 1.4 ? "Orange" : "Red"; }

// ---------- fun facts (edit or add freely) ----------
// Key = the star's name as shown on the map, or "HIP <number>" for stars without a proper name.
// Each entry is a list of short sentences, shown first in the Star tab.
const STAR_FACTS = {
  "Sirius": ["The brightest star in the night sky.", "It has a tiny partner, Sirius B, a white dwarf about the size of Earth. Astronomers predicted it in 1844 from the way Sirius wobbles, and first saw it in 1862.", "In ancient Egypt, Sirius rising just before dawn signalled the yearly flooding of the Nile."],
  "Canopus": ["The second-brightest star in the night sky.", "Spacecraft have used Canopus as a guide star to keep track of which way they are pointing."],
  "Arcturus": ["The brightest star north of the celestial equator.", "Its light was used to switch on the lights of the 1933 Chicago World's Fair."],
  "Vega": ["The first star other than the Sun to be photographed, in 1850.", "About 14,000 years ago it was the north pole star, and it will be again around 13,700 CE.", "It spins so fast that its equator bulges outward."],
  "Rigil Kentaurus": ["Alpha Centauri A, the brightest star of the nearest star system to the Sun.", "Its partner Toliman circles it every 80 years, and the faint red dwarf Proxima Centauri orbits the pair from far out."],
  "Toliman": ["Alpha Centauri B, the smaller, more orange partner of Rigil Kentaurus. The two orbit each other every 80 years."],
  "Proxima Centauri": ["The closest star to the Sun, 4.24 light-years away.", "It has at least one planet, Proxima b, in the zone where liquid water could exist."],
  "Betelgeuse": ["A red supergiant so large that, in place of the Sun, it would reach past the orbit of Mars.", "In 2019 and 2020 it dimmed noticeably, an event called the Great Dimming, caused by a cloud of dust it puffed out.", "It is expected to explode as a supernova within about 100,000 years."],
  "Rigel": ["A blue supergiant tens of thousands of times more luminous than the Sun.", "It is usually the brightest star in Orion, marking the hunter's foot."],
  "Polaris": ["Today's North Star, less than 1° from the north celestial pole.", "It is a Cepheid, a star that pulses in brightness, and it has two smaller companions."],
  "Aldebaran": ["The red eye of Taurus the bull.", "It looks like part of the Hyades star cluster but is less than half as far away.", "The Pioneer 10 spacecraft is heading roughly in its direction, but would take about 2 million years to get there."],
  "Antares": ["A red supergiant at the heart of Scorpius.", "Its name means \"rival of Mars\", because its red colour looks like the planet's."],
  "Spica": ["Actually two hot blue stars, so close that they orbit each other every 4 days."],
  "Capella": ["It looks like one star but is four: two yellow giants orbiting each other, and a distant pair of red dwarfs."],
  "Procyon": ["Its name means \"before the dog\": it rises just before Sirius, the Dog Star.", "Like Sirius, it has a white dwarf companion."],
  "Altair": ["It spins about once every 9 hours, so fast that it is visibly flattened.", "In 2007 it became one of the first stars beyond the Sun to have its surface imaged."],
  "Deneb": ["One of the most luminous stars you can see by eye. Its distance is uncertain: estimates range from about 1,400 to 2,600 light-years.", "It marks the tail of Cygnus the swan and one corner of the Summer Triangle."],
  "Fomalhaut": ["It is surrounded by a wide ring of dust and ice, imaged by the Hubble and James Webb space telescopes."],
  "Pollux": ["The nearest giant star to the Sun.", "A giant planet has been found orbiting it."],
  "Castor": ["It looks like one star but is six: three pairs of stars orbiting together."],
  "Regulus": ["It spins about once every 16 hours, so fast that it is close to flying apart."],
  "Achernar": ["One of the flattest stars known: it spins so fast that its equator is much wider than its poles."],
  "Mira": ["Its name means \"wonderful\". Over about 11 months it brightens and fades, at times becoming too faint to see by eye.", "Ultraviolet images from NASA's GALEX showed a tail of gas about 13 light-years long trailing behind it."],
  "Algol": ["The \"Demon Star\": every 2.87 days it dims for a few hours when its fainter partner passes in front of it."],
  "Thuban": ["The north pole star when the pyramids of Egypt were built, around 2,600 BCE."],
  "Mizar": ["With Alcor beside it, it forms a pair once used as a test of good eyesight.", "Through a telescope Mizar itself splits in two. It was one of the first double stars found with a telescope."],
  "Alcor": ["Mizar's fainter partner, easy to see beside it on a dark night."],
  "Barnard's Star": ["The fastest-moving star across our sky: it crosses the width of the full Moon in about 175 years.", "The closest single star to the Sun, about 6 light-years away."],
  "Kapteyn's Star": ["An ancient star from the Milky Way's halo, possibly torn long ago from the star cluster Omega Centauri."],
  "Alnitak": ["The eastern star of Orion's Belt, close in the sky to the Horsehead and Flame nebulae."],
  "Alnilam": ["The middle star of Orion's Belt, a blue supergiant."],
  "Mintaka": ["The western star of Orion's Belt. It sits almost exactly on the celestial equator, so from anywhere on Earth it rises due east and sets due west."],
  "Bellatrix": ["Its name means \"female warrior\". It marks one of Orion's shoulders."],
  "Dubhe": ["With Merak, one of the two \"pointer\" stars of the Big Dipper: follow them to find Polaris."],
  "Merak": ["With Dubhe, one of the two \"pointer\" stars of the Big Dipper: follow them to find Polaris."],
  "Alkaid": ["The end of the Big Dipper's handle. Unlike the five middle stars it doesn't move with their group, so the Dipper slowly changes shape."],
  "Hadar": ["With Rigil Kentaurus, one of the two \"pointers\" that lead to the Southern Cross."],
  "Acrux": ["The brightest star of the Southern Cross, which appears on the flags of Australia, New Zealand, Brazil, Papua New Guinea and Samoa."],
  "Gacrux": ["The red giant at the top of the Southern Cross, and the nearest red giant to the Sun."],
  "Shaula": ["The stinger at the end of the scorpion's tail."],
  "Eltanin": ["By watching it from 1725 to 1728, James Bradley discovered the aberration of starlight: the first direct proof that Earth moves around the Sun."],
  "Ran": ["Epsilon Eridani, one of the nearest Sun-like stars.", "In 1960 it was one of the first two stars searched for alien radio signals, by Project Ozma."],
  "HIP 8102": ["Tau Ceti, one of the most Sun-like stars nearby.", "In 1960 it was one of the first two stars searched for alien radio signals, by Project Ozma."],
  "HIP 104214": ["61 Cygni was the first star to have its distance measured, by Friedrich Bessel in 1838."],
  "Alpheratz": ["It is shared by two star patterns: the head of Andromeda and a corner of the Great Square of Pegasus."],
  "Alphard": ["Its name means \"the solitary one\": it shines alone in an empty part of the sky."],
  "Denebola": ["Its name means \"the lion's tail\"."],
  "Sadr": ["The centre of the Northern Cross, in Cygnus."],
  "Albireo": ["Through a small telescope it splits into a striking pair: one gold, one blue."],
  "Adhara": ["About 4.7 million years ago it was much closer to us and was the brightest star in Earth's sky."],
  "Peacock": ["It got its name in the late 1930s, when Britain's Royal Air Force needed names for its navigation stars."],
  "Zubenelgenubi": ["Its name means \"southern claw\": it was once part of the scorpion's claws, before Libra became its own constellation."],
  "Vindemiatrix": ["Its name means \"grape harvester\": in ancient times it rose before the Sun around the start of the grape harvest."],
  "Alcyone": ["The brightest star in the Pleiades star cluster."],
  "Atlas": ["In Greek myth Atlas is the father of the Pleiades sisters, and his star sits beside them in the cluster."],
  "Alderamin": ["Earth's slow wobble will make it the north pole star around 7500 CE."],
  "Errai": ["It will be the north pole star around 4000 CE, and it has a giant planet."],
  "Kochab": ["With Pherkad it was called a \"guardian of the pole\", back when it was closer to the north pole than Polaris, around 1000 BCE."],
  "Mirfak": ["The brightest member of the Alpha Persei cluster, a family of young stars moving together."]
};
// When a star's light left it, matched to something happening on Earth (year, what happened). Edit freely.
const LIGHT_EVENTS = [
  [1969, "the first Moon landing"], [1957, "Malaya's independence and the launch of Sputnik, the first satellite"],
  [1903, "the Wright brothers' first powered flight"], [1879, "Edison's first long-lasting light bulb"],
  [1781, "William Herschel's discovery of Uranus"], [1687, "Isaac Newton's laws of motion"],
  [1609, "Galileo first pointing a telescope at the sky"], [1400, "the founding of the Melaka Sultanate"],
  [1054, "the supernova that made the Crab Nebula"], [79, "the eruption of Vesuvius that buried Pompeii"],
  [-221, "the first unification of China under Qin Shi Huang"], [-776, "the first recorded Olympic Games"],
  [-2560, "the building of the Great Pyramid of Giza"], [-3200, "the invention of writing in Mesopotamia"],
  [-9500, "the building of Göbekli Tepe, one of the oldest known temples"], [-15000, "the Lascaux cave paintings"]];
function lightEvent(leftYear) {
  const age = NOW_YEAR - leftYear; if (age < 40) return null;
  let best = null, bd = 0.12;                    // within about ±12% of the light's age
  for (const [y, what] of LIGHT_EVENTS) { const d = Math.abs(Math.log((NOW_YEAR - y) / age)); if (d < bd) { bd = d; best = [y, what]; } }
  return best;
}
const SPECT = (D.spect || []).map(x => x.replace(/\.+$/, ""));   // catalogue types sometimes end in "..."
function starType(sp) {               // spectral type -> plain words, e.g. "G2V" -> main-sequence star
  if (!sp) return null;
  const c = /^D/.test(sp) ? "D" : /III/.test(sp) ? "III" : /IV/.test(sp) ? "IV" : /II/.test(sp) ? "II" : /Ia|Ib|I(?![IV])/.test(sp) ? "I" : /V/.test(sp) ? "V" : null;
  return { D: "a white dwarf, the hot leftover core of a star like the Sun", I: "a supergiant, one of the biggest and brightest kinds of star",
           II: "a bright giant", III: "a giant: a star that has swollen up late in its life", IV: "a subgiant, starting to swell as it ages",
           V: "a main-sequence star, still burning hydrogen in its core like the Sun" }[c] || null;
}
let RANK_NOW = null;                  // brightness rank of every star today (1 = brightest)
function rankToday(i) {
  if (!RANK_NOW) {
    const t = [0, 0, 0, 0, 0], m = new Float32Array(N); for (let k = 0; k < N; k++) m[k] = posAt(k, NOW_YEAR - 2000, t)[4];
    const order = Array.from({ length: N }, (_, k) => k).sort((a, b) => m[a] - m[b]); RANK_NOW = new Int32Array(N); order.forEach((k, r) => RANK_NOW[k] = r + 1);
  }
  return RANK_NOW[i];
}
const ord = n => n + (n % 100 >= 11 && n % 100 <= 13 ? "th" : ["th", "st", "nd", "rd"][n % 10] || "th");
const fmtRAh = ra => { const h = ((ra / 15) % 24 + 24) % 24, H = Math.floor(h), M = Math.floor((h - H) * 60), Sx = ((h - H) * 60 - M) * 60; return `${H}h ${String(M).padStart(2, "0")}m ${Sx.toFixed(1).padStart(4, "0")}s`; };
const fmtDecd = d => { const a = Math.abs(d), Dg = Math.floor(a), M = Math.floor((a - Dg) * 60), Sx = Math.round(((a - Dg) * 60 - M) * 60); return `${d < 0 ? "−" : "+"}${Dg}° ${String(M).padStart(2, "0")}′ ${String(Sx).padStart(2, "0")}″`; };
const toRaDec = p => { const ra = (Math.atan2(p[1], p[0]) * 180 / Math.PI + 360) % 360; return [ra, Math.asin(p[2] / p[3]) * 180 / Math.PI]; };

let lastStarKey = "";
function renderStar() {
  const pane = document.getElementById("starDetail"), list = document.getElementById("starList"), i = state.sel;
  const key = i + ":" + Math.round(state.year);
  if (key === lastStarKey) return; lastStarKey = key;
  pane.hidden = i < 0; list.hidden = i >= 0;
  if (i < 0) { pane.innerHTML = ""; if (state.tab === "pane-star") { computeEpoch(state.year); renderNames(true); } return; }
  const o = i * 8, now = posAt(i, NOW_YEAR - 2000, [0, 0, 0, 0, 0]), then = posAt(i, state.year - 2000, [0, 0, 0, 0, 0]);
  const v = [S[o + 3], S[o + 4], S[o + 5]], speed = Math.hypot(...v) * PCYR2KMS;
  const radial = (then[0] * v[0] + then[1] * v[1] + then[2] * v[2]) / then[3] * PCYR2KMS;
  const cosA = (now[0] * then[0] + now[1] * then[1] + now[2] * then[2]) / (now[3] * then[3]);
  const moved = Math.acos(Math.min(1, cosA)) * 180 / Math.PI;
  const unknownDist = now[3] > 9e4;
  const vv = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  const tc = vv > 0 ? -(S[o] * v[0] + S[o + 1] * v[1] + S[o + 2] * v[2]) / vv : NaN;
  const dmin = vv > 0 ? Math.hypot(S[o] + v[0] * tc, S[o + 1] + v[1] * tc, S[o + 2] + v[2] * tc) * PC2LY : NaN;
  const eye = m => m <= 6 ? `<span class="pill eye">naked eye</span>` : `<span class="pill dim">needs binoculars</span>`;
  const ly = d => unknownDist ? "unknown" : (d * PC2LY < 100 ? (d * PC2LY).toFixed(1) : fmt(d * PC2LY)) + " ly";
  const same = Math.round(state.year) === NOW_YEAR;
  const leftYear = state.year - then[3] * PC2LY;
  const [raN, decN] = toRaDec(now), [raT, decT] = toRaDec(then);

  // fun facts: hand-written ones first, then ones worked out from the star's own numbers
  const facts = [...(STAR_FACTS[proper[i]] || STAR_FACTS["HIP " + HIP[i]] || [])];
  const type = starType(SPECT[i]);
  if (type) facts.push(`Its spectral type is <b>${SPECT[i]}</b>, which makes it ${type}.`);
  const T = 4600 * (1 / (0.92 * S[o + 7] + 1.7) + 1 / (0.92 * S[o + 7] + 0.62));
  if (isFinite(T) && T > 1500) facts.push(`${colorWord(S[o + 7])} light means a surface of about <b>${fmt(Math.round(T / 100) * 100)} K</b> (${fmt(Math.round((T - 273) / 100) * 100)} °C). The Sun's is about 5,800 K.`);
  if (!unknownDist) {
    const L = Math.pow(10, (4.83 - S[o + 6]) / 2.5);
    facts.push(L >= 1.5 && L < 2.5 ? "It gives off about <b>twice</b> as much visible light as the Sun." : L >= 1.5 ? `It gives off about <b>${fmt(L < 20 ? Math.round(L * 10) / 10 : Number(L.toPrecision(2)))} times</b> as much visible light as the Sun.`
      : L >= 0.67 ? "It gives off about as much visible light as the Sun." : 1 / L < 2.5 ? "It gives off only about <b>half</b> as much visible light as the Sun."
      : `It gives off only about <b>1/${fmt(Number((1 / L).toPrecision(2)))}</b> of the Sun's visible light.`);
  }
  const rNow = rankToday(i);
  let rankFact = rNow === 1 || rNow > 2000 ? "" : `The <b>${ord(rNow)}</b> brightest star in the night sky today.`;   // the map only holds stars brighter than mag 7, so deeper ranks would mislead
  if (!same && rNow <= 2000) {
    let rThen = 1; for (let k = 0; k < N; k++) if (MAG[k] < MAG[i]) rThen++;
    const when = `In ${yearLabel(state.year)} it ${state.year > NOW_YEAR ? "will be" : "was"}`;
    if (rThen !== rNow) rankFact += ` ${when} ${rThen === 1 ? "<b>the brightest star</b> in the night sky" : `the <b>${ord(rThen)}</b> brightest`}.`;
    else if (rNow === 1) rankFact = `${when} still the brightest star in the night sky.`;
  }
  if (rankFact) facts.push(rankFact.trim());
  if (!unknownDist && speed >= 1) facts.push(`At ${speed.toFixed(1)} km/s it could fly around the Earth in about <b>${(40075 / speed / 60) < 1 ? Math.round(40075 / speed) + " seconds" : Math.round(40075 / speed / 60) + " minutes"}</b>.`);
  const ev = lightEvent(leftYear);

  pane.innerHTML = `
    <div class="kicker">${conAbbr(i) ? conName(conAbbr(i)) : ""}${HIP[i] ? " · HIP " + HIP[i] : ""}</div>
    <h2>${starName(i)}</h2>
    <p class="note" style="margin-top:4px">${colorWord(S[o + 7])} star${SPECT[i] ? " · type " + SPECT[i] : ""}${unknownDist ? " · distance not measured" : ""}</p>
    <table class="compare">
      <tr><th></th><th>Today</th><th class="then">${yearLabel(state.year)}</th></tr>
      <tr><td>Right Ascension</td><td>${fmtRAh(raN)}<span class="deg">${raN.toFixed(4)}°</span></td><td class="then">${fmtRAh(raT)}<span class="deg">${raT.toFixed(4)}°</span></td></tr>
      <tr><td>Declination</td><td>${fmtDecd(decN)}<span class="deg">${decN.toFixed(4)}°</span></td><td class="then">${fmtDecd(decT)}<span class="deg">${decT.toFixed(4)}°</span></td></tr>
      <tr><td>Distance</td><td>${ly(now[3])}</td><td class="then">${ly(then[3])}</td></tr>
      <tr><td>Brightness (mag)</td><td>${now[4].toFixed(2)}</td><td class="then">${then[4].toFixed(2)}</td></tr>
      <tr><td>Visible?</td><td>${eye(now[4])}</td><td class="then">${eye(then[4])}</td></tr>
    </table>
    <p class="note" style="margin:2px 0 0">Coordinates are J2000 (the standard star-chart grid), moved along the star's own motion.</p>
    <h3 style="margin-top:16px">Fun facts</h3>
    <div class="facts">${facts.map(f => `<div class="fact">${f}</div>`).join("")}</div>
    <h3 style="margin-top:16px">On the move</h3>
    <div class="facts">
      ${unknownDist ? `<div class="fact">With no measured distance, this star stays fixed on the map.</div>` : `
      <div class="fact">Moves at <b>${speed.toFixed(1)} km/s</b> relative to the Sun, ${radial < -0.5 ? `<b>approaching</b> at ${(-radial).toFixed(1)} km/s` : radial > 0.5 ? `<b>moving away</b> at ${radial.toFixed(1)} km/s` : "mostly sideways"}${same ? "" : " in " + yearLabel(state.year)}.</div>
      ${same ? "" : `<div class="fact">Shifted <b>${moved < 1 ? (moved * 60).toFixed(1) + "′" : moved.toFixed(2) + "°"}</b> across the sky between today and ${yearLabel(state.year)}${moved > 0.5 ? ", about " + (moved / 0.5).toFixed(0) + " full-Moon widths" : ""}.</div>`}
      <div class="fact">Light seen from Earth in ${yearLabel(state.year)} left this star around <b>${yearLabel(leftYear)}</b>${ev ? `, around the time of ${ev[1]} (${yearLabel(ev[0])})` : ""}. You always see a star as it was, never as it is.</div>
      ${isFinite(tc) && Math.abs(tc) < 5e6 ? `<div class="fact">Closest to the Sun around <b>${yearLabel(2000 + tc)}</b>, at ${dmin < 100 ? dmin.toFixed(1) : fmt(dmin)} light-years.</div>` : ""}`}
    </div>
    <div class="row" style="margin-top:14px"><button class="btn" id="centreStar">Centre on this star</button><button class="btn" id="spherexStar">Look with SPHEREx</button></div>`;
  pane.insertAdjacentHTML("afterbegin", `<button class="btn backlist" id="backToList">← All stars</button>`);
  document.getElementById("centreStar").onclick = () => flyTo(Math.atan2(P[i * 3 + 1], P[i * 3]), Math.asin(P[i * 3 + 2]), Math.min(state.fov, 60));
  document.getElementById("backToList").onclick = () => { state.sel = -1; lastStarKey = ""; dirty = true; updateHUD(); };
  document.getElementById("spherexStar").onclick = () => {        // SPHEREx images are from 2025-26, so use the star's position then
    const p = posAt(i, 2025.9 - 2000, [0, 0, 0, 0, 0]), [ra, dec] = toRaDec(p);
    showTab("pane-spherex"); runTimelapse({ name: starName(i), ra, dec });
  };
}

// ---------- time dial ----------
const slider = document.getElementById("slider");
const sliderToYear = s => NOW_YEAR + Math.sign(s) * Math.pow(Math.abs(s) / 1000, 2) * SPAN;
const yearToSlider = y => { const d = y - NOW_YEAR; return Math.sign(d) * Math.sqrt(Math.min(1, Math.abs(d) / SPAN)) * 1000; };
slider.addEventListener("input", () => { stopAll(); setYear(sliderToYear(+slider.value)); });
(function ticks() {
  const el = document.getElementById("ticks");
  for (const d of [-200000, -100000, -50000, -10000, 0, 10000, 50000, 100000, 200000]) {
    const s = document.createElement("span"); s.style.left = (50 + yearToSlider(NOW_YEAR + d) / 20) + "%";
    s.textContent = d === 0 ? "now" : (d > 0 ? "+" : "−") + Math.abs(d / 1000) + "k"; if (!d) s.className = "now"; el.appendChild(s);
  }
})();
function setYear(y) { state.year = Math.max(NOW_YEAR - SPAN, Math.min(NOW_YEAR + SPAN, y)); dirty = true; hudDirty = true; }
let hudDirty = true;
const btnF = document.getElementById("playFwd"), btnB = document.getElementById("playBack");
function stopAll() { playing = 0; journey = null; if (typeof endAtlasAnim === "function") { endAtlasAnim(); endM101Anim(); } btnF.classList.remove("on"); btnB.classList.remove("on"); btnF.textContent = "▶"; btnB.textContent = "◀"; }
function play(dir) {
  if (playing === dir) { stopAll(); return; }
  stopAll(); playing = dir; (dir > 0 ? btnF : btnB).classList.add("on"); (dir > 0 ? btnF : btnB).textContent = "❚❚";
}
btnF.onclick = () => play(1); btnB.onclick = () => play(-1);
document.getElementById("today").onclick = () => { stopAll(); setYear(NOW_YEAR); };
function parseYear(txt) {
  const t = txt.trim().toUpperCase().replace(/,/g, ""); const m = /^(-?\d+(?:\.\d+)?)\s*(K)?\s*(BCE|BC|CE|AD)?$/.exec(t);
  if (!m) return null; let n = parseFloat(m[1]) * (m[2] ? 1000 : 1);
  if (m[3] === "BCE" || m[3] === "BC" || n < 0) n = 1 - Math.abs(n);
  return n;
}
const goYear = () => { const y = parseYear(document.getElementById("yearInput").value); const inp = document.getElementById("yearInput");
  if (y === null) { inp.style.borderColor = "var(--past)"; inp.title = "Type a year like 14000, 2600 BCE or 50k"; return; }
  inp.style.borderColor = ""; animateYear(y); };
document.getElementById("goYear").onclick = goYear;
document.getElementById("yearInput").addEventListener("keydown", e => { if (e.key === "Enter") goYear(); });

// ---------- view interaction ----------
const pointers = new Map(); let dragStart = null, pinchDist = 0, moved = 0;
cv.addEventListener("pointerdown", e => {
  cv.setPointerCapture(e.pointerId); pointers.set(e.pointerId, [e.clientX, e.clientY]);
  if (pointers.size === 1) { dragStart = [e.clientX, e.clientY]; moved = 0; cv.classList.add("dragging"); }
  if (pointers.size === 2) { const [a, b] = [...pointers.values()]; pinchDist = Math.hypot(a[0] - b[0], a[1] - b[1]); }
});
cv.addEventListener("pointermove", e => {
  if (!pointers.has(e.pointerId)) return;
  const prev = pointers.get(e.pointerId); pointers.set(e.pointerId, [e.clientX, e.clientY]);
  if (pointers.size === 2) {
    const [a, b] = [...pointers.values()], d = Math.hypot(a[0] - b[0], a[1] - b[1]);
    if (pinchDist) zoomAt(pinchDist / d, (a[0] + b[0]) / 2, (a[1] + b[1]) / 2); pinchDist = d; moved += 10; return;
  }
  const dx = e.clientX - prev[0], dy = e.clientY - prev[1]; moved += Math.abs(dx) + Math.abs(dy);
  if (journey) journey.view = null; if (atlasAnim) atlasAnim.follow = false; if (m101Anim) m101Anim.follow = false;
  state.ra = (state.ra + dx / (scale * Math.max(0.25, Math.cos(state.dec))) + 2 * Math.PI) % (2 * Math.PI);
  state.dec = Math.max(-1.55, Math.min(1.55, state.dec + dy / scale)); dirty = true;
});
function endPointer(e) {
  pointers.delete(e.pointerId); if (pointers.size < 2) pinchDist = 0;
  if (pointers.size === 0) { cv.classList.remove("dragging"); if (moved < 6 && dragStart) pick(e.clientX, e.clientY); dragStart = null; }
}
cv.addEventListener("pointerup", endPointer); cv.addEventListener("pointercancel", endPointer);
// ---------- zoom: toward the mouse, faster steps, buttons, keys, double-click ----------
const FOV_MIN = 0.5, FOV_MAX = 200;
function unproject(x, y) {      // screen point -> sky direction (inverse of the stereographic projection)
  setBasis(); const a = (x - cx) / scale, b = (cy - y) / scale, r2 = a * a + b * b, c = (4 - r2) / (4 + r2), k = (1 + c) / 2;
  return [0, 1, 2].map(i => c * F[i] + a * k * R[i] + b * k * U[i]);
}
const raDecOf = v => [Math.atan2(v[1], v[0]), Math.asin(Math.max(-1, Math.min(1, v[2])))];
function zoomAt(f, x, y) {      // zoom by factor f, keeping the sky point under (x, y) where it is
  if (atlasAnim) atlasAnim.follow = false; if (m101Anim) m101Anim.follow = false; if (journey) journey.view = null;
  const fov = Math.max(FOV_MIN, Math.min(FOV_MAX, state.fov * f)); if (fov === state.fov) return;
  const keep = x == null ? null : unproject(x, y); state.fov = fov;
  if (keep) for (let it = 0; it < 4; it++) {
    const [r0, d0] = raDecOf(keep), [r1, d1] = raDecOf(unproject(x, y));
    state.ra = (state.ra + (((r0 - r1 + 3 * Math.PI) % (2 * Math.PI)) - Math.PI) + 2 * Math.PI) % (2 * Math.PI);
    state.dec = Math.max(-1.55, Math.min(1.55, state.dec + (d0 - d1)));
  }
  dirty = true;
}
function zoomBy(f) { zoomAt(f, null, null); }
function zoomSmooth(f) {      // buttons and keys: a short eased zoom around the centre
  const fov = Math.max(FOV_MIN, Math.min(FOV_MAX, state.fov * f));
  if (atlasAnim) atlasAnim.follow = false; if (m101Anim) m101Anim.follow = false;
  if (reduceMotion) { state.fov = fov; dirty = true; return; }
  journey = journey || {}; journey.view = { t0: performance.now(), dur: 320, ra0: state.ra, dra: 0, dec0: state.dec, dec1: state.dec, fov0: state.fov, fov1: fov };
}
cv.addEventListener("wheel", e => {
  e.preventDefault();
  const px = e.deltaMode === 1 ? e.deltaY * 18 : e.deltaMode === 2 ? e.deltaY * 400 : e.deltaY;      // lines/pages -> pixels
  const k = e.ctrlKey ? 0.012 : 0.0028;      // ctrlKey = touchpad pinch (small deltas), otherwise a mouse wheel or two-finger scroll
  zoomAt(Math.exp(Math.max(-0.6, Math.min(0.6, px * k))), e.clientX, e.clientY);
}, { passive: false });
cv.addEventListener("dblclick", e => { e.preventDefault(); zoomAt(0.5, e.clientX, e.clientY); });
document.getElementById("zoomIn").onclick = () => zoomSmooth(0.5);
document.getElementById("zoomOut").onclick = () => zoomSmooth(2);
window.addEventListener("keydown", e => {
  if (e.target.closest && e.target.closest("input, select, textarea")) return;
  if (e.key === "+" || e.key === "=") zoomSmooth(0.6); else if (e.key === "-" || e.key === "_") zoomSmooth(1 / 0.6);
});
const inCut = (x, y) => cutRect && Math.abs(x - cutRect.x) <= cutRect.half + 4 && Math.abs(y - cutRect.y) <= cutRect.half + 4;
const spotAt = (x, y) => { const h = spotHits.find(h => x >= h.x0 && x <= h.x1 && y >= h.y0 && y <= h.y1); return h ? h.sp : null; };
cv.addEventListener("pointermove", e => { if (!pointers.size) cv.style.cursor = inCut(e.clientX, e.clientY) || spotAt(e.clientX, e.clientY) ? "pointer" : ""; });
function openSpot(sp) {      // a tutorial marker: go there and open its explanation in the SPHEREx tab
  showTab("pane-spherex");
  if (sp.m101) startM101Anim();
  else { const [ra, de] = raDecOf(sp.v); flyTo(ra, de, sp.fov || 10); }
  const h = document.getElementById(sp.sec); if (!h) return;
  if (h.tagName === "DETAILS") h.open = true;
  h.scrollIntoView({ block: "start", behavior: reduceMotion ? "auto" : "smooth" });
  h.classList.remove("flash"); void h.offsetWidth; h.classList.add("flash");
}
function pick(x, y) {
  if (inCut(x, y)) { if (cutBig) document.getElementById("cpPlay").click(); else openCutPop(); return; }
  const sp = spotAt(x, y); if (sp) { openSpot(sp); return; }
  let best = -1, bestScore = 1e9;
  for (let i = 0; i < N; i++) {
    if (!VIS[i]) continue; const d = Math.hypot(SX[i] - x, SY[i] - y); if (d > 16) continue;
    const s = d - SR[i] * 2.5; if (s < bestScore) { bestScore = s; best = i; }
  }
  state.sel = best; dirty = true; lastStarKey = ""; updateHUD();
  if (best >= 0) showTab("pane-star");
}
window.addEventListener("keydown", e => {
  if (e.target.tagName === "INPUT" || e.target.tagName === "SELECT") return;
  if (e.key === "/") { e.preventDefault(); sInput.focus(); return; }
  if (e.code === "Space") { e.preventDefault(); play(playing || 1); }
  if (e.key === "ArrowRight") { stopAll(); setYear(state.year + 1000); }
  if (e.key === "ArrowLeft") { stopAll(); setYear(state.year - 1000); }
});

// ---------- animation ----------
const reduceMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
const ease = t => t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2;
function flyTo(ra, dec, fov) {
  if (reduceMotion) { Object.assign(state, { ra: (ra + 2 * Math.PI) % (2 * Math.PI), dec, fov }); dirty = true; return; }
  let dra = ((ra - state.ra + 3 * Math.PI) % (2 * Math.PI)) - Math.PI;
  journey = journey || {}; journey.view = { t0: performance.now(), dur: 1300, ra0: state.ra, dra, dec0: state.dec, dec1: dec, fov0: state.fov, fov1: fov };
}
function animateYear(y, delay = 0) {
  playing = 0; btnF.classList.remove("on"); btnB.classList.remove("on"); btnF.textContent = "▶"; btnB.textContent = "◀";
  if (reduceMotion) { setYear(y); return; }
  journey = journey || {}; const span = Math.abs(y - state.year);
  journey.time = { t0: performance.now() + delay, dur: Math.min(6000, 1800 + Math.sqrt(span) * 18), y0: state.year, y1: y };
}
let last = performance.now();
function frame(now) {
  const dtSec = Math.min(0.1, (now - last) / 1000); last = now;
  if (playing) {
    setYear(state.year + playing * +document.getElementById("speed").value * dtSec);
    if (state.year >= NOW_YEAR + SPAN || state.year <= NOW_YEAR - SPAN) stopAll();
  }
  if (journey) {
    const v = journey.view;
    if (v) { const t = Math.min(1, (now - v.t0) / v.dur), e = ease(t);
      state.ra = (v.ra0 + v.dra * e + 2 * Math.PI) % (2 * Math.PI); state.dec = v.dec0 + (v.dec1 - v.dec0) * e; state.fov = v.fov0 + (v.fov1 - v.fov0) * e; dirty = true;
      if (t >= 1) journey.view = null; }
    const tm = journey.time;
    if (tm && now >= tm.t0) { const t = Math.min(1, (now - tm.t0) / tm.dur); setYear(tm.y0 + (tm.y1 - tm.y0) * ease(t)); if (t >= 1) journey.time = null; }
    if (!journey.view && !journey.time) journey = null;
  }
  if (atlasAnim) stepAtlasAnim(now, dtSec);
  if (m101Anim) stepM101Anim(dtSec);
  stepCutGrow(dtSec);
  if (dirty) { draw(); dirty = false; placeCutBar(); }
  if (hudDirty) { updateHUD(); hudDirty = false; }
  requestAnimationFrame(frame);
}

// ---------- journeys ----------
const findStar = n => proper.indexOf(n);
const JOURNEYS = [
  { title: "The Big Dipper bends", to: NOW_YEAR + 100000, ra: 185, dec: 56, fov: 55, sel: "Alkaid", trails: true,
    d: "Five of its seven stars travel together as one family. Dubhe and Alkaid head the other way, so the handle and bowl pull apart." },
  { title: "Thuban, pole star of the pyramid builders", to: -2599, ra: 230, dec: 78, fov: 75, sel: "Thuban", poles: true,
    d: "When the Great Pyramid was built around 2600 BCE, Earth's axis pointed near Thuban in Draco, not Polaris." },
  { title: "Vega takes over the north", to: 14000, ra: 275, dec: 62, fov: 85, sel: "Vega", poles: true,
    d: "Earth's wobble swings the north pole toward Vega, the bright star of summer evenings. Polaris will no longer mark north." },
  { title: "Alpha Centauri comes closest", to: null, ra: 218, dec: -61, fov: 60, sel: "Rigil Kentaurus", trails: true,
    d: "Our nearest bright neighbour is closing in. It will peak in brightness and then sweep on through the southern sky." },
  { title: "Sky of the Lascaux cave painters", to: -15000, ra: 72, dec: 18, fov: 80, sel: "Aldebaran",
    d: "About 17,000 years ago people painted animals in Lascaux cave, France. Some researchers read a dot cluster there as the Pleiades." },
  { title: "Interstellar comet 3I/ATLAS", to: NOW_YEAR, ra: 209, dec: -6, fov: 120, sel: "", atlas: true,
    d: "A comet from another star system, photographed by SPHEREx in August and December 2025. Ride along its real track, image by image." },
  { title: "Into the Pinwheel Galaxy", to: NOW_YEAR, ra: 210.8, dec: 54.3, fov: 50, sel: "", m101: true,
    d: "Start at the Big Dipper's handle and zoom 150 times into a galaxy 21 million light-years away, then watch SPHEREx's real images of it." },
  { title: "Arcturus, a visitor from the galaxy's halo", to: NOW_YEAR - 150000, ra: 214, dec: 19, fov: 90, sel: "Arcturus", trails: true,
    d: "Arcturus moves at about 120 km/s. Run the clock back and watch it slide across Boötes, much faster than its neighbours." }
];
(function buildJourneys() {
  const ac = findStar("Rigil Kentaurus");
  if (ac >= 0) { const o = ac * 8, v = [S[o + 3], S[o + 4], S[o + 5]], vv = v[0] ** 2 + v[1] ** 2 + v[2] ** 2;
    JOURNEYS[3].to = 2000 - (S[o] * v[0] + S[o + 1] * v[1] + S[o + 2] * v[2]) / vv; }
  const box = document.getElementById("journeyList");
  JOURNEYS.forEach(j => {
    const b = document.createElement("button"); b.className = "journey";
    b.innerHTML = `<span class="to">→ ${yearLabel(j.to)}</span><h3>${j.title}</h3><span class="d">${j.d}</span>`;
    b.onclick = () => {
      stopAll(); if (j.trails) setLayer("trails", true); if (j.poles) setLayer("poles", true); if (j.atlas) setLayer("atlas", true);
      if (j.atlas) { startAtlasAnim(); hudDirty = true; return; }
      if (j.m101) { startM101Anim(); hudDirty = true; return; }
      const i = findStar(j.sel); state.sel = i; lastStarKey = "";
      flyTo(j.ra * Math.PI / 180, j.dec * Math.PI / 180, j.fov); animateYear(j.to, 700); hudDirty = true;
    };
    box.appendChild(b);
  });
})();

// ---------- layers & search ----------
const LAYER_LABELS = { lines: "Constellation lines", conNames: "Constellation names", names: "Star names", milky: "Milky Way",
  trails: "Motion trails (next 30,000 years)", poles: "Celestial poles and their wobble", grid: "Coordinate grid", spherex: "SPHEREx example spots", atlas: "Comet 3I/ATLAS track (2025)" };
function setLayer(k, v) { state.layers[k] = v; const el = document.getElementById("layer-" + k); if (el) el.checked = v; dirty = true; }
(function buildToggles() {
  const box = document.getElementById("toggles");
  for (const [k, label] of Object.entries(LAYER_LABELS)) {
    const l = document.createElement("label"); l.className = "toggle"; l.htmlFor = "layer-" + k;
    l.innerHTML = `<span>${label}</span><input type="checkbox" id="layer-${k}" ${state.layers[k] ? "checked" : ""}>`;
    l.querySelector("input").onchange = e => setLayer(k, e.target.checked); box.appendChild(l);
  }
})();

// ---------- search bar & star name list ----------
const norm = t => t.normalize("NFD").replace(/[̀-ͯ]/g, "").replace(/[’'`.\-]/g, "").toLowerCase().trim();
const NAMED = proper.map((n, i) => n ? i : -1).filter(i => i >= 0);
const NAMED_KEY = new Map(NAMED.map(i => [i, norm(proper[i])]));
const CONS = Object.keys(conStars).filter(a => conStars[a].length).map(a => ({ abbr: a, name: conName(a), key: norm(conName(a)) }));
function selectStar(i) {
  state.sel = i; lastStarKey = ""; hudDirty = true; dirty = true;
  flyTo(Math.atan2(P[i * 3 + 1], P[i * 3]), Math.asin(P[i * 3 + 2]), Math.min(state.fov, 60));
  showTab("pane-star");
}
function goConstellation(abbr) {
  let x = 0, y = 0, z = 0; for (const s of conStars[abbr]) { x += P[s * 3]; y += P[s * 3 + 1]; z += P[s * 3 + 2]; }
  const n = Math.hypot(x, y, z); flyTo(Math.atan2(y, x), Math.asin(z / n), 60);
}
function findMatches(q) {
  const k = norm(q); if (!k) return [];
  const out = [];
  const hip = /^hip\s*(\d+)$/i.exec(q.trim());
  if (hip) { const i = HIP.indexOf(+hip[1]); if (i >= 0) out.push({ type: "star", i, label: starName(i), rank: 0 }); }
  for (const i of NAMED) { const key = NAMED_KEY.get(i), p = key.indexOf(k); if (p >= 0) out.push({ type: "star", i, label: proper[i], rank: (p === 0 ? 0 : 2) + MAG[i] / 100 }); }
  if ("3i/atlas comet interstellar".includes(k) || k.includes("atlas") || k === "3i") out.push({ type: "atlas", label: "3I/ATLAS", rank: 0.2 });
  for (const c of CONS) { const p = c.key.indexOf(k); if (p >= 0 || norm(c.abbr) === k) out.push({ type: "con", abbr: c.abbr, label: c.name, rank: p === 0 || norm(c.abbr) === k ? 0.5 : 2.5 }); }
  return out.sort((x, y) => x.rank - y.rank).slice(0, 10);
}
const sInput = document.getElementById("starSearch"), sList = document.getElementById("suggest");
let sItems = [], sActive = -1;
function esc(t) { return t.replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c])); }
function highlight(label, q) {
  const t = q.trim(), i = label.toLowerCase().indexOf(t.toLowerCase());
  return i >= 0 && t ? esc(label.slice(0, i)) + "<mark>" + esc(label.slice(i, i + t.length)) + "</mark>" + esc(label.slice(i + t.length)) : esc(label);
}
function showSuggest() {
  const q = sInput.value; sItems = findMatches(q); sActive = sItems.length ? 0 : -1;
  if (!q.trim()) { sList.hidden = true; sInput.setAttribute("aria-expanded", "false"); return; }
  sList.innerHTML = sItems.length ? sItems.map((m, k) => `<li role="option" id="sug-${k}" data-k="${k}" aria-selected="${k === sActive}"><span>${highlight(m.label, q)}</span><span class="what">${m.type === "atlas" ? "interstellar comet" : m.type === "con" ? "constellation" : conName(conAbbr(m.i)) + " · mag " + MAG[m.i].toFixed(1)}</span></li>`).join("")
    : `<li class="none">No star or constellation called “${esc(q.trim())}”</li>`;
  sList.hidden = false; sInput.setAttribute("aria-expanded", "true");
  if (state.tab === "pane-star" && state.sel < 0) renderNames(true);
}
function pickSuggest(k) {
  const m = sItems[k]; if (!m) return;
  sList.hidden = true; sInput.setAttribute("aria-expanded", "false"); sInput.value = m.label;
  if (m.type === "atlas") { document.getElementById("atlasAll").click(); return; }
  m.type === "star" ? selectStar(m.i) : goConstellation(m.abbr);
}
sInput.addEventListener("input", showSuggest);
sInput.addEventListener("focus", () => { if (sInput.value.trim()) showSuggest(); });
sInput.addEventListener("keydown", e => {
  if (e.key === "ArrowDown" || e.key === "ArrowUp") {
    e.preventDefault(); if (!sItems.length) return;
    sActive = (sActive + (e.key === "ArrowDown" ? 1 : -1) + sItems.length) % sItems.length;
    sList.querySelectorAll("li").forEach((li, k) => li.setAttribute("aria-selected", k === sActive));
    sInput.setAttribute("aria-activedescendant", "sug-" + sActive);
  } else if (e.key === "Enter") { e.preventDefault(); pickSuggest(Math.max(0, sActive)); }
  else if (e.key === "Escape") { sList.hidden = true; sInput.setAttribute("aria-expanded", "false"); }
});
sList.addEventListener("mousedown", e => { const li = e.target.closest("li[data-k]"); if (li) { e.preventDefault(); pickSuggest(+li.dataset.k); } });
document.addEventListener("pointerdown", e => { if (!e.target.closest(".searchbar")) { sList.hidden = true; sInput.setAttribute("aria-expanded", "false"); } });

let lastNamesRender = 0, lastNamesKey = "";
function renderNames(force) {
  const now = performance.now();
  const sort = document.getElementById("namesSort").value, q = norm(sInput.value);
  const key = sort + "|" + q + "|" + Math.round(state.year / 50) + "|" + state.sel;
  if (!force && key === lastNamesKey) return;
  if (!force && now - lastNamesRender < 400) { clearTimeout(renderNames.t); renderNames.t = setTimeout(() => renderNames(false), 420); return; }
  lastNamesKey = key; lastNamesRender = now;
  let list = NAMED.filter(i => !q || NAMED_KEY.get(i).includes(q) || norm(conName(conAbbr(i))).includes(q));
  const cmp = { mag: (a, b) => MAG[a] - MAG[b], az: (a, b) => proper[a].localeCompare(proper[b]), dist: (a, b) => DIST[a] - DIST[b],
    con: (a, b) => conName(conAbbr(a)).localeCompare(conName(conAbbr(b))) || MAG[a] - MAG[b] }[sort];
  list.sort(cmp);
  const y = Math.round(state.year) === NOW_YEAR ? "today" : "in " + yearLabel(state.year);
  document.getElementById("namesCount").textContent = `${list.length} named star${list.length === 1 ? "" : "s"}${q ? " matching" : ""} · values ${y}`;
  const ly = d => d > 9e4 ? "—" : d * PC2LY < 100 ? (d * PC2LY).toFixed(1) : fmt(d * PC2LY);
  document.getElementById("nameList").innerHTML =
    `<li class="cols"><span>Name</span><span>Constellation</span><span class="num">Mag</span><span class="num">Light-yr</span></li>` +
    list.map(i => `<li><button data-i="${i}" class="${i === state.sel ? "sel " : ""}${MAG[i] > 6 ? "faint" : ""}" title="${MAG[i] > 6 ? "Too faint to see without binoculars in this year" : ""}"><span class="nm">${esc(proper[i])}</span><span class="cn">${esc(conName(conAbbr(i)))}</span><span class="num">${MAG[i].toFixed(2)}</span><span class="num">${ly(DIST[i])}</span></button></li>`).join("");
}
document.getElementById("nameList").addEventListener("click", e => { const b = e.target.closest("button[data-i]"); if (b) selectStar(+b.dataset.i); });
document.getElementById("namesSort").addEventListener("change", () => renderNames(true));

// ---------- tabs ----------
function showTab(id) {
  state.tab = id;
  document.querySelectorAll(".tabs button").forEach(b => b.setAttribute("aria-selected", b.dataset.pane === id));
  document.querySelectorAll(".pane").forEach(p => p.classList.toggle("on", p.id === id));
  if (id === "pane-star" && state.sel < 0) { computeEpoch(state.year); renderNames(true); }
}
document.querySelectorAll(".tabs button").forEach(b => b.onclick = () => showTab(b.dataset.pane));

// ---------- FITS reader (SPHEREx cutouts) ----------
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
const d2r = Math.PI / 180;
function makeWCS(h) {
  if (!h.CTYPE1 || !String(h.CTYPE1).includes("TAN") || h.CRVAL1 === undefined) return null;
  let cd;
  if (h.CD1_1 !== undefined) cd = [h.CD1_1, h.CD1_2 || 0, h.CD2_1 || 0, h.CD2_2];
  else { const c1 = h.CDELT1 || 1, c2 = h.CDELT2 || 1; cd = [(h.PC1_1 ?? 1) * c1, (h.PC1_2 || 0) * c1, (h.PC2_1 || 0) * c2, (h.PC2_2 ?? 1) * c2]; }
  const det = cd[0] * cd[3] - cd[1] * cd[2], inv = [cd[3] / det, -cd[1] / det, -cd[2] / det, cd[0] / det];
  const sipTerms = (p) => { const t = []; const ord = h[p + "_ORDER"] || 0;
    for (let i = 0; i <= ord; i++) for (let j = 0; j <= ord - i; j++) { const c = h[`${p}_${i}_${j}`]; if (c) t.push([i, j, c]); } return t; };
  const A = sipTerms("A"), B = sipTerms("B");
  const poly = (t, u, v) => t.reduce((s, [i, j, c]) => s + c * Math.pow(u, i) * Math.pow(v, j), 0);
  const ra0 = h.CRVAL1 * d2r, de0 = h.CRVAL2 * d2r;
  return {
    pix2sky(x, y) {               // 1-based FITS pixel → [ra, dec] radians
      const u = x - h.CRPIX1, v = y - h.CRPIX2, uu = u + poly(A, u, v), vv = v + poly(B, u, v);
      const xi = (cd[0] * uu + cd[1] * vv) * d2r, eta = (cd[2] * uu + cd[3] * vv) * d2r;
      const den = Math.cos(de0) - eta * Math.sin(de0);
      return [ra0 + Math.atan2(xi, den), Math.atan((Math.sin(de0) + eta * Math.cos(de0)) / Math.hypot(xi, den))];
    },
    sky2pix(ra, dec) {
      const cc = Math.sin(de0) * Math.sin(dec) + Math.cos(de0) * Math.cos(dec) * Math.cos(ra - ra0);
      const xi = Math.cos(dec) * Math.sin(ra - ra0) / cc / d2r, eta = (Math.cos(de0) * Math.sin(dec) - Math.sin(de0) * Math.cos(dec) * Math.cos(ra - ra0)) / cc / d2r;
      const U0 = inv[0] * xi + inv[1] * eta, V0 = inv[2] * xi + inv[3] * eta; let u = U0, v = V0;
      for (let k = 0; k < 6; k++) { const nu = U0 - poly(A, u, v), nv = V0 - poly(B, u, v); u = nu; v = nv; }
      return [u + h.CRPIX1, v + h.CRPIX2];
    }
  };
}
function obsDate(...hdrs) {
  for (const h of hdrs) {
    for (const k of ["DATE-AVG", "DATE-OBS", "DATE-BEG"]) if (typeof h[k] === "string" && /^\d{4}-/.test(h[k])) return new Date(h[k].endsWith("Z") ? h[k] : h[k] + "Z");
    for (const k of ["MJD-AVG", "MJD-OBS", "MJD-BEG", "MJD_AVG"]) if (typeof h[k] === "number") return new Date((h[k] - 40587) * 86400000);
  }
  return null;
}
function readBinTable(buf, h) {
  const H = h.hdr, dv = new DataView(buf), size = { L: 1, B: 1, I: 2, J: 4, K: 8, A: 1, E: 4, D: 8, C: 8, M: 16, P: 8, Q: 16 };
  const cols = []; let off = 0;
  for (let k = 1; k <= (H.TFIELDS || 0); k++) {
    const m = /^(\d*)([LXBIJKAEDCMPQ])/.exec(String(H["TFORM" + k]).trim()); if (!m) return [];
    const rep = m[1] === "" ? 1 : +m[1], t = m[2];
    cols.push({ name: String(H["TTYPE" + k] || "col" + k).trim(), t, rep, off });
    off += t === "X" ? Math.ceil(rep / 8) : size[t] * rep;
  }
  const rows = [];
  for (let r = 0; r < H.NAXIS2; r++) {
    const base = h.off + r * H.NAXIS1, row = {};
    for (const c of cols) {
      const o = base + c.off;
      if (c.t === "A") { let t = ""; for (let i = 0; i < c.rep; i++) { const ch = dv.getUint8(o + i); if (!ch) break; t += String.fromCharCode(ch); } row[c.name] = t.trim(); }
      else if (c.t === "J") row[c.name] = dv.getInt32(o, false); else if (c.t === "I") row[c.name] = dv.getInt16(o, false);
      else if (c.t === "K") row[c.name] = Number(dv.getBigInt64(o, false)); else if (c.t === "E") row[c.name] = dv.getFloat32(o, false);
      else if (c.t === "D") row[c.name] = dv.getFloat64(o, false); else if (c.t === "B") row[c.name] = dv.getUint8(o);
    }
    rows.push(row);
  }
  return rows;
}
const mjdDate = v => !Number.isFinite(v) ? null : v > 2400000 ? new Date((v - 2440587.5) * 864e5) : v > 30000 ? new Date((v - 40587) * 864e5) : null;
// Every image in a file: a single SPHEREx image/cutout, or IRSA's multi-cutout file (IMAGE1, IMAGE2 … + CUTOUT_INFO).
function loadFITS(buf, name) {
  const hdus = parseFITS(buf);
  const info = hdus.find(x => x.hdr.XTENSION === "BINTABLE" && x.hdr.EXTNAME === "CUTOUT_INFO");
  const rows = info ? readBinTable(buf, info) : [];
  const imgs = hdus.filter(x => x.hdr.NAXIS === 2 && x.hdr.XTENSION !== "BINTABLE" && x.hdr.NAXIS1 > 1);
  let pick = imgs.filter(x => /^IMAGE\d*$/.test(x.hdr.EXTNAME || ""));
  if (!pick.length) pick = imgs.slice(0, 1);
  if (!pick.length) throw new Error(`${name} has no 2-D image inside. Use a SPHEREx spectral image or cutout file.`);
  return pick.slice(0, 300).map(h => {
    const suffix = (h.hdr.EXTNAME || "").replace(/^IMAGE/, ""), im = readImage(buf, h);
    const zodi = imgs.find(x => x.hdr.EXTNAME === "ZODI" + suffix && x.hdr.NAXIS1 === im.w && x.hdr.NAXIS2 === im.h);
    if (zodi) { const z = readImage(buf, zodi).data; for (let i = 0; i < im.data.length; i++) im.data[i] -= z[i]; }
    const row = suffix ? rows.find(r => String(r.cutout_index) === suffix) : null;
    const date = (row && mjdDate(row.observation_date)) || obsDate(h.hdr, hdus[0].hdr);
    return { ...im, name: suffix ? `${name} #${suffix}` : name, wcs: makeWCS(h.hdr), date, wl: row && Number.isFinite(row.central_wavelength) ? row.central_wavelength : null, zodi: !!zodi };
  });
}

// ---------- SPHEREx viewer: blink two dates or play many ----------
const fc = document.getElementById("fitsCanvas"), fx = fc.getContext("2d");
const fitsState = { frames: [], aligned: [], mode: "play", k: 0, sim: true, dmax: 1, stretch: [] };
function percentile(arr, p) { const a = Array.from(arr).filter(Number.isFinite).sort((x, y) => x - y); return a.length ? a[Math.min(a.length - 1, Math.floor(p / 100 * a.length))] : 0; }
function alignAll() {
  const F = fitsState.frames, A = F[0]; fitsState.aligned = [];
  for (const B of F) {
    if (B === A || (!A.wcs && B.w === A.w && B.h === A.h)) { fitsState.aligned.push(B.data); continue; }
    const out = new Float32Array(A.w * A.h).fill(NaN), d = B.data, w = B.w;
    for (let j = 0; j < A.h; j++) for (let i = 0; i < A.w; i++) {
      let x, y;
      if (A.wcs && B.wcs) { const [ra, de] = A.wcs.pix2sky(i + 1, j + 1); [x, y] = B.wcs.sky2pix(ra, de); x -= 1; y -= 1; }
      else { x = i * (B.w / A.w); y = j * (B.h / A.h); }
      const x0 = Math.floor(x), y0 = Math.floor(y); if (x0 < 0 || y0 < 0 || x0 >= B.w - 1 || y0 >= B.h - 1) continue;
      const fx_ = x - x0, fy = y - y0, o = y0 * w + x0;
      out[j * A.w + i] = d[o] * (1 - fx_) * (1 - fy) + d[o + 1] * fx_ * (1 - fy) + d[o + w] * (1 - fx_) * fy + d[o + w + 1] * fx_ * fy;
    }
    fitsState.aligned.push(out);
  }
  restretch();
}
function restretch() {
  const L = fitsState.aligned; if (!L.length) return;
  const p = +document.getElementById("contrast").value;
  fitsState.stretch = L.map(a => fitsState.real === "m101" ? [percentile(a, 2), percentile(a, 99.9)] : [percentile(a, 20), percentile(a, p)]);   // each frame gets its own stretch: wavelengths differ in brightness
  if (L.length > 1) { const a = L[0], b = L[L.length - 1], d = new Float32Array(a.length); for (let i = 0; i < d.length; i++) d[i] = Math.abs(b[i] - a[i]); fitsState.dmax = percentile(d, p) || 1; }
  renderFits();
}
function viewIndex() { const n = fitsState.aligned.length; return { play: fitsState.k % n, first: 0, last: n - 1, frame: fitsState.k % n }[fitsState.mode]; }
function renderFits() {
  const F = fitsState.frames, L = fitsState.aligned; if (!L.length) return;
  const A = F[0], n = L.length;
  if (fc.width !== A.w || fc.height !== A.h) { fc.width = A.w; fc.height = A.h; }
  const img = fx.createImageData(A.w, A.h), px = img.data, diff = fitsState.mode === "diff" && n > 1;
  const k = diff ? -1 : viewIndex(), src = diff ? null : L[k], [lo, hi] = diff ? [0, 1] : fitsState.stretch[k], ash = Math.asinh(10);
  for (let j = 0; j < A.h; j++) for (let i = 0; i < A.w; i++) {
    const q = j * A.w + i, o = ((A.h - 1 - j) * A.w + i) * 4;
    if (diff) {
      const d = L[n - 1][q] - L[0][q], t = Number.isFinite(d) ? Math.max(-1, Math.min(1, d / fitsState.dmax)) : 0;
      if (t >= 0) { px[o] = 95 * t; px[o + 1] = 179 * t; px[o + 2] = 212 * t; } else { px[o] = 235 * -t; px[o + 1] = 154 * -t; px[o + 2] = 92 * -t; }
    } else {
      const v = src[q], u = Math.max(0, (v - lo) / (hi - lo || 1));
      const t = !Number.isFinite(v) ? 0 : fitsState.real === "m101" ? Math.min(1, u) : Math.max(0, Math.min(1, Math.asinh(10 * u) / ash));   // real M101 frames are already log-scaled
      if (fitsState.real === "m101") inferno(Math.pow(t, 0.85), px, o); else px[o] = px[o + 1] = px[o + 2] = 255 * t;
    }
    px[o + 3] = 255;
  }
  fx.putImageData(img, 0, 0);
  const tag = document.getElementById("fitsTag"), f = diff ? null : F[k];
  const ds = d => d ? d.toISOString().slice(0, 10) : "no date";
  tag.textContent = diff ? "Last − First" : `${n > 2 ? (k + 1) + "/" + n : k === 0 ? "First" : "Last"}${f.date || !fitsState.real ? " · " + ds(f.date) : ""}${f.wl ? " · " + f.wl.toFixed(3) + " µm" : ""}`;
  tag.style.color = diff ? "var(--ink)" : k === 0 ? "var(--past)" : k === n - 1 ? "var(--future)" : "var(--ink)";
  document.getElementById("diffLegend").hidden = !diff;
  const st = document.getElementById("simTag"); st.hidden = !fitsState.sim && !fitsState.real;
  st.textContent = fitsState.real === "demo" ? "Demo data (synthetic)" : fitsState.real ? "Real SPHEREx data" : "Simulated";
  st.classList.toggle("real", !!fitsState.real && fitsState.real !== "demo");
  const fr = document.getElementById("fitsFrame"); document.getElementById("frameRow").hidden = n < 3;
  fr.max = n - 1; if (!diff) fr.value = k;
  dirty = true; syncCutPop();
}
function renderMeta() {
  const F = fitsState.frames, el = document.getElementById("fitsMeta"); if (!F.length) { el.innerHTML = ""; return; }
  const A = F[0], Z = F[F.length - 1], dstr = d => d ? d.toISOString().slice(0, 16).replace("T", " ") + " UTC" : "no date in header";
  if (fitsState.real === "m101") {
    el.innerHTML = `<span class="k">Images</span><span>${F.length}, oldest first</span><span class="k">Band</span><span>2 (SPHEREx-D2)</span>
      <span class="k">Wavelength</span><span>1.113–1.618 µm</span><span class="k">Size</span><span>${A.w} × ${A.h} px, 0.1° (6′)</span>
      <span class="k">Data</span><span>Real SPHEREx cutouts (DR1/QR3), taken from the figure in IRSA's cutouts tutorial. Brightness is read back from its colour scale, and the figure gives no dates.</span>`;
    document.getElementById("showOnSky").disabled = false; return;
  }
  let html = F.length > 2 ? `<span class="k">Images</span><span>${F.length}</span>` : "";
  html += `<span class="k">${F.length > 2 ? "First" : "Then"}</span><span>${dstr(A.date)}</span>`;
  if (F.length > 1) html += `<span class="k">${F.length > 2 ? "Last" : "Now"}</span><span>${dstr(Z.date)}</span>`;
  if (F.length > 1 && A.date && Z.date) html += `<span class="k">Span</span><span>${Math.round(Math.abs(Z.date - A.date) / 864e5)} days</span>`;
  const wls = F.map(f => f.wl).filter(Number.isFinite);
  if (wls.length) html += `<span class="k">Wavelength</span><span>${Math.min(...wls).toFixed(2)}–${Math.max(...wls).toFixed(2)} µm</span>`;
  html += `<span class="k">Size</span><span>${A.w} × ${A.h} px${A.zodi ? " · zodiacal glow removed" : ""}</span>`;
  if (fitsState.sim) html += `<span class="k">Data</span><span>Simulated preview, not real SPHEREx images</span>`;
  else if (fitsState.real === "live") html += `<span class="k">Data</span><span>Real SPHEREx cutouts from NASA's IRSA archive (${LIVE.collection})</span>`;
  else if (fitsState.real === "demo") html += `<span class="k">Data</span><span>Synthetic images from the backend's demo mode</span>`;
  else if (F.length > 1) html += `<span class="k">Aligned</span><span>${F.every(f => f.wcs) ? "by sky coordinates" : "by pixel grid (no WCS found)"}</span>`;
  else if (!fitsState.sim) html += `<span class="k">Next</span><span>Add a second date to blink</span>`;
  el.innerHTML = html;
  document.getElementById("showOnSky").disabled = !(A.wcs || fitsState.center);
}
// ---------- enlarged player on the map (click the video square) ----------
// The square grows in place on the sky; a small control bar sits under it and moves with the map.
const cutPop = document.getElementById("cutPop");
let cutBig = false, cutGrow = 0;
function cutBigHalf() {      // half-size of the enlarged square, in px
  if (W > 820) return Math.min(210, 0.3 * H, 0.26 * (W - 392));
  return Math.max(70, Math.min(0.42 * W, (H - 150 - 24 - 165 - 170) / 2));      // phone: the panel is tucked away while enlarged
}
function cutSpot() {      // where the enlarged square goes on screen: in the open part of the map, leaving room for the bar
  const h = cutBigHalf();
  if (W > 820) return [(W - 392) / 2, Math.min(H / 2, H - 112 - 150 - h)];
  return [W / 2, 150 + h + 24];
}
function centreOn(v, sx, sy) {      // view centre that puts sky vector v at screen point (sx, sy), at the current zoom
  setBasis(); const ra = Math.atan2(v[1], v[0]), de = Math.asin(v[2]);
  const E = [-Math.sin(ra), Math.cos(ra), 0], N = [-Math.sin(de) * Math.cos(ra), -Math.sin(de) * Math.sin(ra), Math.cos(de)];
  const dx = (sx - W / 2) / scale, dy = (sy - H / 2) / scale, c = v.map((q, i) => q + E[i] * dx + N[i] * dy), n = Math.hypot(...c);
  return [(Math.atan2(c[1], c[0]) + 2 * Math.PI) % (2 * Math.PI), Math.asin(c[2] / n)];
}
function setFitsMode(m) { fitsState.mode = m; document.querySelectorAll(".modes .btn").forEach(x => x.classList.toggle("on", x.dataset.mode === m)); renderFits(); }
function syncCutPop() {
  if (!cutBig || !fitsState.aligned.length) return;
  const n = fitsState.aligned.length, diff = fitsState.mode === "diff";
  const sl = document.getElementById("cpFrame"); sl.max = n - 1; sl.hidden = n < 2; if (!diff) sl.value = viewIndex();
  document.getElementById("cpPlay").textContent = fitsState.mode === "play" ? "❚❚ Pause" : "▶ Play";
  document.getElementById("cpDiff").classList.toggle("on", diff);
  document.getElementById("cpLegend").hidden = !diff;
  document.getElementById("cpNote").textContent = diff ? "Only what changed between the first and last image. A moving star shows as an orange spot (where it was) next to a blue spot (where it went)."
    : `Image ${viewIndex() + 1} of ${n}. Watch for a dot that shifts position (a moving object) or changes brightness (a variable star).`;
}
function placeCutBar() {      // keep the control bar glued under the square
  const show = cutBig && cutGrow > 0.97 && cutRect;
  if (!show) { if (!cutPop.hidden) cutPop.hidden = true; return; }
  const w = Math.max(300, cutRect.half * 2 + 2), l = Math.max(8, Math.min(W - w - 8, cutRect.x - w / 2));
  cutPop.style.left = l + "px"; cutPop.style.top = (cutRect.y + cutRect.half) + "px"; cutPop.style.width = w + "px";
  if (cutPop.hidden) { cutPop.hidden = false; syncCutPop(); }
}
function openCutPop() {
  if (!fitsState.aligned.length) return;
  if (!cutoutMark && fitsState.center) { const [r, d] = fitsState.center; cutoutMark = radec(r * Math.PI / 180, d * Math.PI / 180); }
  if (!cutoutMark) return;
  cutBig = true; document.body.classList.add("cutbig"); const [sx, sy] = cutSpot(), [ra, de] = centreOn(cutoutMark, sx, sy);
  if (reduceMotion) { state.ra = ra; state.dec = de; } else { journey = journey || {}; journey.view = { t0: performance.now(), dur: 700, ra0: state.ra, dra: ((ra - state.ra + 3 * Math.PI) % (2 * Math.PI)) - Math.PI, dec0: state.dec, dec1: de, fov0: state.fov, fov1: state.fov }; }
  dirty = true;
}
function closeCutPop() { cutBig = false; document.body.classList.remove("cutbig"); cutPop.hidden = true; dirty = true; }
function stepCutGrow(dt) {
  const target = cutBig ? 1 : 0; if (cutGrow === target) return;
  cutGrow = reduceMotion ? target : target > cutGrow ? Math.min(1, cutGrow + dt * 3.5) : Math.max(0, cutGrow - dt * 4.5); dirty = true;
}
document.getElementById("cpClose").onclick = closeCutPop;
document.querySelector("#pane-spherex .viewer").addEventListener("click", openCutPop);      // the small player opens the big one too
window.addEventListener("keydown", e => { if (e.key === "Escape" && cutBig) closeCutPop(); });
document.getElementById("cpPlay").onclick = () => { if (fitsState.mode === "play") { fitsState.k = viewIndex(); setFitsMode("frame"); } else setFitsMode("play"); };
const stepFrame = d => { const n = fitsState.aligned.length; fitsState.k = ((fitsState.mode === "diff" ? fitsState.k : viewIndex()) + d + n) % n; setFitsMode("frame"); };
document.getElementById("cpPrev").onclick = () => stepFrame(-1);
document.getElementById("cpNext").onclick = () => stepFrame(1);
document.getElementById("cpDiff").onclick = () => setFitsMode(fitsState.mode === "diff" ? "play" : "diff");
document.getElementById("cpFrame").addEventListener("input", e => { fitsState.k = +e.target.value; setFitsMode("frame"); });
document.querySelectorAll(".modes .btn").forEach(b => b.onclick = () => {
  fitsState.mode = b.dataset.mode; document.querySelectorAll(".modes .btn").forEach(x => x.classList.toggle("on", x === b)); renderFits();
});
document.getElementById("fitsFrame").addEventListener("input", e => {
  fitsState.mode = "frame"; fitsState.k = +e.target.value; document.querySelectorAll(".modes .btn").forEach(x => x.classList.remove("on")); renderFits();
});
document.getElementById("contrast").addEventListener("input", restretch);
setInterval(() => { if (fitsState.mode === "play" && (state.tab === "pane-spherex" || cutoutMark || cutBig) && fitsState.aligned.length > 1) { fitsState.k = (fitsState.k + 1) % fitsState.aligned.length; renderFits(); } },
  550);
async function handleFiles(files) {
  const loaded = [];
  try { for (const f of [...files]) loaded.push(...loadFITS(await f.arrayBuffer(), f.name)); }
  catch (err) { document.getElementById("fitsMeta").innerHTML = `<span class="k">Error</span><span>${err.message}</span>`; return; }
  if (!loaded.length) return;
  let frames = fitsState.sim ? [] : fitsState.frames;
  if (loaded.length > 1) frames = [];                       // a new set replaces the old one
  frames = frames.concat(loaded).sort((a, b) => (a.date || 0) - (b.date || 0)).slice(0, 300);
  Object.assign(fitsState, { frames, sim: false, real: null, k: 0, label: null, center: null });
  document.querySelectorAll("#tlPresets .chip").forEach(c => c.classList.remove("on"));
  if (A0wcs(frames)) fitsState.fov = Math.min(frames[0].w, frames[0].h) * 6.15 / 3600;
  alignAll(); renderMeta();
  const A = frames[0];
  if (A.wcs) { const [ra, de] = A.wcs.pix2sky(A.w / 2 + 0.5, A.h / 2 + 0.5); cutoutMark = radec(ra, de); dirty = true; }
}
const drop = document.getElementById("drop");
drop.addEventListener("dragover", e => { e.preventDefault(); drop.classList.add("over"); });
drop.addEventListener("dragleave", () => drop.classList.remove("over"));
drop.addEventListener("drop", e => { e.preventDefault(); drop.classList.remove("over"); handleFiles(e.dataTransfer.files); });
document.getElementById("fitsInput").addEventListener("change", e => handleFiles(e.target.files));
document.getElementById("showTutorial").onclick = () => { setLayer("spherex", true); flyTo(TUTORIAL_RA * Math.PI / 180, TUTORIAL_DEC * Math.PI / 180, 25); };
document.getElementById("showCloud").onclick = () => { setLayer("spherex", true); flyTo(CLOUD_RA * Math.PI / 180, CLOUD_DEC * Math.PI / 180, 25); };
document.getElementById("showM101").onclick = () => { setLayer("spherex", true); flyTo(M101_RA * Math.PI / 180, M101_DEC * Math.PI / 180, 25); };
document.getElementById("showOnSky").onclick = () => {
  const A = fitsState.frames[0]; if (!A) return;
  if (fitsState.center) { const [r, d] = fitsState.center; cutoutMark = radec(r * Math.PI / 180, d * Math.PI / 180); flyTo(r * Math.PI / 180, d * Math.PI / 180, 8); return; }
  if (!A.wcs) return;
  const [ra, de] = A.wcs.pix2sky(A.w / 2 + 0.5, A.h / 2 + 0.5); cutoutMark = radec(ra, de); flyTo(ra, de, 6);
};
function A0wcs(frames) { return frames[0] && frames[0].wcs; }
// Real SPHEREx cutouts of M101 (band 2), recovered from the figure in IRSA's cutouts tutorial (BSD 3-Clause, Caltech).
const SPX = window.SPHEREX_DATA;
(function prepAtlas() {      // comet positions from the SPHEREx observation list (obscore_3I_Atlas.vot)
  const a = SPX.atlas, MON = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
  ATLAS = { n: a.mjd.length, mjd: a.mjd, det: a.det, v: a.ra.map((r, i) => radec(r * Math.PI / 180, a.dec[i] * Math.PI / 180)), days: [] };
  let lastDay = "";
  a.mjd.forEach((m, i) => { const d = new Date((m - 40587) * 864e5), key = d.toISOString().slice(0, 10);
    if (key !== lastDay) { ATLAS.days.push({ i, label: d.getUTCDate() + " " + MON[d.getUTCMonth()], key }); lastDay = key; } });
  ATLAS.days.forEach((d, k) => { const nx = ATLAS.days[k + 1], pv = ATLAS.days[k - 1];
    d.end = !pv || !nx || ATLAS.mjd[nx.i] - ATLAS.mjd[d.i] > 30 || ATLAS.mjd[d.i] - ATLAS.mjd[pv.i] > 30; });
})();
(function atlasPanel() {
  const rows = {}; ATLAS.mjd.forEach((m, i) => { const k = new Date((m - 40587) * 864e5).toISOString().slice(0, 10); (rows[k] = rows[k] || { n: 0, b: new Set() }).n++; rows[k].b.add(ATLAS.det[i]); });
  document.getElementById("atlasDays").innerHTML += Object.entries(rows).map(([k, r]) => `<tr><td>${new Date(k + "T00:00:00Z").toLocaleDateString("en-GB", { day: "numeric", month: "short", year: "numeric", timeZone: "UTC" })}</td><td>${r.n}</td><td>${[...r.b].sort().join(", ")}</td></tr>`).join("");
  const go = (i0, i1, fov) => { let x = 0, y = 0, z = 0; for (let i = i0; i <= i1; i++) { x += ATLAS.v[i][0]; y += ATLAS.v[i][1]; z += ATLAS.v[i][2]; }
    const n = Math.hypot(x, y, z); setLayer("atlas", true); if (Math.abs(state.year - NOW_YEAR) > 50) animateYear(NOW_YEAR, 600); flyTo(Math.atan2(y, x), Math.asin(z / n), fov); };
  const split = ATLAS.mjd.findIndex((m, i) => i > 0 && m - ATLAS.mjd[i - 1] > 30);
  document.getElementById("atlasAug").onclick = () => go(0, split - 1, 22);
  document.getElementById("atlasDec").onclick = () => go(split, ATLAS.n - 1, 18);
  document.getElementById("atlasAll").onclick = () => {      // centre between the two ends, nudged so the panel doesn't hide December
    const a = ATLAS.v[0], b = ATLAS.v[ATLAS.n - 1], x = a[0] + b[0], y = a[1] + b[1], z = a[2] + b[2], n = Math.hypot(x, y, z);
    setLayer("atlas", true); if (Math.abs(state.year - NOW_YEAR) > 50) animateYear(NOW_YEAR, 600);
    flyTo(Math.atan2(y, x) - 0.25, Math.asin(z / n), 140); };
  window.goAtlas = go; window.atlasSplit = split;
})();
// ---------- 3I/ATLAS journey animation ----------
// Real positions in time order. Time is squeezed: 1 day of real time = ATLAS_DAY seconds on screen,
// and the 4-month gap with no SPHEREx images is crossed in ATLAS_GAP seconds.
const ATLAS_DAY = 0.75, ATLAS_GAP = 3.2, BAND_UM = ["0.75–1.09", "1.10–1.62", "1.63–2.41", "2.42–3.82", "3.83–4.41", "4.42–5.00"];
let atlasAnim = null;
const ATLAS_T = (() => { const T = [0]; for (let i = 1; i < ATLAS.n; i++) { const dd = ATLAS.mjd[i] - ATLAS.mjd[i - 1];
  T.push(T[i - 1] + (dd > 30 ? ATLAS_GAP : Math.min(1.6, Math.max(0.015, dd * ATLAS_DAY)))); } return T; })();
function slerp(a, b, f) {
  const d = Math.max(-1, Math.min(1, a[0] * b[0] + a[1] * b[1] + a[2] * b[2])), w = Math.acos(d);
  if (w < 1e-9) return a.slice(); const sa = Math.sin((1 - f) * w) / Math.sin(w), sb = Math.sin(f * w) / Math.sin(w);
  return [a[0] * sa + b[0] * sb, a[1] * sa + b[1] * sb, a[2] * sa + b[2] * sb];
}
function sunVec(mjd) {      // low-precision Sun direction (equatorial), good to ~0.01°
  const n = mjd - 51544.5, g = (357.528 + 0.9856003 * n) * Math.PI / 180;
  const lam = (280.460 + 0.9856474 * n + 1.915 * Math.sin(g) + 0.020 * Math.sin(2 * g)) * Math.PI / 180, eps = 23.439 * Math.PI / 180;
  return [Math.cos(lam), Math.cos(eps) * Math.sin(lam), Math.sin(eps) * Math.sin(lam)];
}
const vecRaDec = v => [Math.atan2(v[1], v[0]), Math.asin(Math.max(-1, Math.min(1, v[2])))];
function startAtlasAnim() {
  stopAll(); setLayer("atlas", true); state.sel = -1; lastStarKey = "";
  if (Math.abs(state.year - NOW_YEAR) > 50) animateYear(NOW_YEAR, 0);
  if (reduceMotion) { document.getElementById("atlasAll").click(); return; }
  const [ra, dec] = vecRaDec(ATLAS.v[0]); flyTo(ra, dec, 18);
  atlasAnim = { t: 0, wait: 1.5, i: 0, head: ATLAS.v[0], mjd: ATLAS.mjd[0], started: false, done: false, paused: false, follow: true };
  document.getElementById("atlasCap").hidden = false; document.getElementById("acPause").textContent = "Pause"; updateAtlasCap();
}
function endAtlasAnim() { atlasAnim = null; const c = document.getElementById("atlasCap"); if (c) c.hidden = true; dirty = true; }
function stepAtlasAnim(now, dt) {
  const a = atlasAnim, A = ATLAS, T = ATLAS_T, end = T[A.n - 1];
  if (a.paused || a.done) return;
  if (a.wait > 0) { a.wait -= dt; if (a.wait > 0) return; }
  a.started = true; a.t = Math.min(end, a.t + dt);
  let lo = 0, hi = A.n - 1; while (hi - lo > 1) { const m = (lo + hi) >> 1; if (T[m] <= a.t) lo = m; else hi = m; }
  const f = T[hi] > T[lo] ? (a.t - T[lo]) / (T[hi] - T[lo]) : 1;
  a.i = a.t >= end ? A.n - 1 : lo; a.f = f; a.gap = A.mjd[hi] - A.mjd[lo] > 30 && a.t < end;
  a.head = a.t >= end ? A.v[A.n - 1] : slerp(A.v[lo], A.v[hi], f);
  a.mjd = a.t >= end ? A.mjd[A.n - 1] : A.mjd[lo] + (A.mjd[hi] - A.mjd[lo]) * f;
  if (a.follow && !(journey && journey.view)) {      // camera glides after the comet; zooms out to cross the gap
    const [ra, dec] = vecRaDec(a.head), k = 1 - Math.exp(-dt * (a.gap ? 6 : 3));
    const dra = ((ra - state.ra + 3 * Math.PI) % (2 * Math.PI)) - Math.PI;
    state.ra = (state.ra + dra * k + 2 * Math.PI) % (2 * Math.PI); state.dec += (dec - state.dec) * k;
    const fovT = a.gap ? 18 - 2 * f + 75 * Math.sin(Math.PI * f) : (a.i >= atlasSplit ? 16 : 18);
    state.fov += (fovT - state.fov) * k;
  }
  if (a.t >= end) { a.done = true; document.getElementById("acPause").textContent = "Done"; }
  dirty = true; updateAtlasCap();
}
function drawComet(a) {
  const h = a.head; if (!project(h[0], h[1], h[2])) return;
  const x = PR[0], y = PR[1], s = sunVec(a.mjd), sd = s[0] * h[0] + s[1] * h[1] + s[2] * h[2];
  let tx = -(s[0] - sd * h[0]), ty = -(s[1] - sd * h[1]), tz = -(s[2] - sd * h[2]), tn = Math.hypot(tx, ty, tz);
  if (tn > 1e-6) {      // tail points away from the Sun, like a real comet's
    tx /= tn; ty /= tn; tz /= tn; const e = 0.004, p2 = [h[0] + tx * e, h[1] + ty * e, h[2] + tz * e];
    if (project(p2[0], p2[1], p2[2])) {
      let dx = PR[0] - x, dy = PR[1] - y; const dn = Math.hypot(dx, dy) || 1; dx /= dn; dy /= dn;
      const L = 58, w = 5.5, g = ctx.createLinearGradient(x, y, x + dx * L, y + dy * L);
      g.addColorStop(0, "rgba(200,235,250,0.55)"); g.addColorStop(1, "rgba(200,235,250,0)");
      ctx.fillStyle = g; ctx.beginPath(); ctx.moveTo(x - dy * w * 0.6, y + dx * w * 0.6);
      ctx.lineTo(x + dx * L - dy * w, y + dy * L + dx * w); ctx.lineTo(x + dx * L + dy * w, y + dy * L - dx * w);
      ctx.lineTo(x + dy * w * 0.6, y - dx * w * 0.6); ctx.closePath(); ctx.fill();
    }
  }
  const pulse = 1 + 0.12 * Math.sin(performance.now() / 180), g2 = ctx.createRadialGradient(x, y, 0, x, y, 13 * pulse);
  g2.addColorStop(0, "rgba(255,255,255,1)"); g2.addColorStop(0.25, "rgba(190,230,250,0.85)"); g2.addColorStop(1, "rgba(130,205,230,0)");
  ctx.fillStyle = g2; ctx.beginPath(); ctx.arc(x, y, 13 * pulse, 0, 7); ctx.fill();
  ctx.font = "italic 600 16px 'Cormorant Garamond', Georgia, serif"; ctx.fillStyle = "rgba(130,205,230,0.95)"; ctx.fillText("3I/ATLAS", x + 14, y - 12);
}
function updateAtlasCap() {
  const a = atlasAnim; if (!a) return; const A = ATLAS, d = new Date((a.mjd - 40587) * 864e5); document.getElementById("acTitle").textContent = "3I/ATLAS";
  document.getElementById("acDate").textContent = d.toLocaleString("en-GB", { day: "numeric", month: "short", year: "numeric", hour: "2-digit", minute: "2-digit", timeZone: "UTC" }) + " UTC";
  document.getElementById("acFill").style.width = (100 * a.t / ATLAS_T[A.n - 1]).toFixed(1) + "%";
  const info = document.getElementById("acInfo");
  if (!a.started) info.textContent = "Flying to where SPHEREx first saw it, 1 August 2025…";
  else if (a.done) info.innerHTML = `Journey complete: ${A.n} SPHEREx images, August and December 2025. The comet is now heading back out of the solar system for good.`;
  else if (a.gap) info.innerHTML = "SPHEREx took no pictures for four months (Sep to Nov), so we skip ahead along the dotted line. The real path in between isn't in this data.";
  else { const det = A.det[a.i]; info.innerHTML = `Image ${a.i + 1} of ${A.n} · <span class="sw" style="background:${BAND_COLORS[det - 1]}"></span>detector ${det} (${BAND_UM[det - 1]} µm). Each dot is one picture; the line is the comet's real track.`; }
}
document.getElementById("atlasPlay").onclick = startAtlasAnim;
document.getElementById("acReplay").onclick = () => m101Anim ? startM101Anim() : startAtlasAnim();
document.getElementById("acClose").onclick = () => { endAtlasAnim(); endM101Anim(); };
document.getElementById("acPause").onclick = () => { const a = atlasAnim || m101Anim; if (!a || a.done) return;
  a.paused = !a.paused; if (!a.paused) a.follow = true; document.getElementById("acPause").textContent = a.paused ? "Resume" : "Pause"; };
// ---------- Pinwheel Galaxy (M101) zoom animation ----------
// Wide view at the Big Dipper → zoom in → the real SPHEREx mosaic of the whole galaxy fades in at its true size,
// then crossfades to the 3.3 µm dust map and back. Images: IRSA SPHEREx mosaic tutorial (made by tools/extract_m101_mosaic.py).
const M101_V = radec(M101_RA * Math.PI / 180, M101_DEC * Math.PI / 180), M101_GALAXY = 28.8 / 60;   // degrees
const M101_T = { fly: 1.5, wide: 4, zoom: 6.5, fade: 1.2, star1: 5, cross: 1.4, dust: 5.5, star2: 4 }, M101_FOV0 = 50;
const M101_IMG = {}; for (const k of ["mosaic", "pah"]) { M101_IMG[k] = new Image(); M101_IMG[k].src = SPX.m101disc[k]; }
let m101Anim = null;
const m101Total = () => Object.values(M101_T).reduce((s, x) => s + x, 0);
const m101Wide = () => W > 820;
function m101Spot() {      // where on screen the galaxy should sit, and how big (radius in px), given the panel and caption
  if (m101Wide()) return { dx: -(196 + 0.12 * W), dy: 0, r: Math.min(290, 0.34 * H, 0.3 * (W - 392)) };
  const capB = document.getElementById("atlasCap").getBoundingClientRect().bottom || 160, panT = document.getElementById("panel").getBoundingClientRect().top || H;
  const top = capB + 40, bot = Math.max(top + 120, panT - 10);
  return { dx: 0, dy: (top + bot) / 2 - H / 2, r: Math.min(0.42 * W, (bot - top) / 2) };
}
function m101Fov1() { const ang = SPX.m101disc.deg / 2 * Math.PI / 180; return Math.max(0.3, Math.max(W, H) * ang / m101Spot().r * 180 / Math.PI); }
function m101Center(fov) {      // view centre that puts M101 at m101Spot() on screen
  const sp = m101Spot(), sc = (Math.max(W, H) / 2) / (2 * Math.tan(fov * Math.PI / 720)), ra = Math.atan2(M101_V[1], M101_V[0]), de = M101_DEC * Math.PI / 180;
  const E = [-Math.sin(ra), Math.cos(ra), 0], N = [-Math.sin(de) * Math.cos(ra), -Math.sin(de) * Math.sin(ra), Math.cos(de)];      // the map shows east to the left
  const c = M101_V.map((v, i) => v + E[i] * sp.dx / sc + N[i] * sp.dy / sc), n = Math.hypot(...c);
  return [(Math.atan2(c[1], c[0]) + 2 * Math.PI) % (2 * Math.PI), Math.asin(c[2] / n)];
}
function m101View(fov) { state.fov = fov; [state.ra, state.dec] = m101Center(fov); dirty = true; }
function m101Closeup(on) { document.body.classList.toggle("closeup", on); document.getElementById("atlasCap").classList.toggle("right", on && m101Wide()); }
function startM101Anim() {
  stopAll(); setLayer("spherex", true); state.sel = -1; lastStarKey = "";
  if (Math.abs(state.year - NOW_YEAR) > 50) animateYear(NOW_YEAR, 0);
  document.getElementById("atlasCap").hidden = false; document.getElementById("acPause").textContent = "Pause";
  if (reduceMotion) { m101View(m101Fov1()); m101Anim = { t: m101Total(), img: 1, dust: 0, done: true, follow: false, paused: false }; m101Closeup(true); updateM101Cap(); return; }
  const [ra0, de0] = m101Center(M101_FOV0); flyTo(ra0, de0, M101_FOV0);
  m101Anim = { t: 0, img: 0, dust: 0, done: false, paused: false, follow: true }; m101Closeup(false); updateM101Cap();
}
function endM101Anim() { if (!m101Anim) return; m101Anim = null; m101Closeup(false); document.getElementById("atlasCap").hidden = true; dirty = true; }
function stepM101Anim(dt) {
  const a = m101Anim, T = M101_T; if (a.paused || a.done) { dirty = true; return; }
  a.t = Math.min(m101Total(), a.t + dt);
  const tz = a.t - T.fly - T.wide, ti = tz - T.zoom, fov1 = m101Fov1();
  if (a.follow && tz > 0 && !(journey && journey.view)) m101View(Math.exp(Math.log(M101_FOV0) + (Math.log(fov1) - Math.log(M101_FOV0)) * ease(Math.min(1, tz / T.zoom))));
  a.img = ti <= 0 ? 0 : Math.min(1, ti / T.fade);
  const td = ti - T.fade - T.star1;      // dust crossfade in, hold, crossfade out
  a.dust = td <= 0 ? 0 : td < T.cross ? ease(td / T.cross) : td < T.cross + T.dust ? 1 : Math.max(0, 1 - ease((td - T.cross - T.dust) / T.cross));
  m101Closeup(ti > -1.5);
  if (a.t >= m101Total()) { a.done = true; a.dust = 0; document.getElementById("acPause").textContent = "Done"; }
  dirty = true; updateM101Cap();
}
function drawM101(a) {
  if (!project(M101_V[0], M101_V[1], M101_V[2])) return;
  const x = PR[0], y = PR[1], ra = Math.atan2(M101_V[1], M101_V[0]), de = M101_DEC * Math.PI / 180;
  const E = [-Math.sin(ra), Math.cos(ra), 0], N = [-Math.sin(de) * Math.cos(ra), -Math.sin(de) * Math.sin(ra), Math.cos(de)];
  const at = (dE, dN) => { const r = M101_V.map((c, i) => c + E[i] * dE + N[i] * dN); project(r[0], r[1], r[2]); return [PR[0], PR[1]]; };
  const rad = d => d * Math.PI / 180;
  const gR = Math.hypot(...at(0, rad(M101_GALAXY / 2)).map((v, i) => v - [x, y][i]));
  if (gR > 14 && a.img < 1) {      // dashed outline of the whole galaxy (about 29′), until the real image fills it
    ctx.globalAlpha = 1 - a.img; ctx.strokeStyle = "rgba(124,200,255,0.6)"; ctx.setLineDash([4, 6]); ctx.lineWidth = 1; ctx.beginPath(); ctx.arc(x, y, gR, 0, 7); ctx.stroke(); ctx.setLineDash([]);
    if (gR > 60) { ctx.font = "11.5px 'IBM Plex Sans', sans-serif"; ctx.fillStyle = "rgba(124,200,255,0.85)"; ctx.fillText("Whole galaxy, about 170,000 light-years across", x + gR * 0.42, y - gR * 0.9 - 6); }
    ctx.globalAlpha = 1;
  } else if (gR <= 14) {      // wide view: pulsing ring to point it out
    const r = 10 + 3 * Math.sin(performance.now() / 260); ctx.strokeStyle = "rgba(124,200,255,0.9)"; ctx.lineWidth = 1.4; ctx.beginPath(); ctx.arc(x, y, r, 0, 7); ctx.stroke();
  }
  ctx.font = "italic 600 16px 'Cormorant Garamond', Georgia, serif"; ctx.fillStyle = "#7cc8ff";
  if (a.img > 0) {      // the real mosaic, north up / east left, at its true size on the sky
    const h = rad(SPX.m101disc.deg / 2), [wx, wy] = at(-h, 0), [nx, ny] = at(0, h), I = M101_IMG.mosaic, w = I.naturalWidth || 287, hh = I.naturalHeight || 287;
    const A = 2 * (wx - x) / w, B = 2 * (wy - y) / w, C = 2 * (x - nx) / hh, D = 2 * (y - ny) / hh, R = Math.hypot(wx - x, wy - y);
    const glow = ctx.createRadialGradient(x, y, R * 0.3, x, y, R * 1.25);      // soft warm halo so it stands out
    glow.addColorStop(0, `rgba(250,140,40,${0.22 * a.img})`); glow.addColorStop(1, "rgba(250,140,40,0)");
    ctx.fillStyle = glow; ctx.beginPath(); ctx.arc(x, y, R * 1.25, 0, 7); ctx.fill();
    ctx.save(); ctx.transform(A, B, C, D, x - (A * w + C * hh) / 2, y - (B * w + D * hh) / 2); ctx.imageSmoothingEnabled = true; ctx.imageSmoothingQuality = "high";
    ctx.globalAlpha = a.img * (1 - a.dust); ctx.drawImage(M101_IMG.mosaic, 0, 0);
    if (a.dust > 0) { ctx.globalAlpha = a.img * a.dust; ctx.drawImage(M101_IMG.pah, 0, 0); }
    ctx.restore(); ctx.globalAlpha = 1;
    ctx.textAlign = "center"; ctx.font = "italic 600 18px 'Cormorant Garamond', Georgia, serif"; ctx.fillStyle = "#7cc8ff"; ctx.fillText("Pinwheel Galaxy (M101)", x, y - R - 24);
    ctx.font = "11.5px 'IBM Plex Mono', monospace"; ctx.fillStyle = "rgba(235,228,211,0.85)";
    ctx.fillText(a.dust > 0.5 ? "SPHEREx · 3.3 µm glowing dust" : "SPHEREx · 102 colours, 0.75–5 µm", x, y - R - 8); ctx.textAlign = "start";
  } else ctx.fillText("Pinwheel Galaxy (M101)", x - 150, y - 16);
}
function updateM101Cap() {
  const a = m101Anim; if (!a) return; const T = M101_T;
  document.getElementById("acTitle").textContent = "Pinwheel Galaxy (M101)";
  document.getElementById("acFill").style.width = (100 * a.t / m101Total()).toFixed(1) + "%";
  const tz = a.t - T.fly - T.wide, info = document.getElementById("acInfo"), date = document.getElementById("acDate");
  if (a.t < T.fly) { date.textContent = "Ursa Major"; info.textContent = "Flying to the Big Dipper…"; }
  else if (tz <= 0) { date.textContent = "21 million ly away"; info.textContent = "The galaxy sits just above the Big Dipper's handle. Its light left it 21 million years ago, long before the first humans."; }
  else if (a.img < 1 && !a.done) { date.textContent = "zoom ×" + Math.round(M101_FOV0 / state.fov); info.textContent = "Zooming in. To your eyes it's only a faint smudge; a telescope shows a spiral nearly twice as wide as our Milky Way."; }
  else if (a.done) { date.textContent = "SPHEREx mosaic"; info.innerHTML = "Both maps come from IRSA's SPHEREx mosaic tutorial. Supernovae were seen exploding in this galaxy in 2011 and 2023."; }
  else if (a.dust > 0.5) { date.textContent = "3.3 µm dust"; info.innerHTML = "The same galaxy at 3.3 µm with the starlight taken away. Glowing dust (PAH molecules) lights up the spiral arms, where new stars are being born."; }
  else { date.textContent = "0.75–5 µm"; info.innerHTML = "Real SPHEREx mosaic: all 102 infrared colours added together, at its true size on the sky. This is mostly starlight from billions of stars."; }
}
document.getElementById("m101Play").onclick = startM101Anim;
function loadRealM101(fly) {
  const frames = SPX.m101.map(f => { const u = bytes(f.px), d = new Float32Array(u.length); for (let i = 0; i < u.length; i++) d[i] = u[i];
    return { w: f.w, h: f.h, data: d, wcs: null, date: null, wl: f.wl, zodi: false }; });
  Object.assign(fitsState, { frames, sim: false, real: "m101", k: 0, mode: "play", fov: 0.1, center: [M101_RA, M101_DEC], label: "Pinwheel Galaxy (M101) · real SPHEREx images" });
  document.querySelectorAll(".modes .btn").forEach(x => x.classList.toggle("on", x.dataset.mode === "play"));
  document.querySelectorAll("#tlPresets .chip").forEach(c => c.classList.toggle("on", c.id === "chipM101"));
  if (fly) {
    document.getElementById("tlCoords").value = `${fmtRA(M101_RA)} ${fmtDec(M101_DEC)}`;
    document.getElementById("tlMsg").innerHTML = "<b>M101, the Pinwheel Galaxy</b> · 8 real band 2 images, in the order SPHEREx took them";
    cutoutMark = radec(M101_RA * Math.PI / 180, M101_DEC * Math.PI / 180); dirty = true;
    flyTo(M101_RA * Math.PI / 180, M101_DEC * Math.PI / 180, 8);
    document.querySelector("#pane-spherex .viewer").scrollIntoView({ block: "center", behavior: reduceMotion ? "auto" : "smooth" });
  }
  alignAll(); renderMeta();
}
document.querySelectorAll("img[data-fig]").forEach(img => img.src = SPX.figs[img.dataset.fig]);

// ---------- SPHEREx time-lapse: coordinates in, image series out ----------
const TL_PRESETS = [
  { name: "Barnard's Star", ra: 269.452077, dec: 4.693365, pmra: -801.551, pmdec: 10362.394, sed: [60, 70, 62, 42, 30, 26] },
  { name: "Epsilon Indi B", ra: 331.043833, dec: -56.782708, pmra: 3966.7, pmdec: -2536.2, sed: [6, 14, 10, 3, 11, 16] },
  { name: "WISE 0855−0714", ra: 133.7951321, dec: -7.24515, pmra: -8118.9, pmdec: 679.3, sed: [0, 0, 0.2, 0.8, 3.5, 6] },
  { name: "RW Aurigae", ra: 76.956068, dec: 30.401326, pmra: 4.9, pmdec: -24.6, sed: [28, 34, 40, 46, 42, 38], vary: v => [1, 0.5, 0.85][v] },
  { name: "R Draconis", ra: 248.16768, dec: 66.754984, pmra: -20.07, pmdec: 7.06, sed: [70, 115, 105, 68, 48, 40], vary: (v, mjd) => 1 + 0.6 * Math.sin(2 * Math.PI * (mjd - 60800) / 245) }];
const TL_BANDS = [[0.75, 1.09], [1.10, 1.62], [1.63, 2.41], [2.42, 3.82], [3.83, 4.41], [4.42, 5.00]];
function parseCoords(txt) {
  const t = txt.trim();
  if (!/^[-+\d\s.:hdms°'′"″,]+$/i.test(t) || !/\d/.test(t)) return null;
  const toks = t.replace(/[hdms°'′"″:,]/gi, " ").trim().split(/\s+/);
  const sexa = /[hms:°'′"″]/i.test(t) || toks.length === 6;
  if (!sexa && toks.length === 2) { const ra = +toks[0], dec = +toks[1]; return isFinite(ra) && isFinite(dec) && Math.abs(dec) <= 90 && ra >= 0 && ra < 360 ? { ra, dec } : null; }
  if (toks.length === 6) {
    const [h, m, sec, dd, dm, ds] = toks.map(Number); if ([h, m, sec, dd, dm, ds].some(v => !isFinite(v)) || h >= 24) return null;
    const dec = (toks[3].startsWith("-") ? -1 : 1) * (Math.abs(dd) + dm / 60 + ds / 3600);
    return Math.abs(dec) <= 90 ? { ra: (h + m / 60 + sec / 3600) * 15, dec } : null;
  }
  return null;
}
const fmtRA = ra => { const h = ra / 15, H = Math.floor(h), M = Math.floor((h - H) * 60), S_ = ((h - H) * 60 - M) * 60; return `${String(H).padStart(2, "0")}h${String(M).padStart(2, "0")}m${S_.toFixed(1).padStart(4, "0")}s`; };
const fmtDec = d => { const a = Math.abs(d), D_ = Math.floor(a), M = Math.floor((a - D_) * 60), S_ = ((a - D_) * 60 - M) * 60; return `${d < 0 ? "-" : "+"}${String(D_).padStart(2, "0")}d${String(M).padStart(2, "0")}m${S_.toFixed(0).padStart(2, "0")}s`; };
function simulateSeries(t, fov) {
  const N = Math.max(24, Math.min(120, Math.round(fov * 3600 / 6.15))), sc = fov / N, c = (N - 1) / 2;
  let seed = Math.round(t.ra * 1000 + t.dec * 7919) >>> 0 || 1; const rnd = () => (seed = (seed * 1664525 + 1013904223) >>> 0) / 4294967296;
  const gauss = () => Math.sqrt(-2 * Math.log(rnd() + 1e-12)) * Math.cos(2 * Math.PI * rnd());
  const cen = posOf(t, 2025.9), d2 = Math.PI / 180;
  const toGrid = (ra, dec) => {      // north up, east left
    const a = ra * d2, b = dec * d2, a0 = cen.ra * d2, b0 = cen.dec * d2, cc = Math.sin(b0) * Math.sin(b) + Math.cos(b0) * Math.cos(b) * Math.cos(a - a0);
    const xi = Math.cos(b) * Math.sin(a - a0) / cc / d2, eta = (Math.cos(b0) * Math.sin(b) - Math.sin(b0) * Math.cos(b) * Math.cos(a - a0)) / cc / d2;
    return [c - xi / sc, c + eta / sc]; };
  const stars = []; for (let k = 0; k < Math.round(70 * (fov / 0.1) ** 2); k++) stars.push([rnd() * N, rnd() * N, Math.pow(rnd(), 2.5) * 45 + 0.8, 0.5 + rnd()]);
  const frames = [];
  [60820, 61005, 61190].forEach((mjd0, v) => TL_BANDS.forEach(([lo, hi], b) => {
    const mjd = mjd0 + b * 1.4, year = 2000 + (mjd - 51544.5) / 365.25, data = new Float32Array(N * N);
    for (let i = 0; i < data.length; i++) data[i] = gauss() * 0.35;
    const add = (x, y, f) => { for (let j = Math.max(0, Math.floor(y - 4)); j <= Math.min(N - 1, y + 4); j++) for (let i = Math.max(0, Math.floor(x - 4)); i <= Math.min(N - 1, x + 4); i++)
      data[j * N + i] += f * Math.exp(-((i - x) ** 2 + (j - y) ** 2) / (2 * 0.75 ** 2)); };
    for (const [x, y, f, slope] of stars) add(x, y, f * Math.pow(0.85 + 0.05 * (b + 1), slope * 2));
    if (t.sed) { const p = posOf(t, year), [gx, gy] = toGrid(p.ra, p.dec); add(gx, gy, t.sed[b] * (t.vary ? t.vary(v, mjd) : 1)); }
    if (v === 1 && (b === 2 || b === 3)) add(N * 0.25 + b * 2.5, N * 0.7 - b * 1.2, 14);       // a passing asteroid
    frames.push({ w: N, h: N, data, wcs: null, date: new Date((mjd - 40587) * 864e5), wl: lo + (hi - lo) * (0.3 + 0.4 * rnd()), zodi: false, sim: true });
  }));
  return { frames, center: cen };
}
function posOf(t, year) {
  const dy = year - 2000, dec = t.dec + (t.pmdec || 0) * dy / 3.6e6;
  return { ra: t.ra + (t.pmra || 0) * dy / 3.6e6 / Math.cos(t.dec * Math.PI / 180), dec };
}
async function runTimelapse(t) {
  await backendReady;                      // a click before the backend check finishes waits for it
  if (BACKEND) return runLive(t);
  const fov = +document.getElementById("tlFov").value, { frames, center } = simulateSeries(t, fov);
  Object.assign(fitsState, { frames, sim: true, real: null, k: 0, mode: "play", fov, center: [center.ra, center.dec], label: (t.name || "Your patch") + " · simulated" });
  document.querySelectorAll(".modes .btn").forEach(x => x.classList.toggle("on", x.dataset.mode === "play"));
  document.querySelectorAll("#tlPresets .chip").forEach(c => c.classList.toggle("on", c.textContent === t.name));
  document.getElementById("tlCoords").value = `${fmtRA(center.ra)} ${fmtDec(center.dec)}`;
  const pm = Math.hypot(t.pmra || 0, t.pmdec || 0) / 1000;
  document.getElementById("tlMsg").innerHTML = `${t.name ? "<b>" + t.name + "</b> · " : ""}position for late 2025${pm > 0.05 ? ` · moves ${pm.toFixed(2)}″ a year` : ""}`;
  alignAll(); renderMeta();
  cutoutMark = radec(center.ra * Math.PI / 180, center.dec * Math.PI / 180); dirty = true;
  flyTo(center.ra * Math.PI / 180, center.dec * Math.PI / 180, 8);
  document.querySelector("#pane-spherex .viewer").scrollIntoView({ block: "center", behavior: reduceMotion ? "auto" : "smooth" });
}
document.getElementById("tlPresets").innerHTML = `<button class="chip real" id="chipM101">Pinwheel Galaxy · real data</button>` + TL_PRESETS.map(p => `<button class="chip">${p.name}</button>`).join("");
document.querySelectorAll("#tlPresets .chip:not(.real)").forEach((c, k) => c.onclick = () => runTimelapse(TL_PRESETS[k]));
document.getElementById("chipM101").onclick = () => loadRealM101(true);
const tlGo = () => {
  const txt = document.getElementById("tlCoords").value, preset = TL_PRESETS.find(p => p.name.toLowerCase() === txt.trim().toLowerCase());
  if (preset) return runTimelapse(preset);
  const c = parseCoords(txt);
  if (!c) { document.getElementById("tlMsg").innerHTML = `Couldn't read that. Try <code>08h55m10.8s -07d14m42.5s</code> or <code>133.795 -7.245</code>.`; return; }
  runTimelapse({ name: null, ra: c.ra, dec: c.dec });
};
document.getElementById("tlGo").onclick = tlGo;
document.getElementById("tlCoords").addEventListener("keydown", e => { if (e.key === "Enter") tlGo(); });
document.getElementById("tlFov").addEventListener("change", () => { if ((fitsState.sim || fitsState.real === "live" || fitsState.real === "demo") && fitsState.center) tlGo(); });

// ---------- live SPHEREx images through the backend (backend/app.py) ----------
// When the page is served by the backend, the coordinate box and test targets download real cutouts
// from NASA's IRSA archive instead of simulating. Change these to download more or fewer images.
const LIVE = { collection: "spherex_dr1", maxImages: 24, parallel: 4 };
let BACKEND = false, liveRun = 0;       // BACKEND becomes {demo: true|false} when /api/ping answers
async function checkBackend() {
  try {
    const r = await fetch("api/ping"), j = r.ok && (r.headers.get("content-type") || "").includes("json") ? await r.json() : null;
    BACKEND = j && j.ok ? { demo: !!j.demo } : false;
  } catch { BACKEND = false; }
  if (!BACKEND) return;
  document.querySelector(".tlwarn").innerHTML = BACKEND.demo
    ? "<b>Backend in demo mode.</b> Coordinates and test targets download <b>synthetic</b> images. Start the backend without --demo for real SPHEREx data."
    : "<b>Connected to the backend.</b> Coordinates and test targets download <b>real SPHEREx images</b> from NASA's IRSA archive. <b>Pinwheel Galaxy · real data</b> still plays the tutorial frames.";
  document.getElementById("tlFull").classList.add("ready");
}
async function runLive(t) {
  const fov = +document.getElementById("tlFov").value, c = posOf(t, 2025.9), run = ++liveRun, msg = document.getElementById("tlMsg");
  const say = h => { if (run === liveRun) msg.innerHTML = h; };
  const name = t.name || "Your patch", kind = BACKEND.demo ? "synthetic demo" : "real SPHEREx";
  document.querySelectorAll("#tlPresets .chip").forEach(ch => ch.classList.toggle("on", ch.textContent === t.name));
  document.getElementById("tlCoords").value = `${fmtRA(c.ra)} ${fmtDec(c.dec)}`;
  // clear the old pictures straight away, so nothing else shows under this name while downloading
  Object.assign(fitsState, { frames: [], aligned: [], stretch: [], fov, sim: false, real: null, k: 0, center: [c.ra, c.dec], label: `${name} · downloading…` });
  fc.width = fc.height = 96; fx.fillStyle = "#000"; fx.fillRect(0, 0, 96, 96);
  document.getElementById("fitsTag").textContent = "Downloading…"; document.getElementById("fitsTag").style.color = "var(--muted)";
  document.getElementById("simTag").hidden = true; document.getElementById("frameRow").hidden = true; renderMeta();
  const failName = why => { if (run === liveRun) { fitsState.label = `${name} · no pictures (${why})`; document.getElementById("fitsTag").textContent = "No pictures"; dirty = true; } };
  cutoutMark = radec(c.ra * Math.PI / 180, c.dec * Math.PI / 180); dirty = true;
  flyTo(c.ra * Math.PI / 180, c.dec * Math.PI / 180, 8);
  say("Searching NASA's IRSA archive for SPHEREx images of this spot…");
  try {
    const sr = await fetch(`api/search?ra=${c.ra.toFixed(6)}&dec=${c.dec.toFixed(6)}&collection=${LIVE.collection}`), res = await sr.json();
    if (!sr.ok) throw new Error(res.error || "search failed");
    let rows = res.rows.filter(r => r.access_url && r.t_min);
    if (!rows.length) { failName("none in the archive yet"); return say(`No SPHEREx images cover this spot in ${LIVE.collection} yet.`); }
    if (rows.length > LIVE.maxImages) rows = Array.from({ length: LIVE.maxImages }, (_, k) => rows[Math.round(k * (rows.length - 1) / (LIVE.maxImages - 1))]);
    const frames = [], size = Math.min(0.5, fov * 1.45 + 0.01);
    let next = 0, done = 0, failed = 0, lastErr = "", shownAt = 0;
    const show = () => {      // put what has arrived so far on screen (oldest first); called as pictures come in
      if (run !== liveRun || !frames.length) return;
      const first = !fitsState.frames.length, fr = frames.slice().sort((a, b) => a.date - b.date);
      Object.assign(fitsState, { frames: fr, sim: false, real: BACKEND.demo ? "demo" : "live", fov, center: [c.ra, c.dec], label: `${name} · ${kind}` });
      if (first) { fitsState.k = 0; fitsState.mode = "play"; document.querySelectorAll(".modes .btn").forEach(x => x.classList.toggle("on", x.dataset.mode === "play")); }
      alignAll(); renderMeta(); dirty = true; shownAt = performance.now();
    };
    async function worker() {
      while (next < rows.length && run === liveRun) {
        const r = rows[next++];
        try {
          const resp = await fetch(`api/cutout?url=${encodeURIComponent(r.access_url)}&ra=${c.ra.toFixed(7)}&dec=${c.dec.toFixed(7)}&size=${size.toFixed(4)}`);
          if (!resp.ok) throw new Error((await resp.json()).error || "cutout failed");
          for (const f of loadFITS(await resp.arrayBuffer(), r.obs_id)) {
            f.date = new Date(((r.t_max ? (r.t_min + r.t_max) / 2 : r.t_min) - 40587) * 864e5);
            if (r.em_min && r.em_max) f.wl = (r.em_min + r.em_max) / 2 * 1e6;
            frames.push(f);
          }
          if (frames.length === 1 || performance.now() - shownAt > 1500) show();
        } catch (e) { failed++; lastErr = e.message; }
        done++;
        say(`Downloading ${kind} cutouts… ${done} of ${rows.length}${failed ? ` · ${failed} failed` : ""}`);
      }
    }
    await Promise.all(Array.from({ length: LIVE.parallel }, worker));
    if (run !== liveRun) return;
    if (!frames.length) { failName("download failed"); return say(`Couldn't download any cutouts. ${lastErr}`); }
    show(); frames.sort((a, b) => a.date - b.date);
    const d = x => x.toISOString().slice(0, 10);
    say(`<b>${name}</b> · ${frames.length} ${kind} images, ${d(frames[0].date)} to ${d(frames.at(-1).date)}${failed ? ` (${failed} failed)` : ""}`);
    document.querySelector("#pane-spherex .viewer").scrollIntoView({ block: "center", behavior: reduceMotion ? "auto" : "smooth" });
  } catch (e) { failName("archive error"); say(`Couldn't reach the archive: ${e.message}`); }
}
const backendReady = checkBackend();
loadRealM101(false);
document.getElementById("chipM101").classList.add("on");

// ---------- boot ----------
function start(saved) {
  if (saved && typeof saved.year === "number") {
    Object.assign(state, { year: saved.year, ra: saved.ra, dec: saved.dec, fov: saved.fov, sel: saved.sel ?? -1 });
    if (saved.layers) Object.assign(state.layers, saved.layers);
    Object.entries(state.layers).forEach(([k, v]) => setLayer(k, v));
    if (saved.tab) showTab(saved.tab === "pane-names" ? "pane-star" : saved.tab);      // the Names tab is now part of Stars
  } else {
    state.sel = -1;      // start on the star list
    const q = new URLSearchParams(location.search), ra = parseFloat(q.get("ra")), dec = parseFloat(q.get("dec"));
    if (isFinite(ra) && isFinite(dec)) { state.ra = ra * Math.PI / 180; state.dec = dec * Math.PI / 180; state.fov = 20; cutoutMark = radec(state.ra, state.dec); state.sel = -1; }
  }
  resize(); hudDirty = true; requestAnimationFrame(frame);
}
window.claude?.hot?.snapshot?.(() => ({ year: state.year, ra: state.ra, dec: state.dec, fov: state.fov, sel: state.sel, layers: state.layers, tab: state.tab }));
window.claude?.hot?.ready ? window.claude.hot.ready(start) : start(window.claude?.hot?.data ?? {});
