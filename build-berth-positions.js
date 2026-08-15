'use strict';
/**
 * build-berth-positions.js
 *
 * Outputs:
 *   data/berth-snap.js        — compact window.BERTH_SNAP = {TD:berth → [lat,lon]}
 *                               loaded by live-trains.html to override server coords
 *   data/berth-positions.json — full data for future tooling
 *
 * Sources:
 *   data/berths-geo.json      — NR berth locations (station centroid accuracy)
 *   data/SMART.json           — berth step graph (platform / line enrichment)
 *   data/freight-network.js   — UK rail line geometry (used for track snapping)
 *
 * Run: node build-berth-positions.js
 */

const fs   = require('fs');
const vm   = require('vm');
const path = require('path');

const DATA = path.join(__dirname, 'data');

// ── Schema indices for berths-geo.json rows ──────────────────────────────────
const I = { td:0, berth:1, stanox:3, platform:4, fromLine:7, locationName:11, lat:17, lon:18 };

// ── 1. Load berths-geo.json ──────────────────────────────────────────────────
console.log('[1/5] Loading berths-geo.json...');
const _raw = JSON.parse(fs.readFileSync(path.join(DATA, 'berths-geo.json'), 'utf8'));
const berthsGeo = Array.isArray(_raw) ? _raw : _raw.records;

const berthCoords = new Map(); // `${td}:${berth}` → object
for (const r of berthsGeo) {
  const lat = r[I.lat], lon = r[I.lon];
  if (!lat || !lon) continue;
  const key = `${r[I.td]}:${r[I.berth]}`;
  if (!berthCoords.has(key)) {
    berthCoords.set(key, {
      lat: +lat, lon: +lon,
      name:     r[I.locationName] || '',
      platform: r[I.platform]     ? String(r[I.platform]).trim() : '',
      stanox:   r[I.stanox]       ? String(r[I.stanox]).trim()   : '',
      line:     r[I.fromLine]     ? String(r[I.fromLine]).trim() : '',
      td: r[I.td], berth: r[I.berth],
    });
  }
}
console.log(`  → ${berthCoords.size} unique berths with coordinates`);

// ── 2. Load SMART.json — berth graph + enrich platform / line ────────────────
console.log('[2/5] Loading SMART.json...');
const _smartRaw = JSON.parse(fs.readFileSync(path.join(DATA, 'SMART.json'), 'utf8'));
const smart = Array.isArray(_smartRaw) ? _smartRaw : (_smartRaw.BERTHDATA || Object.values(_smartRaw)[0]);

const fwdGraph = new Map();
const revGraph = new Map();
for (const s of smart) {
  const td = s.TD;
  if (!td || !s.FROMBERTH || !s.TOBERTH) continue;
  const from = `${td}:${s.FROMBERTH}`, to = `${td}:${s.TOBERTH}`;
  if (!fwdGraph.has(from)) fwdGraph.set(from, new Set());
  fwdGraph.get(from).add(to);
  if (!revGraph.has(to)) revGraph.set(to, new Set());
  revGraph.get(to).add(from);
  const rec = berthCoords.get(from);
  if (rec) {
    if (!rec.platform && s.PLATFORM) rec.platform = String(s.PLATFORM).trim();
    if (!rec.line     && s.FROMLINE) rec.line      = String(s.FROMLINE).trim();
  }
}
console.log(`  → ${fwdGraph.size} from-berths in step graph`);

// ── 3. Build spatial grid index from freight-network.js ─────────────────────
console.log('[3/5] Building track snap index from freight-network.js...');
const freightJs  = fs.readFileSync(path.join(DATA, 'freight-network.js'), 'utf8');
const freightCtx = { window: {} };
vm.createContext(freightCtx);
vm.runInContext(freightJs, freightCtx);
const freightGJ = freightCtx.window.FREIGHT_GEOJSON;

const CELL = 0.04; // ~4km grid cells
const trackGrid = new Map(); // `ci,cj` → [{x1,y1,x2,y2}]
let segCount = 0;

function cellsForSeg(x1, y1, x2, y2) {
  const keys = [];
  const ciMin = Math.floor(Math.min(x1,x2)/CELL) - 1;
  const ciMax = Math.floor(Math.max(x1,x2)/CELL) + 1;
  const cjMin = Math.floor(Math.min(y1,y2)/CELL) - 1;
  const cjMax = Math.floor(Math.max(y1,y2)/CELL) + 1;
  for (let ci = ciMin; ci <= ciMax; ci++)
    for (let cj = cjMin; cj <= cjMax; cj++)
      keys.push(`${ci},${cj}`);
  return keys;
}

if (freightGJ && freightGJ.features) {
  for (const f of freightGJ.features) {
    const coords = f.geometry.coordinates;
    for (let i = 0; i < coords.length - 1; i++) {
      const [x1,y1] = coords[i], [x2,y2] = coords[i+1];
      const seg = {x1,y1,x2,y2};
      for (const k of cellsForSeg(x1,y1,x2,y2)) {
        if (!trackGrid.has(k)) trackGrid.set(k, []);
        trackGrid.get(k).push(seg);
      }
      segCount++;
    }
  }
}
console.log(`  → ${segCount} track segments across ${trackGrid.size} grid cells`);

function nearestOnTrack(lat, lon) {
  const MAX_DEG = 0.025; // ~2.5 km snap radius
  const ci0 = Math.floor(lon/CELL), cj0 = Math.floor(lat/CELL);
  let best = null;

  for (let di = -2; di <= 2; di++) {
    for (let dj = -2; dj <= 2; dj++) {
      const segs = trackGrid.get(`${ci0+di},${cj0+dj}`);
      if (!segs) continue;
      for (const {x1,y1,x2,y2} of segs) {
        const dx = x2-x1, dy = y2-y1, lenSq = dx*dx+dy*dy;
        let nx, ny;
        if (lenSq < 1e-14) { nx=x1; ny=y1; }
        else {
          const t = Math.max(0, Math.min(1, ((lon-x1)*dx + (lat-y1)*dy) / lenSq));
          nx = x1+t*dx; ny = y1+t*dy;
        }
        const d = (lon-nx)*(lon-nx) + (lat-ny)*(lat-ny);
        if (!best || d < best.d) best = {lat:ny, lon:nx, d};
      }
    }
  }

  if (best && best.d < MAX_DEG*MAX_DEG) return {lat: best.lat, lon: best.lon, snapped: true};
  return {lat, lon, snapped: false};
}

// ── 4. Compute track bearings ────────────────────────────────────────────────
console.log('[4/5] Computing track bearings...');
const toRad = d => d*Math.PI/180, toDeg = r => r*180/Math.PI;

function geoBearing(lat1,lon1,lat2,lon2) {
  const dL = toRad(lon2-lon1), φ1=toRad(lat1), φ2=toRad(lat2);
  const y = Math.sin(dL)*Math.cos(φ2);
  const x = Math.cos(φ1)*Math.sin(φ2) - Math.sin(φ1)*Math.cos(φ2)*Math.cos(dL);
  return (toDeg(Math.atan2(y,x))+360)%360;
}
function circularMean(angs) {
  let s=0,c=0; for(const a of angs){s+=Math.sin(toRad(a));c+=Math.cos(toRad(a));}
  return (toDeg(Math.atan2(s,c))+360)%360;
}

const berthBearings = new Map();
for (const [key, rec] of berthCoords) {
  const angs = [];
  for (const toKey of (fwdGraph.get(key)||[])) {
    const d = berthCoords.get(toKey);
    if (d) angs.push(geoBearing(rec.lat,rec.lon,d.lat,d.lon));
  }
  if (!angs.length) {
    for (const fromKey of (revGraph.get(key)||[])) {
      const s = berthCoords.get(fromKey);
      if (s) angs.push(geoBearing(s.lat,s.lon,rec.lat,rec.lon));
    }
  }
  if (angs.length) berthBearings.set(key, Math.round(circularMean(angs)));
}
console.log(`  → ${berthBearings.size} berths with bearing`);

// ── 5. Build output files ────────────────────────────────────────────────────
console.log('[5/5] Snapping berths to nearest track and writing output...');

const snapOut  = {};   // for berth-snap.js  → window.BERTH_SNAP
const fullOut  = {};   // for berth-positions.json
let snappedCount = 0, unsnapped = 0;

for (const [key, rec] of berthCoords) {
  const {lat, lon, snapped} = nearestOnTrack(rec.lat, rec.lon);
  if (snapped) snappedCount++; else unsnapped++;

  snapOut[key] = [+lat.toFixed(6), +lon.toFixed(6)];
  fullOut[key]  = {
    lat: +lat.toFixed(6),
    lon: +lon.toFixed(6),
    name:     rec.name     || null,
    platform: rec.platform || null,
    line:     rec.line     || null,
    bearing:  berthBearings.get(key) ?? null,
    source:   snapped ? 'track-snap' : 'berths-geo',
  };
}

fs.writeFileSync(path.join(DATA, 'berth-snap.js'),
  'window.BERTH_SNAP=' + JSON.stringify(snapOut) + ';');

fs.writeFileSync(path.join(DATA, 'berth-positions.json'),
  JSON.stringify(fullOut));

console.log(`  → ${Object.keys(snapOut).length} berths written`);
console.log(`  → ${snappedCount} snapped to track geometry (${unsnapped} kept original)`);
console.log(`\nDone.`);
console.log(`  berth-snap.js        → ${(fs.statSync(path.join(DATA,'berth-snap.js')).size/1024).toFixed(1)} KB`);
console.log(`  berth-positions.json → ${(fs.statSync(path.join(DATA,'berth-positions.json')).size/1024).toFixed(1)} KB`);
