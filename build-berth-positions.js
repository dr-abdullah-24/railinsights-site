'use strict';
/**
 * build-berth-positions.js
 *
 * Outputs data/berth-positions.json
 *   TD:berth → { lat, lon, name, platform, line, bearing, source }
 *
 * Sources:
 *   data/berths-geo.json  — pre-processed berth coordinates (NR + CORPUS)
 *   data/SMART.json       — Network Rail berth step data
 *   OSM Overpass API      — platform geometry for position refinement
 *
 * Run: node build-berth-positions.js
 */

const fs   = require('fs');
const https = require('https');
const path  = require('path');

const DATA = path.join(__dirname, 'data');
const OUT  = path.join(DATA, 'berth-positions.json');

// ── Schema indices for berths-geo.json rows ──────────────────────────────────
// [td, berth, station, stanox, platform, event, route, fromLine, toLine,
//  offsetMetres, sourceDate, locationName, tiploc, crs, nlc,
//  easting, northing, latitude, longitude]
const I = {
  td: 0, berth: 1, stanox: 3, platform: 4,
  fromLine: 7, locationName: 11, lat: 17, lon: 18,
};

// ── 1. Load berths-geo.json ──────────────────────────────────────────────────
console.log('[1/5] Loading berths-geo.json...');
const berthsGeo = JSON.parse(fs.readFileSync(path.join(DATA, 'berths-geo.json'), 'utf8'));

const berthCoords = new Map(); // `${td}:${berth}` → {lat,lon,name,platform,stanox,line}
for (const r of berthsGeo) {
  const lat = r[I.lat];
  const lon = r[I.lon];
  if (!lat || !lon) continue;
  const key = `${r[I.td]}:${r[I.berth]}`;
  if (!berthCoords.has(key)) {
    berthCoords.set(key, {
      lat,
      lon,
      name:     r[I.locationName] || '',
      platform: r[I.platform]     ? String(r[I.platform]).trim() : '',
      stanox:   r[I.stanox]       ? String(r[I.stanox]).trim()   : '',
      line:     r[I.fromLine]     ? String(r[I.fromLine]).trim() : '',
      td:       r[I.td],
      berth:    r[I.berth],
    });
  }
}
console.log(`  → ${berthCoords.size} unique berths with coordinates`);

// ── 2. Load SMART.json — build berth graph + enrich platform/line ────────────
console.log('[2/5] Loading SMART.json...');
const smart = JSON.parse(fs.readFileSync(path.join(DATA, 'SMART.json'), 'utf8'));

const fwdGraph = new Map(); // `${td}:${fromBerth}` → Set<`${td}:${toBerth}`>
const revGraph = new Map(); // `${td}:${toBerth}`   → Set<`${td}:${fromBerth}`>

for (const s of smart) {
  const td = s.TD;
  if (!td || !s.FROMBERTH || !s.TOBERTH) continue;
  const from = `${td}:${s.FROMBERTH}`;
  const to   = `${td}:${s.TOBERTH}`;

  if (!fwdGraph.has(from)) fwdGraph.set(from, new Set());
  fwdGraph.get(from).add(to);
  if (!revGraph.has(to)) revGraph.set(to, new Set());
  revGraph.get(to).add(from);

  // Fill in platform / line where berths-geo left it blank
  const rec = berthCoords.get(from);
  if (rec) {
    if (!rec.platform && s.PLATFORM) rec.platform = String(s.PLATFORM).trim();
    if (!rec.line     && s.FROMLINE) rec.line      = String(s.FROMLINE).trim();
  }
}
console.log(`  → ${fwdGraph.size} from-berths with forward steps in graph`);

// ── 3. Compute track bearings ────────────────────────────────────────────────
console.log('[3/5] Computing track bearings from berth graph...');

const toRad = d => d * Math.PI / 180;
const toDeg = r => r * 180 / Math.PI;

function geodesicBearing(lat1, lon1, lat2, lon2) {
  const dLon = toRad(lon2 - lon1);
  const φ1   = toRad(lat1), φ2 = toRad(lat2);
  const y = Math.sin(dLon) * Math.cos(φ2);
  const x = Math.cos(φ1) * Math.sin(φ2) - Math.sin(φ1) * Math.cos(φ2) * Math.cos(dLon);
  return (toDeg(Math.atan2(y, x)) + 360) % 360;
}

function circularMean(angles) {
  let sinSum = 0, cosSum = 0;
  for (const a of angles) { sinSum += Math.sin(toRad(a)); cosSum += Math.cos(toRad(a)); }
  return (toDeg(Math.atan2(sinSum, cosSum)) + 360) % 360;
}

const berthBearings = new Map();
for (const [key, rec] of berthCoords) {
  const angles = [];

  const fwds = fwdGraph.get(key);
  if (fwds) {
    for (const toKey of fwds) {
      const dest = berthCoords.get(toKey);
      if (dest) angles.push(geodesicBearing(rec.lat, rec.lon, dest.lat, dest.lon));
    }
  }
  // If no forward neighbours use reverse direction (same track, opposite arrow)
  if (angles.length === 0) {
    const revs = revGraph.get(key);
    if (revs) {
      for (const fromKey of revs) {
        const src = berthCoords.get(fromKey);
        if (src) angles.push(geodesicBearing(src.lat, src.lon, rec.lat, rec.lon));
      }
    }
  }

  if (angles.length > 0) {
    berthBearings.set(key, Math.round(circularMean(angles)));
  }
}
console.log(`  → ${berthBearings.size} berths with bearing computed`);

// ── 4. Fetch OSM railway platform geometry ───────────────────────────────────
console.log('[4/5] Querying Overpass API for UK railway platforms...');

const OSM_QUERY = `[out:json][timeout:90];
(
  node["railway"="platform"]["ref"](49.8,-8.7,60.9,1.9);
  way["railway"="platform"]["ref"](49.8,-8.7,60.9,1.9);
  relation["railway"="platform"]["ref"](49.8,-8.7,60.9,1.9);
);
out center;`;

function httpsPost(url, body) {
  return new Promise((resolve, reject) => {
    const u = new URL(url);
    const opts = {
      hostname: u.hostname,
      path:     u.pathname,
      method:   'POST',
      timeout:  120_000,
      headers: {
        'Content-Type':   'application/x-www-form-urlencoded',
        'Content-Length': Buffer.byteLength(body),
        'User-Agent':     'RailInsights/1.0 build-berth-positions.js',
      },
    };
    const req = https.request(opts, res => {
      const chunks = [];
      res.on('data', c => chunks.push(c));
      res.on('end', () => resolve(Buffer.concat(chunks).toString('utf8')));
    });
    req.on('timeout', () => { req.destroy(); reject(new Error('request timed out')); });
    req.on('error', reject);
    req.write(body);
    req.end();
  });
}

function haversine(lat1, lon1, lat2, lon2) {
  const R = 6_371_000;
  const φ1 = toRad(lat1), φ2 = toRad(lat2);
  const Δφ = toRad(lat2 - lat1), Δλ = toRad(lon2 - lon1);
  const a  = Math.sin(Δφ / 2) ** 2 + Math.cos(φ1) * Math.cos(φ2) * Math.sin(Δλ / 2) ** 2;
  return 2 * R * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
}

async function fetchOSMPlatforms() {
  const body = 'data=' + encodeURIComponent(OSM_QUERY);
  const raw  = await httpsPost('https://overpass-api.de/api/interpreter', body);
  const data = JSON.parse(raw);

  // Group by platform ref for fast lookup
  const byRef = new Map();
  for (const el of data.elements) {
    const lat = el.center?.lat ?? el.lat;
    const lon = el.center?.lon ?? el.lon;
    if (!lat || !lon) continue;
    const ref = String(el.tags?.ref || '').trim();
    if (!ref) continue;
    if (!byRef.has(ref)) byRef.set(ref, []);
    byRef.get(ref).push({ lat, lon });
  }
  console.log(`  → ${data.elements.length} OSM elements, ${byRef.size} distinct platform refs`);
  return byRef;
}

// ── 5. Assemble and write output ─────────────────────────────────────────────
async function main() {
  let osmByRef = new Map();
  try {
    osmByRef = await fetchOSMPlatforms();
  } catch (err) {
    console.warn(`  ⚠ OSM fetch failed (${err.message}) — skipping platform position refinement`);
  }

  console.log('[5/5] Building berth-positions.json...');

  const output    = {};
  let osmRefined  = 0;
  const THRESHOLD = 500; // metres — max distance for OSM platform match

  for (const [key, rec] of berthCoords) {
    let lat = rec.lat, lon = rec.lon, source = 'berths-geo';

    if (rec.platform && osmByRef.size) {
      const candidates = osmByRef.get(rec.platform) || [];
      let bestDist = THRESHOLD, bestPt = null;
      for (const pt of candidates) {
        const d = haversine(rec.lat, rec.lon, pt.lat, pt.lon);
        if (d < bestDist) { bestDist = d; bestPt = pt; }
      }
      if (bestPt) { lat = bestPt.lat; lon = bestPt.lon; source = 'osm'; osmRefined++; }
    }

    output[key] = {
      lat,
      lon,
      name:     rec.name     || null,
      platform: rec.platform || null,
      line:     rec.line     || null,
      bearing:  berthBearings.get(key) ?? null,
      source,
    };
  }

  fs.writeFileSync(OUT, JSON.stringify(output));
  console.log(`  → ${Object.keys(output).length} berths written`);
  console.log(`  → ${osmRefined} positions refined to OSM platform geometry`);
  console.log(`\nDone → ${OUT}`);
}

main().catch(err => {
  console.error('[Fatal]', err.message);
  process.exit(1);
});
