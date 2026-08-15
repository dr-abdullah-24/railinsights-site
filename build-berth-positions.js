'use strict';
/**
 * build-berth-positions.js
 *
 * When used as Render's start command:
 *   1. Regenerates berth-snap.js / berth-positions.json if they are missing
 *   2. Hands off to server.js (keeps the process alive)
 *
 * When run manually:
 *   node build-berth-positions.js --build-only          (build then exit)
 *   node build-berth-positions.js --build-only --force  (always regenerate)
 *
 * Snap source: OSM Overpass (same WGS84 coordinates as OpenRailwayMap tiles)
 * Fallback:    freight-network.js (if Overpass is unreachable)
 */

const fs    = require('fs');
const vm    = require('vm');
const https = require('https');
const path  = require('path');

const DATA       = path.join(__dirname, 'data');
const SNAP_OUT   = path.join(DATA, 'berth-snap.js');
const POS_OUT    = path.join(DATA, 'berth-positions.json');
const BUILD_ONLY = process.argv.includes('--build-only');
const FORCE      = process.argv.includes('--force');

// ── Helpers ───────────────────────────────────────────────────────────────────
const toRad = d => d * Math.PI / 180;
const toDeg = r => r * 180 / Math.PI;

function geoBearing(lat1,lon1,lat2,lon2){
  const dL=toRad(lon2-lon1),φ1=toRad(lat1),φ2=toRad(lat2);
  return(toDeg(Math.atan2(Math.sin(dL)*Math.cos(φ2),Math.cos(φ1)*Math.sin(φ2)-Math.sin(φ1)*Math.cos(φ2)*Math.cos(dL)))+360)%360;
}
function circularMean(angs){
  let s=0,c=0;for(const a of angs){s+=Math.sin(toRad(a));c+=Math.cos(toRad(a));}
  return(toDeg(Math.atan2(s,c))+360)%360;
}

// Grid spatial index over line segments
function buildGrid(segments, cellDeg) {
  const grid = new Map();
  for (const seg of segments) {
    const {x1,y1,x2,y2} = seg;
    const ciMin=Math.floor(Math.min(x1,x2)/cellDeg)-1, ciMax=Math.floor(Math.max(x1,x2)/cellDeg)+1;
    const cjMin=Math.floor(Math.min(y1,y2)/cellDeg)-1, cjMax=Math.floor(Math.max(y1,y2)/cellDeg)+1;
    for (let ci=ciMin;ci<=ciMax;ci++) for(let cj=cjMin;cj<=cjMax;cj++){
      const k=`${ci},${cj}`;
      if(!grid.has(k))grid.set(k,[]);
      grid.get(k).push(seg);
    }
  }
  return grid;
}

function nearestOnGrid(grid, cellDeg, lat, lon, maxDeg) {
  const ci0=Math.floor(lon/cellDeg), cj0=Math.floor(lat/cellDeg);
  const r=Math.ceil(maxDeg/cellDeg)+1;
  let best=null;
  for(let di=-r;di<=r;di++) for(let dj=-r;dj<=r;dj++){
    const segs=grid.get(`${ci0+di},${cj0+dj}`);
    if(!segs)continue;
    for(const {x1,y1,x2,y2} of segs){
      const dx=x2-x1,dy=y2-y1,lenSq=dx*dx+dy*dy;
      let nx,ny;
      if(lenSq<1e-14){nx=x1;ny=y1;}
      else{const t=Math.max(0,Math.min(1,((lon-x1)*dx+(lat-y1)*dy)/lenSq));nx=x1+t*dx;ny=y1+t*dy;}
      const d=(lon-nx)*(lon-nx)+(lat-ny)*(lat-ny);
      if(!best||d<best.d)best={lat:ny,lon:nx,d};
    }
  }
  if(best&&best.d<maxDeg*maxDeg)return{lat:best.lat,lon:best.lon,snapped:true};
  return{lat,lon,snapped:false};
}

// ── OSM Overpass fetch (chunked by latitude band) ────────────────────────────
// Full UK in one request is too large; split into 6 ~2° lat bands.
// Uses same WGS84 coordinates as ORM tiles → perfect visual alignment.
const UK_BANDS = [
  [59.0, 61.0], [57.0, 59.0], [55.0, 57.0],
  [53.0, 55.0], [51.0, 53.0], [49.8, 51.0],
];
const UK_LON = [-8.7, 1.9];

function httpsPost(url, body) {
  return new Promise((resolve, reject) => {
    const u = new URL(url);
    const req = https.request({
      hostname: u.hostname, path: u.pathname, method: 'POST', timeout: 60_000,
      headers: {
        'Content-Type':   'application/x-www-form-urlencoded',
        'Content-Length': Buffer.byteLength(body),
        'User-Agent':     'RailInsights/1.0 build-berth-positions.js',
      },
    }, res => {
      const chunks = [];
      res.on('data', c => chunks.push(c));
      res.on('end', () => resolve(Buffer.concat(chunks).toString('utf8')));
    });
    req.on('timeout', () => { req.destroy(); reject(new Error('timeout')); });
    req.on('error', reject);
    req.write(body);
    req.end();
  });
}

async function fetchBand(latMin, latMax) {
  const bbox = `${latMin},${UK_LON[0]},${latMax},${UK_LON[1]}`;
  const q    = `[out:json][timeout:55];way["railway"~"^(rail|light_rail|narrow_gauge)$"](${bbox});out geom;`;
  const raw  = await httpsPost('https://overpass-api.de/api/interpreter', 'data=' + encodeURIComponent(q));
  const data = JSON.parse(raw); // throws if Overpass returned XML error
  const segs = [];
  for (const el of data.elements) {
    if (!el.geometry || el.geometry.length < 2) continue;
    for (let i = 0; i < el.geometry.length - 1; i++) {
      const a = el.geometry[i], b = el.geometry[i+1];
      if (a && b) segs.push({ x1:a.lon, y1:a.lat, x2:b.lon, y2:b.lat });
    }
  }
  return segs;
}

async function fetchOsmSegments() {
  console.log(`  Fetching ${UK_BANDS.length} latitude bands from Overpass (~2 min)...`);
  const allSegs = [];
  for (const [latMin, latMax] of UK_BANDS) {
    process.stdout.write(`  Band ${latMin}–${latMax}°N ... `);
    const segs = await fetchBand(latMin, latMax);
    allSegs.push(...segs);
    console.log(`${segs.length} segments`);
    await new Promise(r => setTimeout(r, 2000)); // polite pause between requests
  }
  console.log(`  → ${allSegs.length} total OSM track segments`);
  return allSegs;
}

// ── Fallback: freight-network.js ──────────────────────────────────────────────
function loadFreightSegments() {
  console.log('  Loading freight-network.js as fallback snap source...');
  const freightJs  = fs.readFileSync(path.join(DATA, 'freight-network.js'), 'utf8');
  const ctx = { window: {} };
  vm.createContext(ctx);
  vm.runInContext(freightJs, ctx);
  const gj = ctx.window.FREIGHT_GEOJSON;
  const segs = [];
  if (gj && gj.features) {
    for (const f of gj.features) {
      const coords = f.geometry.coordinates;
      for (let i = 0; i < coords.length-1; i++)
        segs.push({ x1:coords[i][0], y1:coords[i][1], x2:coords[i+1][0], y2:coords[i+1][1] });
    }
  }
  console.log(`  → ${segs.length} segments from freight-network.js`);
  return segs;
}

// ── Main build ────────────────────────────────────────────────────────────────
async function runBuild() {
  // 1. berth-locations.js — primary source (curated per-location precision)
  //    berths-geo.json uses STANOX-area centroids which can be several km wrong;
  //    berth-locations.js maps berths to named signal/junction locations.
  console.log('[1/5] Loading berth-locations.js...');
  const berthLocJs = fs.readFileSync(path.join(DATA, 'berth-locations.js'), 'utf8');
  const berthLocCtx = { window: {} };
  vm.createContext(berthLocCtx);
  vm.runInContext(berthLocJs, berthLocCtx);
  const berthLocData = berthLocCtx.window.BERTH_DATA;

  const berthCoords = new Map();
  for (const loc of berthLocData.locations) {
    for (const b of (loc.b || [])) {
      const key = `${loc.td}:${b}`;
      if (!berthCoords.has(key)) {
        berthCoords.set(key, {
          lat: +loc.lat, lon: +loc.lon,
          name:     loc.n  || '',
          platform: '',
          stanox:   '',
          line:     '',
          td: loc.td, berth: b,
        });
      }
    }
  }
  console.log(`  → ${berthCoords.size} unique berths`);

  // 2. SMART.json — step graph
  console.log('[2/5] Loading SMART.json...');
  const _smartRaw = JSON.parse(fs.readFileSync(path.join(DATA, 'SMART.json'), 'utf8'));
  const smart = Array.isArray(_smartRaw) ? _smartRaw : (_smartRaw.BERTHDATA || Object.values(_smartRaw)[0]);

  const fwdGraph = new Map(), revGraph = new Map();
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
  console.log(`  → ${fwdGraph.size} from-berths in graph`);

  // 3. Track segments for snapping (OSM preferred, freight fallback)
  console.log('[3/5] Fetching track geometry for snapping...');
  let trackSegs;
  let snapSource;
  try {
    trackSegs  = await fetchOsmSegments();
    snapSource = 'osm';
  } catch (err) {
    console.warn(`  ⚠ Overpass failed (${err.message}) — using freight-network.js`);
    trackSegs  = loadFreightSegments();
    snapSource = 'freight';
  }

  const CELL   = 0.01; // ~1km cells — OSM data is dense enough for fine cells
  const trackGrid = buildGrid(trackSegs, CELL);
  console.log(`  → ${trackSegs.length} segments indexed across ${trackGrid.size} cells [${snapSource}]`);

  // 4. Bearings from SMART graph
  console.log('[4/5] Computing track bearings...');
  const berthBearings = new Map();
  for (const [key,rec] of berthCoords) {
    const angs=[];
    for(const tk of(fwdGraph.get(key)||[])){const d=berthCoords.get(tk);if(d)angs.push(geoBearing(rec.lat,rec.lon,d.lat,d.lon));}
    if(!angs.length)for(const fk of(revGraph.get(key)||[])){const s=berthCoords.get(fk);if(s)angs.push(geoBearing(s.lat,s.lon,rec.lat,rec.lon));}
    if(angs.length)berthBearings.set(key,Math.round(circularMean(angs)));
  }
  console.log(`  → ${berthBearings.size} berths with bearing`);

  // 5. Snap and write
  console.log('[5/5] Snapping berths to nearest track and writing output...');
  const snapOut={}, fullOut={};
  let snapped=0, kept=0;
  for (const [key,rec] of berthCoords) {
    const r = nearestOnGrid(trackGrid, CELL, rec.lat, rec.lon, 0.025);
    if(r.snapped)snapped++;else kept++;
    snapOut[key] = [+r.lat.toFixed(6), +r.lon.toFixed(6)];
    fullOut[key] = {
      lat:+r.lat.toFixed(6), lon:+r.lon.toFixed(6),
      name:rec.name||null, platform:rec.platform||null, line:rec.line||null,
      bearing:berthBearings.get(key)??null,
      source:r.snapped?snapSource+'-snap':'berths-geo',
    };
  }

  fs.writeFileSync(SNAP_OUT, 'window.BERTH_SNAP='+JSON.stringify(snapOut)+';');
  fs.writeFileSync(POS_OUT,  JSON.stringify(fullOut));

  const snapKB = (fs.statSync(SNAP_OUT).size/1024).toFixed(1);
  const posKB  = (fs.statSync(POS_OUT).size/1024).toFixed(1);
  console.log(`  → ${Object.keys(snapOut).length} berths, ${snapped} snapped (${kept} kept original)`);
  console.log(`  berth-snap.js        → ${snapKB} KB  [source: ${snapSource}]`);
  console.log(`  berth-positions.json → ${posKB} KB`);
}

// ── Entry point ───────────────────────────────────────────────────────────────
const outputsExist = fs.existsSync(SNAP_OUT) && fs.existsSync(POS_OUT);

async function main() {
  if (!outputsExist || FORCE) {
    const reason = FORCE ? '--force flag' : 'output files missing';
    console.log(`[berth-positions] Running build (${reason})...`);
    await runBuild();
    console.log('[berth-positions] Build complete.');
  } else {
    console.log('[berth-positions] Snap data present — skipping build.');
  }

  if (BUILD_ONLY) {
    console.log('[berth-positions] Exiting (--build-only).');
    process.exit(0);
  }

  console.log('[berth-positions] Starting railinsights server...');
  require('./server.js');
}

main().catch(err => {
  console.error('[berth-positions] Fatal:', err.message);
  process.exit(1);
});
