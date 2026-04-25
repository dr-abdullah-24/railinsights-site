const fs = require('fs');
const path = require('path');

const root = path.resolve(__dirname, '..');
const dataDir = path.join(root, 'data');

const smartPath = path.join(dataDir, 'SMART.json');
const corpusPath = path.join(dataDir, 'CORPUS.json');
const locPath = path.join(dataDir, 'LOC_Data.txt');
const outPath = path.join(dataDir, 'berths-geo.json');

function clean(value) {
  return String(value ?? '').trim();
}

function toNumber(value) {
  const number = Number(clean(value));
  return Number.isFinite(number) ? number : null;
}

function osgbToWgs84(E, N) {
  const deg = 180 / Math.PI;
  const rad = Math.PI / 180;

  const a = 6377563.396;
  const b = 6356256.909;
  const F0 = 0.9996012717;
  const lat0 = 49 * rad;
  const lon0 = -2 * rad;
  const N0 = -100000;
  const E0 = 400000;
  const e2 = 1 - (b * b) / (a * a);
  const n = (a - b) / (a + b);

  let lat = lat0;
  let M = 0;
  do {
    lat = (N - N0 - M) / (a * F0) + lat;
    const Ma = (1 + n + (5 / 4) * n ** 2 + (5 / 4) * n ** 3) * (lat - lat0);
    const Mb = (3 * n + 3 * n ** 2 + (21 / 8) * n ** 3) * Math.sin(lat - lat0) * Math.cos(lat + lat0);
    const Mc = ((15 / 8) * n ** 2 + (15 / 8) * n ** 3) * Math.sin(2 * (lat - lat0)) * Math.cos(2 * (lat + lat0));
    const Md = (35 / 24) * n ** 3 * Math.sin(3 * (lat - lat0)) * Math.cos(3 * (lat + lat0));
    M = b * F0 * (Ma - Mb + Mc - Md);
  } while (Math.abs(N - N0 - M) >= 0.00001);

  const sinLat = Math.sin(lat);
  const cosLat = Math.cos(lat);
  const tanLat = Math.tan(lat);
  const nu = a * F0 / Math.sqrt(1 - e2 * sinLat ** 2);
  const rho = a * F0 * (1 - e2) / (1 - e2 * sinLat ** 2) ** 1.5;
  const eta2 = nu / rho - 1;
  const dE = E - E0;
  const secLat = 1 / cosLat;

  const VII = tanLat / (2 * rho * nu);
  const VIII = tanLat / (24 * rho * nu ** 3) * (5 + 3 * tanLat ** 2 + eta2 - 9 * tanLat ** 2 * eta2);
  const IX = tanLat / (720 * rho * nu ** 5) * (61 + 90 * tanLat ** 2 + 45 * tanLat ** 4);
  const X = secLat / nu;
  const XI = secLat / (6 * nu ** 3) * (nu / rho + 2 * tanLat ** 2);
  const XII = secLat / (120 * nu ** 5) * (5 + 28 * tanLat ** 2 + 24 * tanLat ** 4);
  const XIIA = secLat / (5040 * nu ** 7) * (61 + 662 * tanLat ** 2 + 1320 * tanLat ** 4 + 720 * tanLat ** 6);

  const latOsgb = lat - VII * dE ** 2 + VIII * dE ** 4 - IX * dE ** 6;
  const lonOsgb = lon0 + X * dE - XI * dE ** 3 + XII * dE ** 5 - XIIA * dE ** 7;

  return helmertOsgb36ToWgs84(latOsgb, lonOsgb, 0).map(v => Number(v.toFixed(6)));
}

function helmertOsgb36ToWgs84(lat, lon, height) {
  const a1 = 6377563.396;
  const b1 = 6356256.909;
  const e21 = 1 - (b1 * b1) / (a1 * a1);
  const nu1 = a1 / Math.sqrt(1 - e21 * Math.sin(lat) ** 2);

  let x1 = (nu1 + height) * Math.cos(lat) * Math.cos(lon);
  let y1 = (nu1 + height) * Math.cos(lat) * Math.sin(lon);
  let z1 = ((1 - e21) * nu1 + height) * Math.sin(lat);

  const tx = 446.448;
  const ty = -125.157;
  const tz = 542.06;
  const s = -20.4894e-6;
  const rx = (0.1502 / 3600) * Math.PI / 180;
  const ry = (0.2470 / 3600) * Math.PI / 180;
  const rz = (0.8421 / 3600) * Math.PI / 180;

  const x2 = tx + x1 * (1 + s) - y1 * rz + z1 * ry;
  const y2 = ty + x1 * rz + y1 * (1 + s) - z1 * rx;
  const z2 = tz - x1 * ry + y1 * rx + z1 * (1 + s);

  const a2 = 6378137.0;
  const b2 = 6356752.3141;
  const e22 = 1 - (b2 * b2) / (a2 * a2);
  const p = Math.sqrt(x2 * x2 + y2 * y2);
  let lat2 = Math.atan2(z2, p * (1 - e22));
  let latPrev;

  do {
    latPrev = lat2;
    const nu2 = a2 / Math.sqrt(1 - e22 * Math.sin(lat2) ** 2);
    lat2 = Math.atan2(z2 + e22 * nu2 * Math.sin(lat2), p);
  } while (Math.abs(lat2 - latPrev) > 1e-12);

  const lon2 = Math.atan2(y2, x2);
  return [lat2 * 180 / Math.PI, lon2 * 180 / Math.PI];
}

function loadLocations() {
  const rows = fs.readFileSync(locPath, 'utf8').split(/\r?\n/).filter(Boolean);
  const byStanox = new Map();
  const byTiploc = new Map();

  for (const row of rows) {
    const parts = row.split('\t');
    const tiploc = clean(parts[2]);
    const name = clean(parts[3]);
    const easting = toNumber(parts[6]);
    const northing = toNumber(parts[7]);
    const stanox = clean(parts[10]);
    const hasGrid = easting > 0 && northing > 0 && easting < 900000 && northing < 1300000;
    const location = {
      tiploc,
      name,
      stanox,
      easting: hasGrid ? easting : null,
      northing: hasGrid ? northing : null,
      latlon: hasGrid ? osgbToWgs84(easting, northing) : null,
    };

    if (stanox) {
      if (!byStanox.has(stanox)) byStanox.set(stanox, []);
      byStanox.get(stanox).push(location);
    }
    if (tiploc) byTiploc.set(tiploc, location);
  }

  return { byStanox, byTiploc, count: rows.length };
}

function loadCorpus() {
  const rows = JSON.parse(fs.readFileSync(corpusPath, 'utf8')).TIPLOCDATA || [];
  const byStanox = new Map();
  for (const row of rows) {
    const stanox = clean(row.STANOX);
    if (!stanox || stanox === '0') continue;
    if (!byStanox.has(stanox)) byStanox.set(stanox, []);
    byStanox.get(stanox).push({
      nlc: clean(row.NLC),
      tiploc: clean(row.TIPLOC),
      crs: clean(row['3ALPHA']),
      uic: clean(row.UIC),
      description: clean(row.NLCDESC),
      shortDescription: clean(row.NLCDESC16),
    });
  }
  return { byStanox, count: rows.length };
}

function addBerth(rows, source) {
  const berth = clean(source.berth);
  if (!berth) return;

  const key = [
    clean(source.td),
    berth,
    clean(source.stanox),
    clean(source.platform),
    clean(source.event),
    clean(source.offset),
  ].join('|');

  if (!rows.has(key)) {
    rows.set(key, {
      td: clean(source.td),
      berth,
      station: clean(source.station),
      stanox: clean(source.stanox),
      platform: clean(source.platform),
      event: clean(source.event),
      stepType: clean(source.stepType),
      route: clean(source.route),
      fromLine: clean(source.fromLine),
      toLine: clean(source.toLine),
      offsetMetres: clean(source.offset),
      sourceDate: clean(source.sourceDate),
      count: 0,
    });
  }
  rows.get(key).count += 1;
}

function main() {
  const smartRows = JSON.parse(fs.readFileSync(smartPath, 'utf8')).BERTHDATA || [];
  const locations = loadLocations();
  const corpus = loadCorpus();
  const rows = new Map();

  for (const row of smartRows) {
    const base = {
      td: row.TD,
      station: row.STANME,
      stanox: row.STANOX,
      platform: row.PLATFORM,
      event: row.EVENT,
      stepType: row.STEPTYPE,
      route: row.ROUTE,
      fromLine: row.FROMLINE,
      toLine: row.TOLINE,
      offset: row.BERTHOFFSET,
      sourceDate: row.COMMENT,
    };
    addBerth(rows, { ...base, berth: row.FROMBERTH });
  }

  const records = [...rows.values()].map(row => {
    const locs = locations.byStanox.get(row.stanox) || [];
    const validLoc = locs.find(loc => loc.latlon);
    const corpusRows = corpus.byStanox.get(row.stanox) || [];
    const corpusPrimary = corpusRows[0] || {};
    return [
      row.td,
      row.berth,
      row.station,
      row.stanox,
      row.platform,
      row.event,
      row.route,
      row.fromLine,
      row.toLine,
      row.offsetMetres,
      row.sourceDate,
      clean(validLoc?.name) || clean(corpusPrimary.description) || row.station,
      clean(validLoc?.tiploc) || clean(corpusPrimary.tiploc),
      clean(corpusPrimary.crs),
      clean(corpusPrimary.nlc),
      validLoc?.easting ?? null,
      validLoc?.northing ?? null,
      validLoc?.latlon?.[0] ?? null,
      validLoc?.latlon?.[1] ?? null,
    ];
  }).sort((a, b) =>
    a[0].localeCompare(b[0]) ||
    a[1].localeCompare(b[1]) ||
    a[12].localeCompare(b[12])
  );

  const payload = {
    generatedAt: new Date().toISOString(),
    sourceFiles: ['data/SMART.json', 'data/LOC_Data.txt', 'data/CORPUS.json'],
    schema: ['td', 'berth', 'station', 'stanox', 'platform', 'event', 'route', 'fromLine', 'toLine', 'offsetMetres', 'sourceDate', 'locationName', 'tiploc', 'crs', 'nlc', 'easting', 'northing', 'latitude', 'longitude'],
    summary: {
      sourceBerthRows: smartRows.length,
      derivedRows: records.length,
      rowsWithCoordinates: records.filter(row => row[17] !== null && row[18] !== null).length,
      locationRows: locations.count,
      corpusRows: corpus.count,
    },
    records,
  };

  fs.writeFileSync(outPath, JSON.stringify(payload));
  console.log(`Wrote ${records.length} berth-location records to ${path.relative(root, outPath)}`);
}

main();
