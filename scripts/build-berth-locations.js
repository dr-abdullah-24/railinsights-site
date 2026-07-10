const fs = require('fs');
const path = require('path');

const root = path.resolve(__dirname, '..');
const geoPath = path.join(root, 'data', 'berths-geo.json');
const outPath = path.join(root, 'data', 'berth-locations.js');

function main() {
  const { schema, records } = JSON.parse(fs.readFileSync(geoPath, 'utf8'));
  const IDX = Object.fromEntries(schema.map((k, i) => [k, i]));

  const locationMap = new Map();

  for (const rec of records) {
    const lat = rec[IDX.latitude];
    const lon = rec[IDX.longitude];
    if (lat === null || lon === null) continue;

    const td = rec[IDX.td];
    const stanox = rec[IDX.stanox];
    const key = `${td}|${stanox}`;

    if (!locationMap.has(key)) {
      locationMap.set(key, {
        lat: Number(lat.toFixed(5)),
        lon: Number(lon.toFixed(5)),
        n: rec[IDX.locationName],
        tip: rec[IDX.tiploc],
        crs: rec[IDX.crs],
        sx: stanox,
        td,
        plt: rec[IDX.platform],
        rt: rec[IDX.route],
        fl: rec[IDX.fromLine],
        tl: rec[IDX.toLine],
        dt: rec[IDX.sourceDate],
        berths: new Set(),
      });
    }

    locationMap.get(key).berths.add(rec[IDX.berth]);
  }

  let totalBerths = 0;
  const locations = [...locationMap.values()].map(({ berths, ...meta }) => {
    const b = [...berths].sort();
    totalBerths += b.length;
    return { ...meta, b };
  });

  fs.writeFileSync(outPath, `window.BERTH_DATA=${JSON.stringify({ locations })}`);
  console.log(`Wrote ${locations.length} locations, ${totalBerths} unique berths to ${path.relative(root, outPath)}`);
}

main();
