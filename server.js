'use strict';
// TD live WebSocket bridge — Confluent Kafka → browser clients
// Run: npm install && npm start
// Serves: ws://localhost:3000  and  http://localhost:3000/status

const { Kafka } = require('kafkajs');
const WebSocket  = require('ws');
const http       = require('http');
const fs         = require('fs');
const vm         = require('vm');

// ── Berth lookup ─────────────────────────────────────────────────────────────
const berthJs = fs.readFileSync('./data/berth-locations.js', 'utf8');
const ctx = { window: {} };
vm.createContext(ctx);
vm.runInContext(berthJs, ctx);

const berthLookup = new Map(); // "TD:berth" → { lat, lon, name }
for (const loc of ctx.window.BERTH_DATA.locations) {
  for (const b of (loc.b || [])) {
    const key = `${loc.td}:${b}`;
    if (!berthLookup.has(key)) {
      berthLookup.set(key, { lat: loc.lat, lon: loc.lon, name: loc.n });
    }
  }
}
console.log(`[Berths] ${berthLookup.size} entries from ${ctx.window.BERTH_DATA.locations.length} locations`);

// ── Train state ───────────────────────────────────────────────────────────────
// headcode → { lat, lon, name, td, berth, type, ts }
const trainState = new Map();

function classifyTrain(hc) {
  if (!hc) return 'special';
  switch (hc[0]) {
    case '1': return 'express';
    case '2': return 'passenger';
    case '3': return 'parcels';
    case '5': return 'ecs';
    case '6': return 'loco';
    case '4': case '7': case '8': case '9': return 'freight';
    default:  return 'special';
  }
}

function isBlank(descr) {
  return !descr || !descr.trim() || descr.trim() === '0000';
}

function handleCA({ area_id, to, descr, ts }) {
  if (isBlank(descr)) return;
  const hc   = descr.trim();
  const loc  = berthLookup.get(`${area_id}:${to}`);
  const prev = trainState.get(hc);
  const n = Number(ts); const stamp = (n > 0 && n < 3e12) ? n : Date.now();

  if (loc) {
    trainState.set(hc, { lat: loc.lat, lon: loc.lon, name: loc.name, td: area_id, berth: to, type: classifyTrain(hc), ts: stamp });
  } else if (prev) {
    // berth not in data — keep old coords, update TD/berth reference
    trainState.set(hc, { ...prev, td: area_id, berth: to, ts: stamp });
  }
}

function handleCB({ descr }) {
  if (!isBlank(descr)) trainState.delete(descr.trim());
}

function handleCC({ area_id, to, descr, ts }) {
  if (isBlank(descr)) return;
  const hc  = descr.trim();
  const loc = berthLookup.get(`${area_id}:${to}`);
  if (loc) {
    trainState.set(hc, { lat: loc.lat, lon: loc.lon, name: loc.name, td: area_id, berth: to, type: classifyTrain(hc), ts: ts ? Number(ts) : Date.now() });
  }
}

// ── HTTP + WebSocket server ───────────────────────────────────────────────────
const httpServer = http.createServer((req, res) => {
  res.setHeader('Access-Control-Allow-Origin', '*');
  res.setHeader('Access-Control-Allow-Methods', 'GET');
  if (req.url === '/status') {
    res.writeHead(200, { 'Content-Type': 'application/json' });
    res.end(JSON.stringify({ trains: trainState.size, berths: berthLookup.size, uptime: Math.round(process.uptime()), clients: wss.clients.size }));
    return;
  }
  res.writeHead(404);
  res.end('Not found');
});

const wss = new WebSocket.Server({ server: httpServer });

function buildPayload() {
  return JSON.stringify({ type: 'state', trains: Object.fromEntries(trainState), count: trainState.size, ts: Date.now() });
}

let broadcastPending = false;
function scheduleBroadcast() {
  if (broadcastPending || wss.clients.size === 0) return;
  broadcastPending = true;
  setTimeout(() => {
    broadcastPending = false;
    const payload = buildPayload();
    for (const client of wss.clients) {
      if (client.readyState === WebSocket.OPEN) {
        try { client.send(payload); } catch {}
      }
    }
  }, 400); // max ~2.5 broadcasts/sec regardless of message volume
}

wss.on('connection', (ws, req) => {
  console.log(`[WS] + client (${wss.clients.size} total) from ${req.socket.remoteAddress}`);
  try { ws.send(buildPayload()); } catch {}
  ws.on('close', () => console.log(`[WS] - client (${wss.clients.size} remaining)`));
  ws.on('error', () => {});
});

const PORT = process.env.PORT || 3000;
httpServer.listen(PORT, () => {
  console.log(`[HTTP] http://localhost:${PORT}/status`);
  console.log(`[WS]   ws://localhost:${PORT}`);
});

// ── Kafka consumer ────────────────────────────────────────────────────────────
const kafka = new Kafka({
  clientId: 'railinsights-live-td',
  brokers: ['pkc-z3p1v0.europe-west2.gcp.confluent.cloud:9092'],
  ssl: true,
  sasl: { mechanism: 'plain', username: 'Y2RFWHM6GEJLEBKG', password: 'cflt7XBjbZNiqZhujMYJL6nvRog+c26a+++1DGqBdItzLK7m/sLJuCBga02d6S5w' },
  retry: { retries: 20, initialRetryTime: 3000, factor: 1.5, maxRetryTime: 60000 }
});

const consumer = kafka.consumer({
  groupId: 'SC-17cd3bcc-d000-4ffb-b5d6-2fe7a4742228',
  sessionTimeout: 30000,
  heartbeatInterval: 3000,
  readUncommitted: false
});

let msgCount = 0, logAt = Date.now();

async function startConsumer() {
  console.log('[Kafka] Connecting to Confluent Cloud...');
  await consumer.connect();
  await consumer.subscribe({ topic: 'TD_ALL_SIG_AREA', fromBeginning: false });
  console.log('[Kafka] Subscribed — waiting for TD messages...');

  await consumer.run({
    eachMessage: async ({ message }) => {
      msgCount++;
      const now = Date.now();
      if (now - logAt > 15_000) {
        console.log(`[Kafka] ${msgCount} msgs/15s | trains tracked: ${trainState.size} | WS clients: ${wss.clients.size}`);
        msgCount = 0; logAt = now;
      }

      try {
        const arr = JSON.parse(message.value.toString());
        let changed = false;
        for (const item of arr) {
          if      (item.CA_MSG) { handleCA(item.CA_MSG); changed = true; }
          else if (item.CB_MSG) { handleCB(item.CB_MSG); changed = true; }
          else if (item.CC_MSG) { handleCC(item.CC_MSG); changed = true; }
          // CT_MSG = heartbeat, ignore
        }
        if (changed) scheduleBroadcast();
      } catch { /* malformed JSON — skip silently */ }
    }
  });
}

startConsumer().catch(err => {
  console.error('[Kafka] Fatal error:', err.message);
  process.exit(1);
});

const shutdown = async () => {
  console.log('\n[Server] Shutting down gracefully...');
  await consumer.disconnect().catch(() => {});
  httpServer.close(() => process.exit(0));
  setTimeout(() => process.exit(0), 5000);
};
process.on('SIGTERM', shutdown);
process.on('SIGINT',  shutdown);
