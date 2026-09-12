
const assert=require('node:assert/strict'),fs=require('node:fs'),vm=require('node:vm');
const RailPosition=require('../assets/rail-position.js'),{analyse,distance}=require('../assets/congestion-model.js');
const ctx={window:{}};vm.runInNewContext(fs.readFileSync('data/berth-quality.js','utf8'),ctx);const quality=ctx.window.BERTH_QUALITY;
const now=Date.now(),base={td:'EX',berth:'D172',lat:55.50173,lon:3.47534,delay:30,delayStatus:'LATE',ts:now};
const whiteball=RailPosition.resolve(base,quality);
assert.ok(distance(whiteball,{lat:50.95365,lon:-3.29586})<.01);
assert.equal(whiteball.positionStatus,'track-estimate');
const maidenhead=RailPosition.resolve({...base,td:'D6',berth:'0569'},quality);
assert.equal(maidenhead.positionStatus,'track-estimate');
assert.ok(distance(maidenhead,{lat:51.51838,lon:-.72281})<.05);
const missing=RailPosition.resolve({...base,berth:'UNKNOWN'},quality);
assert.equal(missing.lat,null);assert.equal(missing.lon,null);assert.equal(RailPosition.canMap(missing),false);
assert.equal(RailPosition.canMap(RailPosition.resolve(base,undefined)),false);
const s=analyse({missing,whiteball},{type:'all',threshold:0,age:10},now);
assert.equal(s.late.length,2);assert.equal(s.clusters.length,0);
assert.ok(s.late.some(t=>t.positionStatus==='unavailable'));
const src=fs.readFileSync('server.js','utf8');
const begin=src.indexOf('function updateTDPosition('),end=src.indexOf('function handleTRUST(',begin);
const server={RailPosition,positionQuality:quality,trainState:new Map(),berthLookup:new Map(),isBlank:v=>!v,classifyTrain:()=> 'express'};
vm.createContext(server);vm.runInContext(src.slice(begin,end),server);
server.handleCA({area_id:'EX',to:'D172',descr:'1C77',ts:now});assert.equal(server.trainState.get('1C77').positionStatus,'track-estimate');
server.handleCA({area_id:'EX',to:'UNKNOWN',descr:'1C77',ts:now+1});assert.equal(server.trainState.get('1C77').lat,null);
server.handleCC({area_id:'EX',to:'D172',descr:'1C77',ts:now+2});assert.ok(server.trainState.get('1C77').lon<0);
server.handleCB({descr:'1C77'});assert.equal(server.trainState.size,0);
for(const p of Object.values(quality.positions)){
 if(p.status==='unavailable')assert.equal(p.lat,undefined);
 else {assert.ok(Number.isFinite(p.lat)&&Number.isFinite(p.lon));if(p.status==='route-estimate')assert.ok(p.offsetMetres<=250);}
}
console.log('Passed: Whiteball/Maidenhead corrections, no stale-coordinate reuse, unavailable locations retained in list, conservative offsets.');
