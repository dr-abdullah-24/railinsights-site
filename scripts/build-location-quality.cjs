// Deterministic, offline audit. Never silently turn a large correction into a train position.
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const root=path.resolve(__dirname,'..'),read=n=>JSON.parse(fs.readFileSync(path.join(root,'data',n),'utf8'));
const {distance}=require('../assets/congestion-model.js');
const ctx={window:{}};vm.runInNewContext(fs.readFileSync(path.join(root,'data/berth-locations.js'),'utf8'),ctx);
const locations=ctx.window.BERTH_DATA.locations,candidates=read('berth-positions.json'),overrides=read('location-overrides.json'),osm=read('railway-location-checks.json');
function project(p,a,b){const scale=Math.cos(p.lat*Math.PI/180),dx=(b.lon-a.lon)*scale,dy=b.lat-a.lat;
 const u=Math.max(0,Math.min(1,(((p.lon-a.lon)*scale)*dx+(p.lat-a.lat)*dy)/(dx*dx+dy*dy||1)));
 return {lat:a.lat+u*(b.lat-a.lat),lon:a.lon+u*(b.lon-a.lon)};
}
function corrected(o){let best;
 for(const w of osm.elements){if(w.tags.name!==o.route||(o.osmWay&&String(w.id)!==o.osmWay))continue;
 for(let i=1;i<w.geometry.length;i++){const p=project(o,w.geometry[i-1],w.geometry[i]),d=distance(o,p)*1000;if(!best||d<best.d)best={...p,d,way:String(w.id)};}}
 if(!best||best.d>100)throw Error('Override does not match its documented railway');return best;
}
const groups=new Map();for(const loc of locations)for(const b of loc.b||[]){const key=loc.td+':'+b;if(!groups.has(key))groups.set(key,[]);groups.get(key).push(loc);}
const positions={},summary={};
for(const [key,locs]of groups){const loc=locs[0],o=overrides[loc.tip],p=candidates[key];let r;
 if(locs.some(other=>distance(loc,other)>.5))r={status:'unavailable',reason:'Berth has conflicting named locations'};
 else if(o){const q=corrected(o);r={lat:+q.lat.toFixed(6),lon:+q.lon.toFixed(6),status:'track-estimate',label:loc.tip==='TIVITNL'?'Mapped tunnel-area estimate':'Mapped track estimate',source:'OpenStreetMap',sourceUrl:'https://www.openstreetmap.org/way/'+q.way,reason:o.reason,offsetMetres:Math.round(distance(loc,q)*1000)};}
 else if(p&&['osm-snap','freight-snap'].includes(p.source)&&p.lat>=49.8&&p.lat<=61&&p.lon>=-8.7&&p.lon<=2&&distance(loc,p)*1000<=250){r={lat:p.lat,lon:p.lon,status:p.source==='osm-snap'?'track-estimate':'route-estimate',source:p.source,offsetMetres:Math.round(distance(loc,p)*1000),reason:p.source==='osm-snap'?'Nearby mapped track; exact berth and line unverified':'Approximate local route geometry; exact track and berth unverified'};}
 else r={status:'unavailable',reason:'No supported railway coordinate within 250 metres'};
 r.name=loc.n;positions[key]=r;summary[r.status]=(summary[r.status]||0)+1;
}
const data={version:1,maxAutomaticOffsetMetres:250,summary,positions};
fs.writeFileSync(path.join(root,'data/berth-quality.js'),'window.BERTH_QUALITY='+JSON.stringify(data)+';\n');
fs.writeFileSync(path.join(root,'data/location-audit.json'),JSON.stringify({summary,unavailable:Object.entries(positions).filter(([,p])=>p.status==='unavailable').map(([berth,p])=>({berth,...p}))},null,2)+'\n');
console.log(summary);
