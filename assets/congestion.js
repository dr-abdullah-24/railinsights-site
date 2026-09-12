'use strict';
(() => {
const $=id=>document.getElementById(id), esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const options=()=>({type:$('type').value,threshold:+$('threshold').value,age:+$('age').value});
let raw={},received=0,selected=null,state,connected=false,retry=2000,map,markers,circles;
const baselines=new Map(), feedParam=new URLSearchParams(location.search).get('feed');
const wsUrl=feedParam==='mock'?'ws://localhost:3001':feedParam==='live'?'wss://railinsights-site-5jfi.onrender.com':['localhost','127.0.0.1',''].includes(location.hostname)?'ws://localhost:3000':'wss://railinsights-site-5jfi.onrender.com';
if(['mock','live'].includes(feedParam))document.querySelectorAll('a[href="analytics.html"]').forEach(a=>a.search='?feed='+feedParam);
const ageText=ts=>Number.isFinite(ts)?Math.max(0,Math.floor((Date.now()-ts)/60000))+' min':'Unknown';
const identity=t=>t.id+'|'+(t.uid||'');
function change(t){const b=baselines.get(identity(t));if(!b||t.delayTs<=b.ts)return '—';const d=t.delay-b.delay;return(d>0?'+':'')+d+' min';}
function fit(trains){if(map&&trains.length)map.fitBounds(trains.map(t=>[t.lat,t.lon]),{padding:[35,35],maxZoom:13});}
function select(id){selected=id;render();const c=state.clusters.find(c=>c.id===id);if(c&&map)map.fitBounds(L.latLng(c.seed.lat,c.seed.lon).toBounds(6000),{padding:[25,25]});}
function initMap(){
 if(!window.L){$('map').innerHTML='<p class="empty">Map library unavailable. Cluster analysis and service details remain available.</p>';$('basemap').disabled=true;$('resetMap').disabled=true;return;}
 map=L.map('map',{zoomAnimation:false,fadeAnimation:false}).setView([54,-2.5],6);
 const satellite=L.tileLayer('https://services.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',{maxZoom:19,attribution:'Imagery © Esri, Maxar, Earthstar Geographics, and the GIS User Community'}).addTo(map);
 const street=L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png',{maxZoom:19,attribution:'© <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'});
 for(const layer of [satellite,street])layer.on('tileerror',()=>{$('imageryStatus').textContent='Some map tiles could not load. Try switching basemap. Service observations remain available.';});
 markers=L.layerGroup().addTo(map);circles=L.layerGroup().addTo(map);
 $('basemap').onchange=()=>{map.removeLayer(satellite);map.removeLayer(street);($('basemap').value==='satellite'?satellite:street).addTo(map);$('imageryStatus').textContent=$('basemap').value==='satellite'?'Esri satellite/aerial mosaic · capture date varies and is not supplied here · not a live satellite feed.':'OpenStreetMap geographic context · train locations are TD berth observations.';};
 $('resetMap').onclick=()=>fit(state?.late||[]);
}
function renderMap(){if(!map)return;markers.clearLayers();circles.clearLayers();
 for(const c of state.clusters)L.circle([c.seed.lat,c.seed.lon],{radius:3000,color:c.id===selected?'#d9ff43':'#f5b74f',weight:2,fillOpacity:.1}).addTo(circles).on('click',()=>select(c.id));
 for(const t of state.late)L.circleMarker([t.lat,t.lon],{radius:6,color:'#10130f',weight:1,fillColor:t.delay>15?'#ff776a':'#f5b74f',fillOpacity:.95}).bindPopup(`<b>${esc(t.id)} · +${esc(t.delay)} min</b><br>${esc(t.name||t.td||'Unknown location')}<br>Position: ${esc(ageText(t.ts))} ago<br>Delay report: ${esc(ageText(t.delayTs))} ago<br>TD berth estimate · not GPS`).addTo(markers);
}
function render(){
 const now=Date.now(),stale=!received||now-received>60000;
 state=CongestionModel.analyse(raw,options(),now);
 if(stale)state={fresh:[],known:[],late:[],clusters:[],excluded:0};
 $('feed').textContent=!received?(connected?'Waiting for data':'Connecting / reconnecting…'):!connected?'Disconnected · last snapshot '+ageText(received)+' ago':stale?'Feed stale · observations withheld':'Receiving observations';
 $('feed').style.color=connected&&!stale?'var(--signal)':'var(--amber)';
 $('updated').textContent=received?'Snapshot '+new Date(received).toLocaleTimeString('en-GB')+' · refreshes automatically':'Waiting for a feed snapshot';
 for(const [id,value]of Object.entries({positions:state.fresh.length,coverage:state.fresh.length?Math.round(100*state.known.length/state.fresh.length)+'%':'—',delayed:state.late.length,clusters:state.clusters.length}))$(id).textContent=stale?'—':value;
 $('excluded').textContent=state.excluded+' stale / invalid positions excluded';$('thresholdNote').textContent='More than '+options().threshold+' minutes late';
 if(selected&&!state.clusters.some(c=>c.id===selected))selected=null;
 const missingTime=!stale&&state.fresh.some(t=>typeof t.delay==='number'&&!t.delayTs);
 $('ranking').className=state.clusters.length?'':'empty';
 $('ranking').innerHTML=state.clusters.length?state.clusters.map((c,i)=>`<button class="cluster" data-id="${esc(c.id)}" aria-pressed="${c.id===selected}"><small>${String(i+1).padStart(2,'0')} / ${c.severe?'Severe delays present':'Delay concentration'}</small><strong>${esc(c.seed.name||c.seed.td||'Unnamed location')}</strong><div class="stats"><span>${c.members.length} services</span><span>${c.total.toFixed(0)} min summed</span></div><small>Mean ${c.mean.toFixed(1)} min · ${c.severe} over 15 min</small></button>`).join(''):stale?'Waiting for a current feed snapshot. No sample services are shown.':missingTime?'Delay timestamps unavailable. Deploy the updated bridge to enable trustworthy delay analysis.':!state.known.length?'No fresh delay reports in this selection. This does not mean the railway is running on time.':'No groups of three delayed services within 3 km under these filters.';
 const cluster=state.clusters.find(c=>c.id===selected),rows=cluster?cluster.members:state.late;
 $('detailTitle').textContent=cluster?(cluster.seed.name||cluster.seed.td||'Selected cluster'):'Affected services';
 $('detailSummary').textContent=cluster?`${rows.length} delayed services · ${cluster.total.toFixed(0)} minutes summed current lateness · mean ${cluster.mean.toFixed(1)} min.`:'All delayed services matching the current filters. Select a concentration to narrow the investigation.';
 $('clearSelection').hidden=!cluster;
 $('assessment').textContent=cluster?`${cluster.severe} of ${rows.length} services are over 15 minutes late. Check shared route and platform constraints before interpreting this geographic cluster as a congestion bottleneck.`:'A cluster highlights co-located lateness. It does not establish the cause of delay.';
 $('services').innerHTML=rows.length?rows.map(t=>`<tr><td><b>${esc(t.id)}</b><small>${esc(t.uid||'UID unavailable')}</small></td><td>${esc(t.name||t.td||'Unknown')}<small>${esc(t.td)} / ${esc(t.berth)}</small></td><td>${esc(t.type||'Unknown')}</td><td>+${esc(t.delay)} min</td><td>${ageText(t.ts)}</td><td>${ageText(t.delayTs)}</td><td>${esc(change(t))}</td></tr>`).join(''):'<tr><td colspan="7">'+(stale?'Waiting for a current feed snapshot.':'No delayed services with fresh delay reports match these filters.')+'</td></tr>';
 $('export').disabled=!rows.length;renderMap();
}
$('ranking').onclick=e=>{const b=e.target.closest('[data-id]');if(b)select(b.dataset.id);};
$('clearSelection').onclick=()=>{selected=null;render();};
for(const id of ['type','threshold','age'])$(id).onchange=()=>{selected=null;render();};
$('export').onclick=()=>{
 const cluster=state.clusters.find(c=>c.id===selected),rows=cluster?cluster.members:state.late;
 const cell=v=>'"'+String(v??'').replace(/^[=+@-]/,"'$&").replace(/"/g,'""')+'"';
 const lines=[['Snapshot UTC','Service','UID','Location','TD','Berth','Type','Lateness minutes','Position UTC','Delay report UTC','Session change','Threshold minutes','Max position age minutes']];
 rows.forEach(t=>lines.push([new Date(received).toISOString(),t.id,t.uid,t.name,t.td,t.berth,t.type,t.delay,new Date(t.ts).toISOString(),new Date(t.delayTs).toISOString(),change(t),options().threshold,options().age]));
 const url=URL.createObjectURL(new Blob(['\uFEFF'+lines.map(r=>r.map(cell).join(',')).join('\r\n')],{type:'text/csv;charset=utf-8'}));const a=document.createElement('a');a.href=url;a.download='rail-congestion-'+new Date().toISOString().slice(0,10)+'.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
};
function connect(){let ws;try{ws=new WebSocket(wsUrl);}catch{setTimeout(connect,retry);return;}
 ws.onopen=()=>{connected=true;render();};
 ws.onmessage=event=>{let msg;try{msg=JSON.parse(event.data);}catch{return;}if(msg.type!=='state'||!msg.trains||typeof msg.trains!=='object'||Array.isArray(msg.trains))return;
 raw=msg.trains;received=Date.now();retry=2000;const now=Date.now();
 for(const [id,t]of Object.entries(raw)){if(!t||typeof t!=='object')continue;const key=identity({...t,id});if(!baselines.has(key)&&Number.isFinite(t.delay)&&t.delayStatus==='LATE'&&Number.isFinite(t.delayTs)&&now-t.delayTs<=900000&&t.delayTs<=now+30000)baselines.set(key,{delay:t.delay,ts:t.delayTs});}
 for(const [key,b]of baselines)if(now-b.ts>3600000)baselines.delete(key);render();};
 ws.onerror=()=>ws.close();ws.onclose=()=>{connected=false;render();setTimeout(connect,retry);retry=Math.min(30000,retry*1.5);};
}
initMap();render();connect();setInterval(render,10000);
})();
