 'use strict';
(() => {
const $=id=>document.getElementById(id);
const esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const options=()=>({type:$('type').value,threshold:+$('threshold').value,age:+$('age').value});
const ageText=ts=>Number.isFinite(ts)?Math.max(0,Math.floor((Date.now()-ts)/60000))+'m ago':'Age unknown';
const canMap=RailPosition.canMap;
const identity=t=>t.id+'|'+(t.uid||'');
let raw={},received=0,selected=null,following=null,state,connected=false,retry=2000,map,circles,firstFit=false,lastFollowPosition='';
const markerMap=new Map();
const feedParam=new URLSearchParams(location.search).get('feed');
const wsUrl=feedParam==='mock'?'ws://localhost:3001':feedParam==='live'?'wss://railinsights-site-5jfi.onrender.com':['localhost','127.0.0.1',''].includes(location.hostname)?'ws://localhost:3000':'wss://railinsights-site-5jfi.onrender.com';
if(['mock','live'].includes(feedParam))document.querySelectorAll('a[href="analytics.html"]').forEach(a=>a.search='?feed='+feedParam);
function rows(){const c=state.clusters.find(c=>c.id===selected),q=$('search').value.trim().toLowerCase();return (c?c.members:state.late).filter(t=>[t.id,t.name,t.uid,t.td].some(v=>String(v||'').toLowerCase().includes(q)));}
function fit(trains){trains=trains.filter(canMap);if(map&&trains.length)map.fitBounds(trains.map(t=>[t.lat,t.lon]),{padding:[35,35],maxZoom:13});}
function follow(id){const t=state.fresh.find(t=>t.id===id);if(!t||!canMap(t))return;following=identity(t);lastFollowPosition='';render();if(map){map.setView([t.lat,t.lon],Math.max(map.getZoom(),14));if(matchMedia('(max-width:700px)').matches)$('map').scrollIntoView({behavior:'smooth',block:'center'});}}
function select(id){selected=id;following=null;render();const c=state.clusters.find(c=>c.id===id);if(c&&map)map.fitBounds(L.latLng(c.seed.lat,c.seed.lon).toBounds(6000),{padding:[25,25]});}
function initMap(){
 if(!window.L){$('map').innerHTML='<p class="empty">Map unavailable. Delayed services are listed alongside.</p>';$('basemap').disabled=true;$('resetMap').disabled=true;return;}
 map=L.map('map',{zoomAnimation:false,fadeAnimation:false}).setView([54,-2.5],6);
 map.attributionControl.addAttribution('Position geometry © <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors');
 const satellite=L.tileLayer('https://services.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',{maxZoom:19,attribution:'Imagery © Esri, Maxar, Earthstar Geographics, and the GIS User Community'}).addTo(map);
 const street=L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png',{maxZoom:19,attribution:'© <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'});
 for(const layer of [satellite,street])layer.on('tileerror',()=>{$('imageryStatus').textContent='Some map tiles unavailable · try another basemap';});
 circles=L.layerGroup().addTo(map);
 $('basemap').onchange=()=>{map.removeLayer(satellite);map.removeLayer(street);($('basemap').value==='satellite'?satellite:street).addTo(map);$('imageryStatus').textContent=$('basemap').value==='satellite'?'Satellite/aerial background · not live imagery':'OpenStreetMap · TD berth positions';};
 $('resetMap').onclick=()=>{following=null;render();fit(rows());};
}
function renderMap(tracked){if(!map)return;
 circles.clearLayers();
 for(const c of state.clusters)L.circle([c.seed.lat,c.seed.lon],{radius:3000,color:c.id===selected?'#d9ff43':'#f5b74f',weight:1,fillOpacity:.06}).addTo(circles).on('click',()=>select(c.id));
 const visible=rows().filter(canMap);if(tracked&&!visible.some(t=>identity(t)===following))visible.push(tracked);
 const ids=new Set(visible.map(identity));
 for(const [key,m]of markerMap)if(!ids.has(key)){map.removeLayer(m);markerMap.delete(key);}
 for(const t of visible){const key=identity(t),active=key===following;let m=markerMap.get(key);
  if(!m){m=L.circleMarker([t.lat,t.lon]).addTo(map).on('click',()=>follow(t.id));m.bindTooltip('',{direction:'top'});markerMap.set(key,m);}
  m.setLatLng([t.lat,t.lon]);m.setRadius(active?10:6);m.setStyle({color:active?'#d9ff43':'#10130f',weight:active?3:1,fillColor:t.delayStatus==='LATE'?(t.delay>15?'#ff776a':'#f5b74f'):'#75df99',fillOpacity:t.positionStatus==='route-estimate'?.22:1});
  m.setTooltipContent(`${esc(t.id)} · ${t.delayStatus==='LATE'?'+'+esc(t.delay)+'m':esc(t.delayStatus||'Delay unknown')} · ${esc(t.name||t.td)}<br>Delay report: ${esc(ageText(t.delayTs))}<br>${esc(t.positionLabel)}`);
  if(active){m.bringToFront();m.openTooltip();}
 }
 if(tracked){const pos=tracked.lat+','+tracked.lon;if(pos!==lastFollowPosition){map.panTo([tracked.lat,tracked.lon],{animate:false});lastFollowPosition=pos;}}
 if(!firstFit&&visible.length){fit(visible);firstFit=true;}
}
function render(){
 const stale=!received||Date.now()-received>60000;
 state=CongestionModel.analyse(raw,options());
 if(stale)state={fresh:[],known:[],reported:[],unverified:[],late:[],clusters:[],excluded:0};
 $('feed').textContent=!received?(connected?'Waiting for data':'Connecting…'):!connected?'Disconnected':stale?'Feed stale':'Receiving observations';
 $('feed').style.color=connected&&!stale?'var(--signal)':'var(--amber)';
 $('updated').textContent=received?'Updated '+new Date(received).toLocaleTimeString('en-GB'):'Waiting for feed';
 $('delayed').textContent=stale?'—':state.reported.length?state.late.length:'—';$('clusters').textContent=stale?'—':state.clusters.length;
 const unknown=state.late.filter(t=>t.delayTs==null).length;
 $('dataNotice').textContent=stale?'Waiting for current train observations.':unknown?`${unknown} delayed services have no report time — shown with age unknown.`:!state.reported.length?'Delay reports unavailable for these services.':!connected?'Connection lost — showing the last snapshot.':'';
 const unlocated=state.late.filter(t=>!canMap(t)).length;
 if(unlocated)$('dataNotice').textContent+=` ${unlocated} services have unverified locations (list only).`;
 $('dataNotice').hidden=!$('dataNotice').textContent;
 $('coverageNote').textContent=`${state.fresh.filter(canMap).length} plotted positions; ${state.known.length} fresh delay reports; ${state.unverified.length} reports with age unknown; ${state.excluded} stale or invalid positions excluded.`;
 if(selected&&!state.clusters.some(c=>c.id===selected))selected=null;
 const c=state.clusters.find(c=>c.id===selected),list=rows(),tracked=state.fresh.find(t=>identity(t)===following&&canMap(t));
 $('scopeName').textContent=c?'Cluster: '+(c.seed.name||c.seed.td):'';$('clearSelection').hidden=!c;$('scope').hidden=!c;
 $('listCount').textContent=list.length;
 $('services').innerHTML=list.length?list.map(t=>`<button class="service" data-service="${esc(t.id)}" aria-pressed="${identity(t)===following}" ${canMap(t)?'':'disabled'}><span class="service-top"><b>${esc(t.id)}</b><strong>+${esc(t.delay)}m</strong></span><span class="location">${esc(t.name||t.td||'Unknown location')}</span><small>${esc(t.positionLabel)} · ${ageText(t.ts)}</small><span class="service-bottom"><small class="${t.delayTs==null?'amber':''}">Delay: ${ageText(t.delayTs)}</small><span class="follow-label">${!canMap(t)?'Not plotted':identity(t)===following?'Following ●':'Follow on map →'}</span></span></button>`).join(''):`<p class="empty">${stale?'Waiting for train observations.':!state.reported.length?'No delay reports available.':$('search').value?'No matching service. Try its headcode or location.':'No services match these filters.'}</p>`;
 $('stopFollow').hidden=!following;
 $('trackingText').textContent=following?(tracked?`Following ${tracked.id} · ${tracked.positionLabel} · ${tracked.name||tracked.td||'Unknown location'} · Position ${ageText(tracked.ts)} · ${tracked.delayStatus==='LATE'?'Reported +'+tracked.delay+'m':'Status: '+(tracked.delayStatus||'Unknown')} · Delay ${ageText(tracked.delayTs)}`:'Follow paused — service no longer has a position in this selection.'):'Select a service to follow it on the map.';
 $('tracking').classList.toggle('active',!!tracked);
 $('clusterCount').textContent=state.clusters.length;
 $('ranking').innerHTML=state.clusters.length?state.clusters.map(c=>`<button class="cluster" data-id="${esc(c.id)}" aria-pressed="${c.id===selected}"><b>${esc(c.seed.name||c.seed.td)}</b><span>${c.members.length} delayed · mean +${c.mean.toFixed(0)}m →</span></button>`).join(''):'<p class="empty">No clusters under these filters.</p>';
 $('export').disabled=!list.length;renderMap(tracked);
}
$('services').onclick=e=>{const b=e.target.closest('[data-service]');if(b)follow(b.dataset.service);};
$('ranking').onclick=e=>{const b=e.target.closest('[data-id]');if(b)select(b.dataset.id);};
$('clearSelection').onclick=()=>{selected=null;render();};
$('stopFollow').onclick=()=>{following=null;render();};
$('search').oninput=()=>render();
for(const id of ['type','threshold','age'])$(id).onchange=()=>{selected=null;render();};
$('export').onclick=()=>{
 const cell=v=>'"'+String(v??'').replace(/^[=+@-]/,"'$&").replace(/"/g,'""')+'"';
 const lines=[['Snapshot UTC','Service','UID','Location','TD','Berth','Type','Reported lateness minutes','Position UTC','Delay report UTC','Report freshness','Position quality','Latitude','Longitude']];
 rows().forEach(t=>lines.push([new Date(received).toISOString(),t.id,t.uid,t.name,t.td,t.berth,t.type,t.delay,new Date(t.ts).toISOString(),t.delayTs==null?'':new Date(t.delayTs).toISOString(),t.delayTs==null?'Age unknown':'Within 15 minutes',t.positionLabel,t.lat,t.lon]));
 const url=URL.createObjectURL(new Blob(['\uFEFF'+lines.map(r=>r.map(cell).join(',')).join('\r\n')],{type:'text/csv;charset=utf-8'}));const a=document.createElement('a');a.href=url;a.download='rail-congestion-'+new Date().toISOString().slice(0,10)+'.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
};
function connect(){let ws;try{ws=new WebSocket(wsUrl);}catch{setTimeout(connect,retry);return;}
 ws.onopen=()=>{connected=true;render();};
 ws.onmessage=event=>{let msg;try{msg=JSON.parse(event.data);}catch{return;}if(msg.type!=='state'||!msg.trains||typeof msg.trains!=='object'||Array.isArray(msg.trains))return;raw=Object.fromEntries(Object.entries(msg.trains).filter(([,t])=>t&&typeof t==='object').map(([id,t])=>[id,RailPosition.resolve(t,window.BERTH_QUALITY)]));received=Date.now();retry=2000;render();};
 ws.onerror=()=>ws.close();ws.onclose=()=>{connected=false;render();setTimeout(connect,retry);retry=Math.min(30000,retry*1.5);};
}
initMap();render();connect();setInterval(render,10000);
})();
