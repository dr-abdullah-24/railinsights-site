(function(root) {
  'use strict';
  const validTime = (value, now, age) => typeof value === 'number' && Number.isFinite(value) && value > 0 && value <= now + 30000 && now - value <= age;
  function distance(a,b) {
    const r=Math.PI/180, dlat=(b.lat-a.lat)*r, dlon=(b.lon-a.lon)*r;
    const h=Math.sin(dlat/2)**2+Math.cos(a.lat*r)*Math.cos(b.lat*r)*Math.sin(dlon/2)**2;
    return 6371*2*Math.asin(Math.sqrt(Math.min(1,h)));
  }
  function analyse(raw, options, now=Date.now()) {
    const all=Object.entries(raw).filter(([,t])=>t && typeof t==='object' && !t.cancelled).map(([id,t])=>({...t,id})).filter(t=> options.type==='all' || (options.type==='passenger' ? ['express','passenger'].includes(t.type) : options.type==='freight' ? t.type==='freight' : !['express','passenger','freight'].includes(t.type)));
    const fresh=all.filter(t=>typeof t.lat==='number' && typeof t.lon==='number' && Number.isFinite(t.lat) && Number.isFinite(t.lon) && Math.abs(t.lat)<=90 && Math.abs(t.lon)<=180 && validTime(t.ts,now,options.age*60000));
    const known=fresh.filter(t=>typeof t.delay==='number' && Number.isFinite(t.delay) && ['LATE','EARLY','ON TIME'].includes(t.delayStatus) && validTime(t.delayTs,now,900000));
    const reported=fresh.filter(t=>typeof t.delay==='number' && Number.isFinite(t.delay) && ['LATE','EARLY','ON TIME'].includes(t.delayStatus) && (t.delayTs == null || validTime(t.delayTs,now,900000)));
    const unverified=reported.filter(t=>t.delayTs == null);
    const late=reported.filter(t=>t.delayStatus==='LATE' && t.delay>options.threshold).sort((a,b)=>b.delay-a.delay || a.id.localeCompare(b.id));
    const assigned=new Set(), clusters=[];
    for(const seed of late) {
      if(assigned.has(seed.id)) continue;
      const members=late.filter(t=>!assigned.has(t.id) && distance(seed,t)<=3);
      if(members.length<3) continue;
      members.forEach(t=>assigned.add(t.id));
      const total=members.reduce((sum,t)=>sum+t.delay,0);
      clusters.push({id:seed.id,seed,members,total,mean:total/members.length,severe:members.filter(t=>t.delay>15).length});
    }
    clusters.sort((a,b)=>b.total-a.total || a.id.localeCompare(b.id));
    return {fresh,known,reported,unverified,late,clusters,excluded:all.length-fresh.length};
  }
  const api={analyse,distance};
  if(typeof module!=='undefined' && module.exports) module.exports=api;
  else root.CongestionModel=api;
})(typeof window!=='undefined'?window:globalThis);
