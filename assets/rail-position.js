(function(root){
 'use strict';
 const labels={'track-estimate':'Mapped track estimate','route-estimate':'Approximate route position',unavailable:'Location unverified'};
 function resolve(train,data){
  const key=train.td+':'+train.berth,p=data?.positions?.[key];
  if(!p||p.status==='unavailable'||!Number.isFinite(p.lat)||!Number.isFinite(p.lon))return {...train,lat:null,lon:null,positionStatus:'unavailable',positionLabel:labels.unavailable,positionReason:p?.reason||'Berth has no checked coordinate'};
  return {...train,lat:p.lat,lon:p.lon,name:p.name||train.name,positionStatus:p.status,positionLabel:p.label||labels[p.status]||labels.unavailable,positionSource:p.source,positionReason:p.reason};
 }
 function canMap(t){return t?.positionStatus!=='unavailable'&&Number.isFinite(t?.lat)&&Number.isFinite(t?.lon);}
 const api={resolve,canMap};if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.RailPosition=api;
})(typeof window!=='undefined'?window:globalThis);
