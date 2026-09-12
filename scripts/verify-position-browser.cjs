
const {chromium}=require(process.env.PLAYWRIGHT_PACKAGE),assert=require('node:assert/strict');
(async()=>{
 const browser=await chromium.launch({executablePath:'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true});
 try {
 const page=await browser.newPage({viewport:{width:1280,height:900}}),errors=[];page.on('pageerror',e=>errors.push(e.message));let ws;
 await page.routeWebSocket('ws://localhost:3001',socket=>ws=socket);
 await page.goto('http://127.0.0.1:8081/congestion.html?feed=mock',{waitUntil:'networkidle'});
 const t=(td,berth)=>({td,berth,lat:55.50173,lon:3.47534,ts:Date.now(),delay:26,delayStatus:'LATE',type:'express'});
 ws.send(JSON.stringify({type:'state',trains:{'1C77':t('EX','D172'),'1C81':t('D6','0569')}}));
 await page.locator('[data-service="1C77"]').click();
 assert.ok((await page.locator('#trackingText').textContent()).includes('Whiteball'));
 await page.waitForTimeout(1200);await page.screenshot({path:'whiteball-position-check.jpg',quality:40});
 await page.locator('[data-service="1C81"]').click();
 assert.ok((await page.locator('#trackingText').textContent()).includes('Maidenhead'));
 await page.waitForTimeout(1200);await page.screenshot({path:'maidenhead-position-check.jpg',quality:40});
 await page.goto('http://127.0.0.1:8081/live-trains.html?feed=mock',{waitUntil:'networkidle'});
 ws.send(JSON.stringify({type:'state',trains:{'1C77':t('EX','D172')}}));
 await page.waitForFunction(()=>typeof trainData!=='undefined'&&trainData['1C77']);
 assert.ok(await page.evaluate(()=>trainData['1C77'].lat<51&&trainMk.has('1C77')));
 ws.send(JSON.stringify({type:'state',trains:{'1C77':t('EX','UNKN')}}));
 await page.waitForFunction(()=>!trainMk.has('1C77'));
 assert.equal(await page.evaluate(()=>visibleTrains.has('1C77')),false);
 ws.send(JSON.stringify({type:'state',trains:{'1C77':t('EX','D172')}}));
 await page.waitForFunction(()=>trainMk.has('1C77'));
 assert.equal(await page.evaluate(()=>visibleTrains.has('1C77')),true);
 assert.equal(errors.length,0,errors.join('\n'));
 console.log('Passed: corrected Whiteball/Maidenhead map tracking; Network Twin correction before UK filter; marker removal/reappearance.');
 } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exit(1);});
