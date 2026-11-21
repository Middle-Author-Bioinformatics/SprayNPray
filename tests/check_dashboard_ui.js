/* Execute a generated offline dashboard against a minimal DOM, without dependencies.
   Usage: node tests/check_dashboard_ui.js /path/to/dashboard.html
   This checks application logic/events, not browser layout or rendering. */
'use strict';
const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const html=fs.readFileSync(process.argv[2],'utf8');
class Element {
 constructor(tag='div'){this.tag=tag;this.children=[];this.attributes={};this.events={};this.value='';this.checked=false;this.style={};this.dataset={};this._text='';this._html='';const classes=new Set();this.classList={contains:c=>classes.has(c),toggle(c,force){const add=force??!classes.has(c);if(add)classes.add(c);else classes.delete(c);return add;}};}
 set textContent(v){this._text=String(v);}get textContent(){return this._text;}
 set innerHTML(v){this._html=String(v);}get innerHTML(){return this._html;}
 setAttribute(k,v){this.attributes[k]=String(v);} appendChild(c){this.children.push(c);return c;}
 replaceChildren(...children){this.children=children;}addEventListener(k,fn){this.events[k]=fn;}
 querySelector(tag){return this.children.find(c=>c.tag===tag)||null;}
 getBBox(){return {x:0,y:0,width:1000,height:100};}getBoundingClientRect(){return {left:0,top:0};}
 getContext(){return new Proxy({},{get:(_,key)=>()=>{}});}click(){this.events.click?.({target:this});}
}
const elements=new Map([...html.matchAll(/id="([^"]+)"/g)].map(m=>{const el=new Element();el.id=m[1];return [m[1],el];}));
for(const match of html.matchAll(/<[^>]+\bid="([^"]+)"[^>]*>/g))for(const c of (match[0].match(/\bclass="([^"]+)"/)?.[1]??'').split(/\s+/).filter(Boolean))elements.get(match[1]).classList.toggle(c,true);
elements.get('statusFilter').value='all';
elements.get('payload').textContent=html.match(/<script id="payload" type="application\/json">([\s\S]*?)<\/script>/)[1];
const tabs=[...html.matchAll(/<button\b([^>]*\bdata-tab="([^"]+)"[^>]*)>/g)].map(m=>{const el=new Element('button');el.dataset.tab=m[2];for(const c of (m[1].match(/\bclass="([^"]+)"/)?.[1]??'').split(/\s+/).filter(Boolean))el.classList.toggle(c,true);return el;});
const document={getElementById:id=>{assert(elements.has(id),'Unknown DOM ID '+id);return elements.get(id);},querySelectorAll:s=>s==='[data-tab]'?tabs:s==='.panel'?tabs.map(t=>elements.get(t.dataset.tab)):[],createElementNS:(_,tag)=>new Element(tag),createElement:tag=>new Element(tag)};
const context=vm.createContext({document,window:{devicePixelRatio:1},innerWidth:1440,innerHeight:1000,URL:{createObjectURL:()=>'',revokeObjectURL(){}},Blob,XMLSerializer:class{serializeToString(){return '<svg/>';}},setTimeout:fn=>fn(),console});
const script=[...html.matchAll(/<script>([\s\S]*?)<\/script>/g)].map(m=>m[1]).join('\n');vm.runInContext(script,context);
const get=id=>elements.get(id), event=(id,k,e={})=>{assert(get(id).events[k]);get(id).events[k](e);};
assert.deepEqual(tabs.map(t=>t.dataset.tab),['assemblies','ani','phylogeny']);
assert(!elements.has('overview'));assert(!elements.has('finding'));
const activeTabs=()=>tabs.filter(t=>t.classList.contains('active')).map(t=>t.dataset.tab);
const activePanels=()=>tabs.map(t=>t.dataset.tab).filter(id=>get(id).classList.contains('active'));
assert.deepEqual(activeTabs(),['assemblies']);assert.deepEqual(activePanels(),['assemblies']);
for(const tab of [tabs[1],tabs[2],tabs[0]]){tab.click();assert.deepEqual(activeTabs(),[tab.dataset.tab]);assert.deepEqual(activePanels(),[tab.dataset.tab]);}
assert.equal(get('rowCount').textContent,'5 samples shown');
assert.equal(vm.runInContext('new Set(Object.values(groups).filter(g=>g!=="Unresolved")).size',context),2);
assert.match(get('aniRange').textContent,/97\.400/);
assert(html.indexOf('href="data/assemblies.csv"')>html.indexOf('<section id="assemblies"'));
assert(html.indexOf('href="data/assemblies.csv"')<html.indexOf('<section id="ani"'));
assert(html.indexOf('href="data/sample_ani.csv"')>html.indexOf('<section id="ani"'));
assert(html.indexOf('id="exportGroups"')>html.indexOf('<section id="ani"'));
assert(!get('assemblyBody').innerHTML.includes('<script>alert(1)</script>'));
get('sampleSearch').value='sample_flye';event('sampleSearch','input');assert.equal(get('rowCount').textContent,'1 samples shown');
event('clearFilters','click');assert.equal(get('rowCount').textContent,'5 samples shown');
get('statusFilter').value='outside';event('statusFilter','change');assert.equal(get('rowCount').textContent,'2 samples shown');
get('statusFilter').value='all';event('statusFilter','change');
get('threshold').value='95';event('applyThreshold','click');assert.equal(vm.runInContext('new Set(Object.values(groups).filter(g=>g!=="Unresolved")).size',context),2);
get('threshold').value='101';event('applyThreshold','click');assert.match(get('thresholdError').textContent,/0 to 100/);
get('threshold').value='99';event('applyThreshold','click');
const svg=()=>get('treeViewport').querySelector('svg');
assert(svg());const sampleCircles=svg().children.filter(n=>n.tag==='circle');assert.equal(sampleCircles.length,3);sampleCircles[0].click();assert.match(get('treeDetail').innerHTML,/sample_a/);assert.match(get('treeDetail').innerHTML,/98\.50/);assert(!get('treeDetail').innerHTML.includes('<script>alert(1)</script>'));
get('includeReferences').checked=true;event('includeReferences','change');assert.equal(svg().children.filter(n=>n.tag==='circle').length,4);
event('zoomIn','click');assert.equal(Number(svg().attributes.width),1375);event('resetTree','click');assert.equal(Number(svg().attributes.width),1100);
event('heatmap','click',{clientX:200,clientY:185});assert.match(get('pairDetail').textContent,/99\.90000%/);
event('exportGroups','click');event('saveTree','click');
const py=JSON.parse(get('payload').textContent).samples;const jsGroups=vm.runInContext('groups',context);for(const s of py)assert.equal(jsGroups[s.sample],s.ani_group);
console.log('Offline dashboard logic passed: filters, groups, missing ANI, click details, reference toggle, zoom and downloads. Browser layout not tested.');
