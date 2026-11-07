/* Ejecutar con: node tests/math-tests.js */
global.window=global;
require('../js/engineering-core.js');
require('../js/fields-bolivia.js');
const E=global.IRMAEngine;
const tests=[];
function test(name,fn){try{fn();tests.push({name,status:'PASS'});}catch(e){tests.push({name,status:'FAIL',error:e.message});}}
const near=(a,b,tol=1e-7)=>{if(Math.abs(a-b)>tol)throw new Error(`${a} != ${b}`)};
test('Vogel: q=qmax cuando Pwf=0',()=>near(E.vogel({Pr:3000,qmax:1000,Pwf:0}).q,1000));
test('Vogel: q=0 cuando Pwf=Pr',()=>near(E.vogel({Pr:3000,qmax:1000,Pwf:3000}).q,0));
test('Vogel rechaza Pwf>Pr',()=>{let ok=false;try{E.vogel({Pr:3000,qmax:1000,Pwf:3100})}catch(e){ok=true}if(!ok)throw Error('No rechazó Pwf>Pr')});
test('Minimum curvature: pozo vertical',()=>near(E.minimumCurvature([{md:0,inc:0,azi:0},{md:100,inc:0,azi:0}]).final.tvd,100));
test('Minimum curvature: horizontal hacia Este',()=>near(E.minimumCurvature([{md:0,inc:90,azi:90},{md:100,inc:90,azi:90}]).final.east,100));
test('Minimum curvature: DLS nulo para línea recta',()=>near(E.minimumCurvature([{md:0,inc:30,azi:45},{md:100,inc:30,azi:45}]).final.dls,0));
test('PLT: contribuciones líquidas suman 100%',()=>{const r=E.plt([{top:100,base:120,qo:10,qw:5,qg:1},{top:130,base:150,qo:20,qw:5,qg:2}]);near(r.rows.reduce((s,x)=>s+x.liquidContribution,0),100)});
test('PLT: no acepta base menor que tope',()=>{let ok=false;try{E.plt([{top:200,base:100,qo:1,qw:1,qg:0}])}catch(e){ok=true}if(!ok)throw Error('No rechazó intervalo invertido')});
test('Duong: declinación con a positiva',()=>{const r=E.duong({q1:100,m:.7,a:.01,days:365,abandonmentRate:5});if(!(r.finalRate<100&&r.cumulative>0))throw Error('Pronóstico no decreciente');});
test('Corey: fw en extremos 0 y 1',()=>{const r=E.fractionalFlow({swi:.2,sor:.2,krw0:.3,kro0:.9,nw:2,no:2,muw:1,muo:3});near(r.rows[0].fw,0);near(r.rows.at(-1).fw,1)});
test('Gradiente monofásico: caída hidrostática positiva',()=>{const r=E.verticalLiquid({depth:1800,diameter:62,roughness:.045,rateBpd:500,density:850,viscosity:4,whp:20});if(!(r.hydrostaticBar>0&&r.frictionFactor>0&&r.bottomholePressure>20))throw Error('Gradiente inválido');if(r.frictionFactor>.2)throw Error('Factor Darcy inesperadamente alto');});
test('Gradiente monofásico: rechaza régimen transicional',()=>{let ok=false;try{E.verticalLiquid({depth:1800,diameter:62,roughness:.045,rateBpd:350,density:850,viscosity:4,whp:20})}catch(e){ok=/transicional/.test(e.message)}if(!ok)throw Error('No rechazó régimen transicional')});
test('Peng-Robinson: Z de metano dentro de rango físico',()=>{const r=E.pengRobinson({T:300,P:5000,composition:[{component:'methane',x:1}]});if(!(r.vapor.Z>r.B&&r.vapor.Z<2))throw Error('Z no físico: '+r.vapor.Z)});
test('Peng-Robinson: mezcla multicomponente produce raíz física',()=>{const r=E.pengRobinson({T:350,P:12000,composition:[{component:'methane',x:.78},{component:'ethane',x:.08},{component:'propane',x:.04},{component:'co2',x:.02},{component:'nitrogen',x:.08}]});if(!(r.vapor.Z>r.B&&r.vapor.Z<2))throw Error('Z mezcla inválido')});
test('Peng-Robinson: composición incorrecta se rechaza',()=>{let ok=false;try{E.pengRobinson({T:300,P:5000,composition:[{component:'methane',x:.7}]})}catch(e){ok=true}if(!ok)throw Error('No rechazó suma molar distinta de 1')});
test('Catálogo: contiene 52 campos',()=>{if(IRMAFields.all().length!==52)throw Error('Conteo '+IRMAFields.all().length)});
test('Catálogo: búsqueda por campo',()=>{if(IRMAFields.search('Margarita').length!==1)throw Error('Búsqueda no devuelve Margarita')});
console.table(tests);const failed=tests.filter(t=>t.status==='FAIL');console.log(`${tests.length-failed.length}/${tests.length} pruebas superadas`);if(failed.length)process.exit(1);
