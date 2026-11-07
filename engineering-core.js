/* IRMA-X Engineering Core — modelos explícitos, unidades SI/campo y validaciones.
   Los modelos empíricos solo deben usarse dentro de su dominio de aplicación. */
(function (g) {
  'use strict';
  const finite = (x, n) => Number.isFinite(Number(x)) && (n === undefined || Number(x) >= n);
  const assert = (ok, msg) => { if (!ok) throw new Error(msg); };
  const rad = d => d * Math.PI / 180;
  const clamp = (x,a,b) => Math.min(b,Math.max(a,x));

  function vogel({Pr, qmax, Pwf}) {
    [Pr,qmax,Pwf].forEach((v,i)=>assert(finite(v,0) && (i!==1 || Number(v)>0), 'Vogel: Pr y Pwf deben ser ≥ 0 y qmax > 0.'));
    assert(Pwf <= Pr && Pr > 0, 'Vogel requiere 0 ≤ Pwf ≤ Pr y Pr > 0.');
    const r=Pwf/Pr, q=qmax*(1-0.2*r-0.8*r*r);
    return {q:Math.max(0,q), equation:'q/qmax = 1 − 0.2(Pwf/Pr) − 0.8(Pwf/Pr)²', ratio:r, model:'Vogel (petróleo saturado)', warning:'Aplicable a petróleo bajo condiciones compatibles con Vogel; no usar para gas ni como sustituto de una prueba calibrada.'};
  }
  function vogelCurve(Pr,qmax,n=41){return Array.from({length:n},(_,i)=>{const Pwf=Pr*i/(n-1);return {Pwf,...vogel({Pr,qmax,Pwf})};});}

  // Minimum Curvature Method. MD in m, inclination/azimuth in degrees; N/E/TVD in m.
  function minimumCurvature(stations, origin={north:0,east:0,tvd:0}, dlsUnit='deg/30m') {
    assert(Array.isArray(stations)&&stations.length>=2,'Se requieren al menos dos estaciones de survey.');
    let N=Number(origin.north||0),E=Number(origin.east||0),T=Number(origin.tvd||0);
    let prev=null, length=0;
    const out=[];
    stations.forEach((s,i)=>{
      const md=Number(s.md), inc=Number(s.inc), azi=Number(s.azi);
      assert([md,inc,azi].every(Number.isFinite),'Survey contiene valores no numéricos.');
      assert(inc>=0&&inc<=180,'Inclinación debe estar entre 0 y 180°.');
      assert(md>=0 && (!prev || md>prev.md),'MD debe ser positiva y estrictamente creciente.');
      if(!prev){out.push({station:i+1,md,inc,azi,north:N,east:E,tvd:T,dls:0,dogleg:0,closure:Math.hypot(N,E)});prev={md,inc,azi};return;}
      const dMD=md-prev.md, i1=rad(prev.inc),i2=rad(inc),a1=rad(prev.azi),a2=rad(azi);
      const cosDL=clamp(Math.cos(i2-i1)-Math.sin(i1)*Math.sin(i2)*(1-Math.cos(a2-a1)),-1,1);
      const dogleg=Math.acos(cosDL), rf=dogleg<1e-10?1:(2/dogleg)*Math.tan(dogleg/2);
      const dN=dMD/2*(Math.sin(i1)*Math.cos(a1)+Math.sin(i2)*Math.cos(a2))*rf;
      const dE=dMD/2*(Math.sin(i1)*Math.sin(a1)+Math.sin(i2)*Math.sin(a2))*rf;
      const dT=dMD/2*(Math.cos(i1)+Math.cos(i2))*rf;
      N+=dN;E+=dE;T+=dT;length+=dMD;
      const dls=dMD>0?dogleg*180/Math.PI*(dlsUnit==='deg/100ft'?30.48:dlsUnit==='deg/30m'?30:1)/dMD:0;
      out.push({station:i+1,md,inc,azi,north:N,east:E,tvd:T,dN,dE,dT,dogleg:dogleg*180/Math.PI,dls,closure:Math.hypot(N,E)});
      prev={md,inc,azi};
    });
    return {stations:out,totalMeasuredDepth:stations.at(-1).md,trajectoryLength:length,final:out.at(-1),dlsUnit,equation:'RF = (2/β) tan(β/2); ΔN = ΔMD/2(sin I₁ cos A₁ + sin I₂ cos A₂)RF; análogo para E; ΔTVD = ΔMD/2(cos I₁ + cos I₂)RF'};
  }

  function plt(intervals, source='measured') {
    assert(Array.isArray(intervals)&&intervals.length>0,'Añade al menos un intervalo.');
    const rows=intervals.map((x,i)=>{
      const top=Number(x.top),base=Number(x.base),qo=Number(x.qo||0),qw=Number(x.qw||0),qg=Number(x.qg||0);
      assert([top,base,qo,qw,qg].every(Number.isFinite)&&top>=0&&base>top&&qo>=0&&qw>=0&&qg>=0,`Intervalo ${i+1}: revisa profundidades y caudales no negativos.`);
      const qt=qo+qw, watercut=qt>0?100*qw/qt:0;
      return {name:x.name||`Intervalo ${i+1}`,top,base,qo,qw,qg,qliq:qt,watercut};
    });
    const totals=rows.reduce((a,r)=>({qo:a.qo+r.qo,qw:a.qw+r.qw,qg:a.qg+r.qg,qliq:a.qliq+r.qliq}),{qo:0,qw:0,qg:0,qliq:0});
    rows.forEach(r=>{r.liquidContribution=totals.qliq?100*r.qliq/totals.qliq:0;r.oilContribution=totals.qo?100*r.qo/totals.qo:0;r.gasContribution=totals.qg?100*r.qg/totals.qg:0;r.flag=r.qliq===0?'Sin aporte líquido':r.watercut>=80?'Alto corte de agua':r.liquidContribution<5?'Baja contribución':'Aporte activo';});
    return {rows,totals,source,classification:source==='measured'?'Datos ingresados como mediciones; verificar contra calibración y condiciones de prueba.':'Estimación sintética; no constituye un registro PLT.',equation:'q_acumulado(z) = Σ q_intervalo; contribuciónᵢ = qᵢ / Σq × 100'};
  }

  // Duong: t is elapsed producing time in days; q1 is rate at t=1 day.
  function duong({q1,m,a,days,abandonmentRate=0,points=100}) {
    assert(finite(q1,0)&&Number(q1)>0&&finite(m)&&Number(m)>0&&Number(m)<1&&finite(a)&&Number(a)>=0,'Duong: q1 > 0, 0 < m < 1 y a ≥ 0 (a se define como magnitud positiva de declinación).');
    assert(finite(days,1)&&Number(days)>=1,'Horizonte de pronóstico ≥ 1 día.');
    const n=Math.max(2,Math.min(1000,Math.floor(points))), out=[];
    const rate=t=>q1*Math.pow(t,-m)*Math.exp((-a/(1-m))*(Math.pow(t,1-m)-1));
    let cum=0,prevT=1,prevQ=rate(1);
    for(let i=0;i<n;i++){const t=1+(days-1)*i/(n-1),q=rate(t);if(i)cum+=(prevQ+q)*0.5*(t-prevT);out.push({day:t,rate:q,cumulative:cum});prevT=t;prevQ=q;}
    let eur=null;if(Number(abandonmentRate)>0){let lo=1,hi=days;if(rate(hi)>abandonmentRate)eur=null;else{for(let k=0;k<60;k++){let mid=(lo+hi)/2;if(rate(mid)>abandonmentRate)lo=mid;else hi=mid;}if(rate(hi)<=abandonmentRate){const tAb=hi,N=120,dt=(tAb-1)/N;let acc=0,pt=rate(1);for(let j=1;j<=N;j++){const qt=rate(1+j*dt);acc+=(pt+qt)*0.5*dt;pt=qt;}eur=acc;}else eur=null;}}
    return {points:out,qAtOneDay:q1,finalRate:out.at(-1).rate,cumulative:cum,eur,model:'Duong',equation:'q(t)=q₁ t^(−m) exp[−a/(1−m)(t^(1−m)−1)] (a ≥ 0)',warning:'Modelo de flujo transitorio; requiere ajuste a historia real, tiempo consistente y validación del régimen de flujo. EUR no se declara si no se alcanza el límite de abandono.'};
  }

  function fractionalFlow({swi,sor,krw0,kro0,nw,no,muw,muo,points=51}) {
    assert([swi,sor,krw0,kro0,nw,no,muw,muo].every(Number.isFinite),'Completa todos los parámetros Corey.');
    assert(swi>=0&&sor>=0&&swi+sor<1&&krw0>=0&&kro0>=0&&nw>0&&no>0&&muw>0&&muo>0,'Parámetros Corey/viscosidades fuera de rango.');
    const denom=1-swi-sor,n=Math.max(2,Math.min(201,Math.floor(points)));
    const rows=Array.from({length:n},(_,i)=>{const sw=swi+denom*i/(n-1),swn=(sw-swi)/denom;const krw=krw0*Math.pow(swn,nw),kro=kro0*Math.pow(1-swn,no),lw=krw/muw,lo=kro/muo,fw=lw+lo?lw/(lw+lo):0;return {sw,swn,krw,kro,lambdaW:lw,lambdaO:lo,lambdaT:lw+lo,fw,fo:1-fw};});
    return {rows,equation:'f_w = (k_rw/μ_w) / [(k_rw/μ_w)+(k_ro/μ_o)]; k_rw=k_rw0(Sw*)^nw; k_ro=k_ro0(1−Sw*)^no',warning:'Modelo Corey idealizado; requiere curvas de permeabilidad relativa o calibración de exponentes y endpoints.'};
  }

  // Single-phase liquid vertical pressure drop (Darcy-Weisbach), deliberately NOT a multiphase VLP correlation.
  function verticalLiquid({depth,diameter,roughness,rateBpd,density,viscosity,whp}) {
    assert([depth,diameter,roughness,rateBpd,density,viscosity,whp].every(Number.isFinite),'Completa parámetros VLP monofásico.');
    assert(depth>0&&diameter>0&&roughness>=0&&rateBpd>=0&&density>0&&viscosity>0&&whp>=0,'Parámetros VLP fuera de rango.');
    const Q=rateBpd*0.158987294928/86400,D=diameter/1000,A=Math.PI*D*D/4,v=Q/A,Re=density*v*D/(viscosity/1000);
    let f,regime;if(Re<1e-12){f=0;regime='sin flujo';}else if(Re<2300){f=64/Re;regime='laminar';}else if(Re<4000){throw new Error('Reynolds en zona transicional (2300–4000): el factor de fricción no se calcula automáticamente por incertidumbre del régimen. Ajusta caudal/diámetro o aporta una correlación validada.');}else{const rr=(roughness/1000)/D;f=0.25/Math.pow(Math.log10(rr/3.7+5.74/Math.pow(Re,0.9)),2);for(let k=0;k<30;k++){const inv=-2*Math.log10(rr/3.7+2.51/(Re*Math.sqrt(f)));const next=1/(inv*inv);if(Math.abs(next-f)<1e-10){f=next;break;}f=next;}regime='turbulento (Colebrook iterativo)';}
    const dpHyd=density*9.80665*depth/1e5,dpFric=f*(depth/D)*(density*v*v/2)/1e5,dpAcc=0,dp=dpHyd+dpFric+dpAcc;
    return {bottomholePressure:whp+dp,velocity:v,reynolds:Re,frictionFactor:f,flowRegime:regime,hydrostaticBar:dpHyd,frictionBar:dpFric,accelerationBar:dpAcc,totalDropBar:dp,equation:'ΔP = ρgL + f(L/D)(ρv²/2); flujo monofásico líquido, vertical, propiedades constantes.',warning:'NO es Beggs & Brill ni Hagedorn & Brown. No representa flujo multifásico; no usar como VLP de gas/aceite/agua.'};
  }

  // Peng–Robinson EOS, conventional quadratic mixing rule and phi coefficients.
  const COMPONENTS={
    methane:{name:'Metano',Tc:190.564,Pc:4599.2,omega:0.01142,M:16.043},ethane:{name:'Etano',Tc:305.32,Pc:4872.2,omega:0.0995,M:30.070},propane:{name:'Propano',Tc:369.83,Pc:4248.0,omega:0.1523,M:44.097},
    isobutane:{name:'i-Butano',Tc:407.81,Pc:3629.0,omega:0.184, M:58.124},nButane:{name:'n-Butano',Tc:425.12,Pc:3796.0,omega:0.2002,M:58.124},
    nitrogen:{name:'Nitrógeno',Tc:126.192,Pc:3395.8,omega:0.0372,M:28.014},co2:{name:'CO₂',Tc:304.1282,Pc:7377.3,omega:0.22394,M:44.010},h2s:{name:'H₂S',Tc:373.10,Pc:8960.0,omega:0.1005,M:34.081},
    nPentane:{name:'n-Pentano',Tc:469.7,Pc:3370.0,omega:0.251,M:72.151},nHexane:{name:'n-Hexano',Tc:507.6,Pc:3025.0,omega:0.301,M:86.178}
  };
  function cubicRealRoots(a,b,c){ // x^3 + a x^2 + b x + c
    const p=b-a*a/3,q=2*a*a*a/27-a*b/3+c,disc=q*q/4+p*p*p/27,shift=-a/3;
    if(disc>1e-14){const u=Math.cbrt(-q/2+Math.sqrt(disc)),v=Math.cbrt(-q/2-Math.sqrt(disc));return [u+v+shift];}
    if(Math.abs(p)<1e-14)return [shift];
    const r=2*Math.sqrt(-p/3),theta=Math.acos(clamp((3*q/(2*p))*Math.sqrt(-3/p),-1,1))/3;
    return [r*Math.cos(theta)+shift,r*Math.cos(theta-2*Math.PI/3)+shift,r*Math.cos(theta-4*Math.PI/3)+shift].sort((x,y)=>x-y);
  }
  function pengRobinson({T,P,composition,kij={}}){
    assert(finite(T,0)&&finite(P,0)&&T>0&&P>0,'Peng-Robinson: T y P deben ser positivos en K y kPa.');
    assert(Array.isArray(composition)&&composition.length>0,'Introduce composición molar.');
    const sum=composition.reduce((s,c)=>s+Number(c.x),0);assert(Math.abs(sum-1)<1e-6,'Las fracciones molares deben sumar 1.');
    const comps=composition.map(c=>{const p=COMPONENTS[c.component];assert(p,`Componente no soportado: ${c.component}`);assert(finite(c.x,0),'Fracción molar inválida.');const Tr=T/p.Tc,k=0.37464+1.54226*p.omega-0.26992*p.omega*p.omega,alpha=Math.pow(1+k*(1-Math.sqrt(Tr)),2),a=0.45724* Math.pow(8.314462618,2)*p.Tc*p.Tc/(p.Pc*1000)*alpha,b=0.07780*8.314462618*p.Tc/(p.Pc*1000);return {...p,x:Number(c.x),a,b,alpha};});
    const R=8.314462618,am=comps.reduce((s,ci,i)=>s+comps.reduce((z,cj,j)=>z+ci.x*cj.x*Math.sqrt(ci.a*cj.a)*(1-Number(kij[`${composition[i].component}|${composition[j].component}`]||0)),0),0),bm=comps.reduce((s,c)=>s+c.x*c.b,0);
    const PPa=P*1000,A=am*PPa/Math.pow(R*T,2),B=bm*PPa/(R*T),roots=cubicRealRoots(-(1-B),A-3*B*B-2*B,-(A*B-B*B-B*B*B));
    const valid=roots.filter(z=>z>B&&Number.isFinite(z));assert(valid.length,'Peng-Robinson: no se obtuvo raíz física Z > B.');
    const phiFor=Z=>comps.map((ci,i)=>{const sumA=comps.reduce((s,cj,j)=>s+comps[j].x*Math.sqrt(ci.a*cj.a)*(1-Number(kij[`${composition[i].component}|${composition[j].component}`]||0)),0);const bi=ci.b/bm,term=2*sumA/am-bi;const L=Math.log((Z+(1+Math.SQRT2)*B)/(Z+(1-Math.SQRT2)*B));const lnphi=bi*(Z-1)-Math.log(Z-B)-A/(2*Math.SQRT2*B)*term*L;return {component:composition[i].component,phi:Math.exp(lnphi),lnPhi:lnphi};});
    const vaporZ=valid.at(-1),liquidZ=valid[0],MW=comps.reduce((s,c)=>s+c.x*c.M,0),rhoV=PPa*MW/(R*T*vaporZ)/1000;
    return {A,B,roots:valid,vapor:{Z:vaporZ,phi:phiFor(vaporZ)},liquid:{Z:liquidZ,phi:phiFor(liquidZ)},molecularWeight:MW, vaporDensityKgM3:rhoV,equation:'P = RT/(v−b) − aα/[v(v+b)+b(v−b)]',warning:'PR clásico con regla cuadrática de mezcla; kij no suministrado se toma 0. Densidad reportada solo para vapor idealizado por EOS. Equilibrio flash, propiedades líquidas y caracterización de fracciones pesadas requieren datos y algoritmos adicionales.'};
  }

  g.IRMAEngine={vogel,vogelCurve,minimumCurvature,plt,duong,fractionalFlow,verticalLiquid,pengRobinson,COMPONENTS};
})(window);
