const $ = id => document.getElementById(id);

const els = {
  coatTime:$('coatTime'), beadCharge:$('beadCharge'), cellCharge:$('cellCharge'),
  mbConc:$('mbConc'), agarPct:$('agarPct'), stiffness:$('stiffness'),
  optStiffness:$('optStiffness'), agarProxy:$('agarProxy'), d0:$('d0'), ed:$('ed'),
  kd:$('kd'), hillN:$('hillN'), stiffWidth:$('stiffWidth'), obsTime:$('obsTime'),
  runBtn:$('runBtn'), sweepBtn:$('sweepBtn'), presetBtn:$('presetBtn'),
  canvas:$('simCanvas'), chart:$('chartCanvas'), status:$('status'),
  penetration:$('penetration'), avgConc:$('avgConc'), binding:$('binding'), deff:$('deff'),
  interpret:$('interpretText'), stiffWidthVal:$('stiffWidthVal'), obsTimeVal:$('obsTimeVal'),
  sweepTable:$('sweepTable')
};

const ctx = els.canvas.getContext('2d');
const cctx = els.chart.getContext('2d');

function clamp(x,a,b){return Math.max(a,Math.min(b,x));}
function lerp(a,b,t){return a+(b-a)*t;}

function stiffnessFromAgar(p){
  // 모델용 calibration curve. 실제 탄성률 측정값이 있으면 직접 입력하세요.
  // 0% -> 0.5 kPa, 2% -> 약 12.5 kPa.
  return 0.5 + 6.0 * p;
}

function coatingFactor(t){
  // 사용자의 0/10/20/40분 결과에서 10분이 가장 좋다는 관찰을 반영한
  // "현상론적" bell-shaped profile. 측정 전하값이 있으면 beadCharge에 직접 입력.
  return Math.exp(-Math.pow((t-10)/13.0,2));
}

function effectiveAttraction(zCell,zBead){
  // 반대 전하일 때만 양의 인력. geometric mean은 mV 단위의 유효 스케일.
  if(zCell*zBead >= 0) return 0;
  return Math.sqrt(Math.abs(zCell*zBead));
}

function hill(A,kd,n){
  if(A<=0) return 0;
  return Math.pow(A,n)/(Math.pow(kd,n)+Math.pow(A,n));
}

function stiffnessFactor(E,Eopt,width){
  if(E<=0 || Eopt<=0) return 0;
  const x=E/Eopt;
  // E/E* · exp(1-E/E*) : E=E*에서 최대 1
  return clamp(x*Math.exp(1-x),0,1) ** (1/Math.max(width,0.1));
}

function getParams(agarOverride=null){
  const t=+els.coatTime.value;
  const agar=agarOverride===null ? +els.agarPct.value : agarOverride;
  const rawE=+els.stiffness.value;
  const E=els.agarProxy.checked ? stiffnessFromAgar(agar) : rawE;
  const coat=coatingFactor(t);

  // 사용자가 입력한 비드 전위를 "코팅 정도"에 따라 조정.
  // 0분에는 표면 개질 효과가 작고 10분 부근에서 최대라는 실험 관찰을 반영.
  const zBead=+els.beadCharge.value * (0.35 + 0.65*coat);
  const zCell=+els.cellCharge.value;
  const A=effectiveAttraction(zCell,zBead);
  const theta=hill(A,+els.kd.value,+els.hillN.value);
  const sFactor=stiffnessFactor(E,+els.optStiffness.value,+els.stiffWidth.value);
  const D0=+els.d0.value;
  const Ed=+els.ed.value;
  const Deff=D0/(1+Math.pow(E/Ed,2));

  return {t,agar,E,coat,zBead,zCell,A,theta,sFactor,D0,Ed,Deff,obs:+els.obsTime.value};
}

function simulateDiffusion(p){
  // 2D 원형 비드 단면. Fick 2차 방정식을 explicit finite difference로 계산.
  const N=61, radius=2.0; // mm
  const dx=radius/(N/2);
  const D=p.Deff;
  let dt=0.18*dx*dx/(4*D);
  dt=clamp(dt,0.005,0.2);
  const steps=Math.ceil(p.obs/dt);
  dt=p.obs/steps;

  let a=new Float32Array(N*N), b=new Float32Array(N*N);
  const inside=new Uint8Array(N*N);
  const cx=(N-1)/2, cy=(N-1)/2, r2=(N/2-1)*(N/2-1);
  for(let y=0;y<N;y++){
    for(let x=0;x<N;x++){
      const i=y*N+x;
      const dxp=x-cx, dyp=y-cy;
      inside[i]=(dxp*dxp+dyp*dyp<=r2)?1:0;
    }
  }
  const alpha=D*dt/(dx*dx);
  const safeAlpha=Math.min(alpha,0.24);

  for(let step=0;step<steps;step++){
    for(let y=1;y<N-1;y++){
      for(let x=1;x<N-1;x++){
        const i=y*N+x;
        if(!inside[i]) {b[i]=0;continue;}
        const nearBoundary = !inside[i-1]||!inside[i+1]||!inside[i-N]||!inside[i+N];
        if(nearBoundary){b[i]=1;continue;}
        const lap=a[i-1]+a[i+1]+a[i-N]+a[i+N]-4*a[i];
        b[i]=clamp(a[i]+safeAlpha*lap,0,1);
      }
    }
    // fixed concentration at outer boundary
    for(let i=0;i<N;i++){ b[i*N]=inside[i*N]?1:0; b[i*N+N-1]=inside[i*N+N-1]?1:0; }
    for(let x=0;x<N;x++){ b[x]=inside[x]?1:0; b[(N-1)*N+x]=inside[(N-1)*N+x]?1:0; }
    const tmp=a;a=b;b=tmp;
  }

  // Binding is represented as a local retention/readout layer on top of Fick concentration.
  // This intentionally separates transport (Fick) from apparent binding (Hill-type).
  let sum=0, weighted=0, count=0, penetrated=0;
  for(let i=0;i<N*N;i++){
    if(!inside[i]) continue;
    const c=a[i];
    const boundReadout=c*(0.35+0.65*p.theta)*(.65+0.35*p.sFactor);
    sum+=c; weighted+=boundReadout; count++;
    if(c>0.10) penetrated++;
  }
  const avg=sum/count;
  const readout=weighted/count;
  const depthScore=penetrated/count;
  // Overall index: transport + binding + stiffness hypothesis.
  const index=clamp(100*(0.55*depthScore+0.30*avg+0.15*readout),0,100);
  return {N,a,inside,avg,readout,depthScore,index};
}

function colorFor(v){
  // v=0 -> pale, v=1 -> deep blue
  v=clamp(v,0,1);
  const r=Math.round(235*(1-v)+35*v);
  const g=Math.round(247*(1-v)+104*v);
  const b=Math.round(250*(1-v)+180*v);
  return `rgb(${r},${g},${b})`;
}

function drawSimulation(sim,p){
  const can=els.canvas, W=can.width,H=can.height;
  ctx.clearRect(0,0,W,H);
  ctx.fillStyle='#f8fafc';ctx.fillRect(0,0,W,H);

  const N=sim.N, pad=42, size=Math.min(W,H)-2*pad, cell=size/N;
  for(let y=0;y<N;y++){
    for(let x=0;x<N;x++){
      const i=y*N+x;
      if(!sim.inside[i]) continue;
      const c=sim.a[i];
      const visible=clamp(c*(0.35+0.65*p.theta)*(0.55+0.45*p.sFactor),0,1);
      ctx.fillStyle=colorFor(visible);
      ctx.fillRect(pad+x*cell,pad+y*cell,cell+.6,cell+.6);
    }
  }

  // outer boundary
  ctx.strokeStyle='#25313c';ctx.lineWidth=3;
  ctx.beginPath();ctx.arc(W/2,H/2,size/2,0,Math.PI*2);ctx.stroke();

  // center marker
  ctx.fillStyle='rgba(255,255,255,.75)';
  ctx.beginPath();ctx.arc(W/2,H/2,6,0,Math.PI*2);ctx.fill();
  ctx.fillStyle='#35414c';ctx.font='12px Segoe UI';
  ctx.textAlign='center';ctx.fillText('비드/조직 중심',W/2,H/2+26);

  ctx.fillStyle='#55616d';ctx.font='12px Segoe UI';
  ctx.textAlign='left';ctx.fillText(`관찰 시간 ${p.obs} min`,pad,25);
  ctx.textAlign='right';ctx.fillText(`D_eff = ${p.Deff.toFixed(4)} mm²/min`,W-pad,25);
}

function updateMetrics(p,sim){
  els.penetration.textContent=sim.index.toFixed(1);
  els.avgConc.textContent=(sim.avg*100).toFixed(1);
  els.binding.textContent=(p.theta*100).toFixed(1)+'%';
  els.deff.textContent=p.Deff.toFixed(4);

  let chargeText;
  if(p.zCell*p.zBead<0) chargeText=`반대 전하라서 정전기적 인력이 생기며, 유효 인력 스케일 A=${p.A.toFixed(1)} mV입니다.`;
  else chargeText=`같은 부호의 전하라 정전기적 인력 항이 0으로 처리됩니다.`;

  els.interpret.innerHTML =
    `<b>${p.agar}% agar, ${p.t}분 코팅, 조직 강성 ${p.E.toFixed(2)} kPa</b>에서 `+
    `${chargeText} Hill-type 결합률은 <b>${(p.theta*100).toFixed(1)}%</b>로 계산됩니다. `+
    `Fick 확산만 보면 유효 확산계수는 <b>${p.Deff.toFixed(4)} mm²/min</b>이며, `+
    `최종 유효 침투 지수는 <b>${sim.index.toFixed(1)}/100</b>입니다. `+
    `강성이 최적점보다 너무 높아지면 확산 저하 때문에 전체 지수가 다시 낮아질 수 있도록 설계했습니다.`;
}

function drawChart(data){
  const W=els.chart.width,H=els.chart.height;
  cctx.clearRect(0,0,W,H);cctx.fillStyle='#fff';cctx.fillRect(0,0,W,H);
  const left=65,right=30,top=35,bottom=55;
  const pw=W-left-right,ph=H-top-bottom;
  cctx.strokeStyle='#e3e8ed';cctx.lineWidth=1;
  for(let y=0;y<=5;y++){
    const yy=top+ph*y/5;cctx.beginPath();cctx.moveTo(left,yy);cctx.lineTo(W-right,yy);cctx.stroke();
    cctx.fillStyle='#7b858f';cctx.font='11px Segoe UI';cctx.textAlign='right';
    cctx.fillText((100-20*y).toString(),left-8,yy+4);
  }
  cctx.strokeStyle='#59636d';cctx.beginPath();cctx.moveTo(left,top);cctx.lineTo(left,top+ph);cctx.lineTo(W-right,top+ph);cctx.stroke();

  const maxX=data.length-1;
  cctx.strokeStyle='#315caa';cctx.lineWidth=3;cctx.beginPath();
  data.forEach((d,i)=>{
    const x=left+(maxX?pw*i/maxX:pw/2);
    const y=top+ph*(1-d.index/100);
    if(i===0)cctx.moveTo(x,y);else cctx.lineTo(x,y);
  });cctx.stroke();

  data.forEach((d,i)=>{
    const x=left+(maxX?pw*i/maxX:pw/2), y=top+ph*(1-d.index/100);
    cctx.fillStyle='#fff';cctx.strokeStyle='#315caa';cctx.lineWidth=2;
    cctx.beginPath();cctx.arc(x,y,5,0,Math.PI*2);cctx.fill();cctx.stroke();
    cctx.fillStyle='#4e5964';cctx.font='11px Segoe UI';cctx.textAlign='center';
    cctx.fillText(d.label,x,top+ph+22);
  });
  cctx.fillStyle='#4e5964';cctx.font='12px Segoe UI';cctx.textAlign='left';
  cctx.fillText('유효 침투 지수 (0~100)',left,18);
}

function makeTable(data){
  let h='<table><thead><tr><th>한천</th><th>모델 강성</th><th>D_eff</th><th>Hill 결합률</th><th>침투 지수</th></tr></thead><tbody>';
  for(const d of data) h+=`<tr><td>${d.label}</td><td>${d.E.toFixed(2)} kPa</td><td>${d.Deff.toFixed(4)}</td><td>${(d.theta*100).toFixed(1)}%</td><td><b>${d.index.toFixed(1)}</b></td></tr>`;
  h+='</tbody></table>';
  els.sweepTable.innerHTML=h;
}

function run(){
  els.status.textContent='계산 중…';els.status.className='status';
  requestAnimationFrame(()=>{
    const p=getParams();
    const sim=simulateDiffusion(p);
    drawSimulation(sim,p);updateMetrics(p,sim);
    els.status.textContent='완료';els.status.className='status ok';
  });
}

function sweep(){
  els.status.textContent='조건 비교 중…';els.status.className='status';
  requestAnimationFrame(()=>{
    const arr=[0,0.5,1,1.5,2].map(a=>{
      const p=getParams(a),sim=simulateDiffusion(p);
      return {label:a+'%',E:p.E,Deff:p.Deff,theta:p.theta,index:sim.index};
    });
    drawChart(arr);makeTable(arr);
    els.status.textContent='5조건 비교 완료';els.status.className='status ok';
  });
}

els.stiffWidth.addEventListener('input',()=>els.stiffWidthVal.textContent=(+els.stiffWidth.value).toFixed(2));
els.obsTime.addEventListener('input',()=>els.obsTimeVal.textContent=els.obsTime.value);
els.runBtn.addEventListener('click',run);
els.sweepBtn.addEventListener('click',sweep);

els.presetBtn.addEventListener('click',()=>{
  els.coatTime.value='10';
  els.beadCharge.value='-25';
  els.cellCharge.value='-30';
  els.mbConc.value='10';
  els.agarPct.value='2';
  els.stiffness.value='8';
  els.optStiffness.value='20';
  els.agarProxy.checked=true;
  els.d0.value='0.020';
  els.ed.value='20';
  els.kd.value='25';
  els.hillN.value='2';
  els.stiffWidth.value='1';
  els.obsTime.value='20';
  els.stiffWidthVal.textContent='1.00';els.obsTimeVal.textContent='20';
  run();
});

els.presetBtn.click();
