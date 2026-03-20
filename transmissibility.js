// ─────────────────────────────────────────────────────────────────────────────
// transmissibility.js  —  iMSS Group
// Interactive tip transmissibility of a graded piezoelectric metastructure
// Ref: AlShaqaq & Erturk, Smart Materials and Structures 30, 015029 (2021)
// ─────────────────────────────────────────────────────────────────────────────

(function () {

  const canvas = document.getElementById('tr-canvas');
  if (!canvas) return;
  const ctx = canvas.getContext('2d');

  // ── Physical parameters (SI) ──────────────────────────────────────────────
  const hs = 1e-4, bs = 1e-2, rho_s = 2700, cs = 69e9;
  const hp = 3e-4,             rho_p = 7750,
        e31 = -12.3, c11 = 61e9, e33 = 13.3e-9;
  const L  = 0.1;
  const NM = 20;
  const Q  = 25;
  const SZ = NM + Q;

  const m      = 2*rho_p*bs*hp + rho_s*bs*hs;
  const EI_b   = (2*bs/3)*(cs*hs**3/8 + c11*((hp+hs/2)**3 - hs**3/8));
  const Cp_hat = e33*bs/(2*hp);
  const TH     = e31*bs/(2*hp)*((hp+hs/2)**2 - hs**2/4);
  const sqrtmL = Math.sqrt(m*L);

  const xql = Array.from({length:Q}, (_,k) => k/Q);
  const xqr = Array.from({length:Q}, (_,k) => (k+1)/Q);
  const Cpq  = xql.map((_,k) => Cp_hat*L*(xqr[k]-xql[k]));
  const Cp   = Cpq[0];

  // ── Eigenvalues ───────────────────────────────────────────────────────────
  const BETA_EXACT = [1.875104068711961,4.694091132974175,
                      7.854757438237613,10.99554073487547,14.13716839104647];
  function beta(n){ return n<5?BETA_EXACT[n]:(2*(n+1)-1)*Math.PI/2; }
  function omgn_fn(n){ const b=beta(n); return b*b*Math.sqrt(EI_b/(m*L**4)); }

  function modeAndDeriv(n,x){
    const B=beta(n);
    const sig=(Math.sin(B)-Math.sinh(B))/(Math.cos(B)+Math.cosh(B));
    const phi =  Math.cos(B*x)-Math.cosh(B*x)+sig*(Math.sin(B*x)-Math.sinh(B*x));
    const dphi=(B/L)*(-Math.sin(B*x)-Math.sinh(B*x)+sig*(Math.cos(B*x)-Math.cosh(B*x)));
    return {phi,dphi};
  }
  function modeAvg(n){
    const B=beta(n);
    return 2*(Math.sin(B)-Math.sinh(B))/(B*(Math.cos(B)+Math.cosh(B)));
  }

  // ── Precomputed quantities ────────────────────────────────────────────────
  const omgn   = Array.from({length:NM},(_,n)=>omgn_fn(n));
  const omg1   = omgn[0];
  const phiTip = Array.from({length:NM},(_,n)=>modeAndDeriv(n,1.0).phi);
  const phiAvg = Array.from({length:NM},(_,n)=>modeAvg(n));

  const dphi_diff = Array.from({length:Q},(_,k)=>
    Array.from({length:NM},(_,n)=>
      modeAndDeriv(n,xqr[k]).dphi-modeAndDeriv(n,xql[k]).dphi));

  const coupling = Array.from({length:Q},(_,k)=>
    Array.from({length:NM},(_,n)=>TH*dphi_diff[k][n]/sqrtmL));

  // ── Grading law ───────────────────────────────────────────────────────────
  function shuntFreq(omgt,delta,p){
    return Array.from({length:Q},(_,k)=>{
      if(delta===0||p===0) return omgt;
      return omgt+delta-2*delta*(k/(Q-1))**p;
    });
  }

  // ── Solver ────────────────────────────────────────────────────────────────
  const A_re=new Float64Array(SZ*SZ), A_im=new Float64Array(SZ*SZ);
  const b_re=new Float64Array(SZ),    b_im=new Float64Array(SZ);

  function solveGE(){
    const n=SZ;
    for(let col=0;col<n;col++){
      let maxV=0,pivRow=col;
      for(let row=col;row<n;row++){
        const r=A_re[row*n+col],i=A_im[row*n+col],v=r*r+i*i;
        if(v>maxV){maxV=v;pivRow=row;}
      }
      if(pivRow!==col){
        for(let j=0;j<n;j++){
          let t=A_re[col*n+j];A_re[col*n+j]=A_re[pivRow*n+j];A_re[pivRow*n+j]=t;
          t=A_im[col*n+j];A_im[col*n+j]=A_im[pivRow*n+j];A_im[pivRow*n+j]=t;
        }
        let t=b_re[col];b_re[col]=b_re[pivRow];b_re[pivRow]=t;
        t=b_im[col];b_im[col]=b_im[pivRow];b_im[pivRow]=t;
      }
      const pr=A_re[col*n+col],pi=A_im[col*n+col],pd=pr*pr+pi*pi;
      for(let row=col+1;row<n;row++){
        const fr=A_re[row*n+col],fi=A_im[row*n+col];
        const facR=(fr*pr+fi*pi)/pd,facI=(fi*pr-fr*pi)/pd;
        for(let j=col;j<n;j++){
          const ar=A_re[col*n+j],ai=A_im[col*n+j];
          A_re[row*n+j]-=facR*ar-facI*ai;
          A_im[row*n+j]-=facR*ai+facI*ar;
        }
        const br=b_re[col],bi=b_im[col];
        b_re[row]-=facR*br-facI*bi;
        b_im[row]-=facR*bi+facI*br;
      }
    }
    for(let i=n-1;i>=0;i--){
      for(let j=i+1;j<n;j++){
        const ar=A_re[i*n+j],ai=A_im[i*n+j];
        b_re[i]-=ar*b_re[j]-ai*b_im[j];
        b_im[i]-=ar*b_im[j]+ai*b_re[j];
      }
      const pr=A_re[i*n+i],pi=A_im[i*n+i],pd=pr*pr+pi*pi;
      const xr=b_re[i],xi=b_im[i];
      b_re[i]=(xr*pr+xi*pi)/pd;
      b_im[i]=(xi*pr-xr*pi)/pd;
    }
  }

  function computeAtFreq(omg,omge,tau,zeta_r,omgt){
    const os2=omg*omg;
    A_re.fill(0); A_im.fill(0); b_re.fill(0); b_im.fill(0);
    for(let n=0;n<NM;n++){
      const wn=omgn[n];
      A_re[n*SZ+n]=wn*wn-os2;
      A_im[n*SZ+n]=2*zeta_r*wn*omg;
      b_re[n]=sqrtmL*os2*phiAvg[n];
    }
    for(let k=0;k<Q;k++){
      const we=omge[k];
      const rD=Cpq[k]*omgt/(tau*Cp);
      A_re[(NM+k)*SZ+(NM+k)]=rD;
      A_im[(NM+k)*SZ+(NM+k)]=omg-we*we/omg;
      for(let n=0;n<NM;n++){
        const c=coupling[k][n];
        A_re[n*SZ+(NM+k)]=-c;
        A_im[(NM+k)*SZ+n]=omg*c/Cp;
      }
    }
    solveGE();
    let dre=1,dim=0;
    for(let n=0;n<NM;n++){
      dre+=phiTip[n]*b_re[n]/sqrtmL;
      dim+=phiTip[n]*b_im[n]/sqrtmL;
    }
    return Math.sqrt(dre*dre+dim*dim);
  }

  function computeSCatFreq(omg,zeta_r){
    const os2=omg*omg;
    let dre=1,dim=0;
    for(let n=0;n<NM;n++){
      const wn=omgn[n],rhs=sqrtmL*os2*phiAvg[n];
      const dr=wn*wn-os2,di=2*zeta_r*wn*omg,dsq=dr*dr+di*di;
      dre+=phiTip[n]*(rhs*dr/dsq)/sqrtmL;
      dim+=phiTip[n]*(-rhs*di/dsq)/sqrtmL;
    }
    return Math.sqrt(dre*dre+dim*dim);
  }

  function computeTR(omge,tau,zeta_r,omg_arr,omgt){
    const out=new Float64Array(omg_arr.length);
    for(let fi=0;fi<omg_arr.length;fi++)
      out[fi]=computeAtFreq(omg_arr[fi],omge,tau,zeta_r,omgt);
    return out;
  }

  function computeSC(zeta_r,omg_arr){
    const out=new Float64Array(omg_arr.length);
    for(let fi=0;fi<omg_arr.length;fi++)
      out[fi]=computeSCatFreq(omg_arr[fi],zeta_r);
    return out;
  }

  // ── Bandgap width: largest CONTINUOUS segment where TR_graded < 0.1 ───────
  // If a peak breaks continuity inside the main region → that segment = 0
  function computeBandgapWidth(omg_norm, tr_grad, threshold){
    let bestWidth = 0;
    let segStart  = null;

    for(let i=0; i<tr_grad.length; i++){
      if(tr_grad[i] < threshold){
        if(segStart === null) segStart = omg_norm[i];
      } else {
        if(segStart !== null){
          const w = omg_norm[i-1] - segStart;
          if(w > bestWidth) bestWidth = w;
          segStart = null;
        }
      }
    }
    // close last segment
    if(segStart !== null){
      const w = omg_norm[omg_norm.length-1] - segStart;
      if(w > bestWidth) bestWidth = w;
    }
    return bestWidth;
  }

  // ── Canvas coordinates ────────────────────────────────────────────────────
  const LOG_MIN=-5, LOG_MAX=3;
  const OMG_MIN=25, OMG_MAX=45;
  const THRESH=0.1;

  function toPixel(omg_norm, tr_val, P, PW, PH){
    const px=P.l+(omg_norm-OMG_MIN)/(OMG_MAX-OMG_MIN)*PW;
    const logV=Math.log10(Math.max(tr_val,1e-10));
    const py=P.t+PH*(LOG_MAX-logV)/(LOG_MAX-LOG_MIN);
    return {px,py};
  }

  // ── Draw shaded area below threshold for graded curve ────────────────────
  function drawShading(omg_norm, tr_grad, P, PW, PH){
    const {py: yThr} = toPixel(0, THRESH, P, PW, PH);
    let inShade=false;

    for(let i=0; i<omg_norm.length; i++){
      const {px,py}=toPixel(omg_norm[i], tr_grad[i], P, PW, PH);
      const clampedPy=Math.max(P.t, Math.min(P.t+PH, py));

      if(tr_grad[i]<THRESH){
        if(!inShade){
          ctx.beginPath();
          ctx.moveTo(px, yThr);
          inShade=true;
        }
        ctx.lineTo(px, clampedPy);
      } else {
        if(inShade){
          ctx.lineTo(px, yThr);
          ctx.closePath();
          ctx.fillStyle='rgba(201,168,76,0.18)';
          ctx.fill();
          inShade=false;
        }
      }
    }
    if(inShade){
      const {px}=toPixel(omg_norm[omg_norm.length-1], 0, P, PW, PH);
      ctx.lineTo(px, yThr);
      ctx.closePath();
      ctx.fillStyle='rgba(201,168,76,0.18)';
      ctx.fill();
    }
  }

  function drawCurve(omg_norm, tr_vals, color, lineWidth, dashed, P, PW, PH){
    ctx.beginPath();
    ctx.strokeStyle=color;
    ctx.lineWidth=lineWidth;
    if(dashed) ctx.setLineDash([5,5]); else ctx.setLineDash([]);
    let first=true;
    for(let i=0;i<omg_norm.length;i++){
      const {px,py}=toPixel(omg_norm[i],tr_vals[i],P,PW,PH);
      if(py<P.t-2||py>P.t+PH+2){first=true;continue;}
      first?(ctx.moveTo(px,py),first=false):ctx.lineTo(px,py);
    }
    ctx.stroke();
    ctx.setLineDash([]);
  }

  function drawPlot(omg_norm, tr_sc, tr_uni, tr_grad, params){
    const W  =canvas.offsetWidth||600;
    const dpr=window.devicePixelRatio||1;
    canvas.width =W  *dpr;
    canvas.height=460*dpr;
    ctx.scale(dpr,dpr);
    const DW=W, DH=460;
    const P={t:28,r:28,b:56,l:72};
    const PW=DW-P.l-P.r, PH=DH-P.t-P.b;

    ctx.clearRect(0,0,DW,DH);

    // Grid
    ctx.strokeStyle='rgba(255,255,255,0.05)'; ctx.lineWidth=1;
    for(let g=LOG_MIN;g<=LOG_MAX;g++){
      const {py}=toPixel(0,Math.pow(10,g),P,PW,PH);
      ctx.beginPath();ctx.moveTo(P.l,py);ctx.lineTo(P.l+PW,py);ctx.stroke();
    }
    for(let g=OMG_MIN;g<=OMG_MAX;g+=5){
      const px=P.l+(g-OMG_MIN)/(OMG_MAX-OMG_MIN)*PW;
      ctx.beginPath();ctx.moveTo(px,P.t);ctx.lineTo(px,P.t+PH);ctx.stroke();
    }

    // Minor decade lines
    for(let g=LOG_MIN;g<LOG_MAX;g++){
      for(let m=2;m<=9;m++){
        const {py}=toPixel(0,m*Math.pow(10,g),P,PW,PH);
        ctx.strokeStyle='rgba(255,255,255,0.03)'; ctx.lineWidth=0.5;
        ctx.beginPath();ctx.moveTo(P.l,py);ctx.lineTo(P.l+PW,py);ctx.stroke();
      }
    }

    // Threshold line
    ctx.setLineDash([4,4]);
    ctx.strokeStyle='rgba(201,168,76,0.4)'; ctx.lineWidth=1;
    const {py:yThr}=toPixel(0,THRESH,P,PW,PH);
    ctx.beginPath();ctx.moveTo(P.l,yThr);ctx.lineTo(P.l+PW,yThr);ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle='rgba(201,168,76,0.55)';
    ctx.font='10px Poppins,sans-serif'; ctx.textAlign='left';
    ctx.fillText('|TR| = 0.1',P.l+5,yThr-4);

    // Shaded bandgap area
    drawShading(omg_norm, tr_grad, P, PW, PH);

    // Curves (back to front: SC → uniform → graded)
    drawCurve(omg_norm,tr_sc,  'rgba(255,255,255,0.3)',  1.5,true, P,PW,PH);
    drawCurve(omg_norm,tr_uni, 'rgba(255,255,255,0.8)',  2.0,false,P,PW,PH);
    drawCurve(omg_norm,tr_grad,'#C9A84C',                2.5,false,P,PW,PH);

    // Axes
    ctx.strokeStyle='rgba(255,255,255,0.2)'; ctx.lineWidth=1;
    ctx.beginPath();
    ctx.moveTo(P.l,P.t);ctx.lineTo(P.l,P.t+PH);ctx.lineTo(P.l+PW,P.t+PH);
    ctx.stroke();

    // x-axis label & ticks
    ctx.fillStyle='rgba(255,255,255,0.5)';
    ctx.font='13px Poppins,sans-serif'; ctx.textAlign='center';
    ctx.fillText('Frequency  \u03C9/\u03C9\u2081',P.l+PW/2,DH-8);
    ctx.font='11px Poppins,sans-serif';
    for(let g=OMG_MIN;g<=OMG_MAX;g+=5){
      const px=P.l+(g-OMG_MIN)/(OMG_MAX-OMG_MIN)*PW;
      ctx.fillText(g,px,P.t+PH+17);
    }

    // y-axis label & ticks
    ctx.save(); ctx.translate(14,P.t+PH/2); ctx.rotate(-Math.PI/2);
    ctx.font='13px Poppins,sans-serif';
    ctx.fillText('|TR| (tip transmissibility)',0,0); ctx.restore();
    ctx.font='11px Poppins,sans-serif'; ctx.textAlign='right';
    const yLabels=['10\u207B\u2075','10\u207B\u2074','10\u207B\u00B3',
                   '10\u207B\u00B2','10\u207B\u00B9','10\u2070',
                   '10\u00B9','10\u00B2','10\u00B3'];
    for(let g=LOG_MIN;g<=LOG_MAX;g++){
      const {py}=toPixel(0,Math.pow(10,g),P,PW,PH);
      ctx.fillStyle='rgba(255,255,255,0.4)';
      ctx.fillText(yLabels[g-LOG_MIN],P.l-8,py+3);
    }

    // ── Bandgap width ───────────────────────────────────────────────────────
    const bgWidth=computeBandgapWidth(omg_norm, tr_grad, THRESH);
    const bgText = bgWidth>0
      ? 'Bandgap width: \u0394\u03C9/\u03C9\u2081 = '+bgWidth.toFixed(2)
      : 'Bandgap: discontinuous (width = 0)';

    // badge background
    const badgeX=P.l+8, badgeY=P.t+8;
    ctx.font='11px Poppins,sans-serif';
    const tw=ctx.measureText(bgText).width;
    ctx.fillStyle='rgba(11,37,69,0.75)';
    ctx.fillRect(badgeX-6, badgeY-13, tw+18, 20);
    ctx.strokeStyle='rgba(201,168,76,0.35)'; ctx.lineWidth=0.75;
    ctx.strokeRect(badgeX-6, badgeY-13, tw+18, 20);
    ctx.fillStyle= bgWidth>0 ? '#C9A84C' : 'rgba(255,255,255,0.4)';
    ctx.textAlign='left';
    ctx.fillText(bgText, badgeX+3, badgeY+1);

    // ── Legend (inside plot, top-right, all items visible) ─────────────────
    const legItems=[
      {label:'Short circuit',                              color:'rgba(255,255,255,0.4)', dash:true},
      {label:'Uniform  (p = 0)',                           color:'rgba(255,255,255,0.9)', dash:false},
      {label:`Graded   (p = ${params.p%1===0?params.p.toFixed(0):params.p.toFixed(2)})`, color:'#C9A84C', dash:false},
    ];
    const legLineW=22, legPad=10, legRowH=20;
    const legTextW=Math.max(...legItems.map(it=>{
      ctx.font='11px Poppins,sans-serif';
      return ctx.measureText(it.label).width;
    }));
    const legW=legLineW+8+legTextW+legPad*2;
    const legH=legItems.length*legRowH+legPad*1.5;
    const legX=P.l+PW-legW-8;
    const legY=P.t+8;

    ctx.fillStyle='rgba(11,37,69,0.82)';
    ctx.fillRect(legX,legY,legW,legH);
    ctx.strokeStyle='rgba(255,255,255,0.1)'; ctx.lineWidth=0.5;
    ctx.strokeRect(legX,legY,legW,legH);

    legItems.forEach((item,i)=>{
      const lx=legX+legPad;
      const ly=legY+legPad+i*legRowH+7;
      ctx.strokeStyle=item.color; ctx.lineWidth=item.dash?1.5:2;
      if(item.dash) ctx.setLineDash([4,4]); else ctx.setLineDash([]);
      ctx.beginPath(); ctx.moveTo(lx,ly); ctx.lineTo(lx+legLineW,ly); ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle='rgba(255,255,255,0.6)';
      ctx.font='11px Poppins,sans-serif'; ctx.textAlign='left';
      ctx.fillText(item.label,lx+legLineW+8,ly+3.5);
    });
  }

  // ── Debounced run ─────────────────────────────────────────────────────────
  let debTimer=null;
  let params={p:1, delta:3, zeta:0.001, tau:500, omgt_r:35};

  function run(){
    const omgt=params.omgt_r*omg1;
    const NF=500;   // high-resolution frequency sweep
    const omg_arr=Array.from({length:NF},(_,i)=>(OMG_MIN+i*(OMG_MAX-OMG_MIN)/(NF-1))*omg1);
    const omg_norm=omg_arr.map(o=>o/omg1);

    const tr_sc  =computeSC(params.zeta,omg_arr);
    const omge_u =shuntFreq(omgt,params.delta*omg1,0);
    const tr_uni =computeTR(omge_u,params.tau,params.zeta,omg_arr,omgt);
    const omge_g =shuntFreq(omgt,params.delta*omg1,params.p);
    const tr_grad=computeTR(omge_g,params.tau,params.zeta,omg_arr,omgt);

    drawPlot(omg_norm,tr_sc,tr_uni,tr_grad,params);
  }

  function schedule(){
    clearTimeout(debTimer);
    debTimer=setTimeout(run,80);
  }

  function wire(sliderId,valId,decimals,key){
    const el=document.getElementById(sliderId);
    if(!el) return;
    el.addEventListener('input',e=>{
      params[key]=+e.target.value;
      const disp=document.getElementById(valId);
      if(disp) disp.textContent=(+e.target.value).toFixed(decimals);
      schedule();
    });
  }

  wire('tr-p',    'tr-p-val',    2,'p');
  wire('tr-delta','tr-delta-val',1,'delta');
  wire('tr-zeta', 'tr-zeta-val', 3,'zeta');
  wire('tr-tau',  'tr-tau-val',  0,'tau');
  wire('tr-omgt', 'tr-omgt-val', 0,'omgt_r');

  run();
  window.addEventListener('resize',schedule);
  window.trRedraw=run;

})();
