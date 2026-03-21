// ─────────────────────────────────────────────────────────────────────────────
// transmissibility.js  —  iMSS Group
// Two-panel interactive figure:
//   Top:    Shunt frequency grading profile ωt,k/ω₁ vs unit cell index k
//   Bottom: Tip transmissibility |TR| vs ω/ω₁
// Ref: AlShaqaq & Erturk, Smart Materials and Structures 30, 015029 (2021)
// ─────────────────────────────────────────────────────────────────────────────

(function () {

  const canvas = document.getElementById('tr-canvas');
  if (!canvas) return;
  const ctx = canvas.getContext('2d');

  // ── Physical parameters (SI) ──────────────────────────────────────────────
  const hs=1e-4, bs=1e-2, rho_s=2700, cs=69e9;
  const hp=3e-4, rho_p=7750, e31=-12.3, c11=61e9, e33=13.3e-9;
  const L=0.1, NM=20, Q=25, SZ=NM+Q;

  const m      = 2*rho_p*bs*hp + rho_s*bs*hs;
  const EI_b   = (2*bs/3)*(cs*hs**3/8 + c11*((hp+hs/2)**3 - hs**3/8));
  const Cp_hat = e33*bs/(2*hp);
  const TH     = e31*bs/(2*hp)*((hp+hs/2)**2 - hs**2/4);
  const sqrtmL = Math.sqrt(m*L);

  const xql = Array.from({length:Q},(_,k)=>k/Q);
  const xqr = Array.from({length:Q},(_,k)=>(k+1)/Q);
  const Cpq = xql.map((_,k)=>Cp_hat*L*(xqr[k]-xql[k]));
  const Cp  = Cpq[0];

  // ── Eigenvalues ───────────────────────────────────────────────────────────
  const BETA=[1.875104068711961,4.694091132974175,
              7.854757438237613,10.99554073487547,14.13716839104647];
  function beta(n){ return n<5?BETA[n]:(2*(n+1)-1)*Math.PI/2; }
  function omgn_fn(n){ const b=beta(n); return b*b*Math.sqrt(EI_b/(m*L**4)); }

  function modeAndDeriv(n,x){
    const B=beta(n), sig=(Math.sin(B)-Math.sinh(B))/(Math.cos(B)+Math.cosh(B));
    return {
      phi:  Math.cos(B*x)-Math.cosh(B*x)+sig*(Math.sin(B*x)-Math.sinh(B*x)),
      dphi: (B/L)*(-Math.sin(B*x)-Math.sinh(B*x)+sig*(Math.cos(B*x)-Math.cosh(B*x)))
    };
  }
  function modeAvg(n){
    const B=beta(n);
    return 2*(Math.sin(B)-Math.sinh(B))/(B*(Math.cos(B)+Math.cosh(B)));
  }

  // ── Precomputed quantities ────────────────────────────────────────────────
  const omgn   = Array.from({length:NM},(_,n)=>omgn_fn(n));
  const omg1   = omgn[0];
  const phiTip = Array.from({length:NM},(_,n)=>modeAndDeriv(n,1).phi);
  const phiAvg = Array.from({length:NM},(_,n)=>modeAvg(n));
  const coupling = Array.from({length:Q},(_,k)=>
    Array.from({length:NM},(_,n)=>
      TH*(modeAndDeriv(n,xqr[k]).dphi-modeAndDeriv(n,xql[k]).dphi)/sqrtmL));

  // ── Grading law (returns NORMALIZED frequencies ωt,k/ω₁) ─────────────────
  function shuntFreqNorm(omgt_r, delta, p){
    return Array.from({length:Q},(_,k)=>{
      if(delta===0 || p===0) return omgt_r;
      return omgt_r + delta - 2*delta*(k/(Q-1))**p;
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

  function computeAtFreq(omg,omge_norm,tau,zeta_r,omgt_r){
    const omgt=omgt_r*omg1;
    const os2=omg*omg;
    A_re.fill(0);A_im.fill(0);b_re.fill(0);b_im.fill(0);
    for(let n=0;n<NM;n++){
      const wn=omgn[n];
      A_re[n*SZ+n]=wn*wn-os2; A_im[n*SZ+n]=2*zeta_r*wn*omg;
      b_re[n]=sqrtmL*os2*phiAvg[n];
    }
    for(let k=0;k<Q;k++){
      const we=omge_norm[k]*omg1;
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
    for(let n=0;n<NM;n++){dre+=phiTip[n]*b_re[n]/sqrtmL;dim+=phiTip[n]*b_im[n]/sqrtmL;}
    return Math.sqrt(dre*dre+dim*dim);
  }

  function computeSCatFreq(omg,zeta_r){
    const os2=omg*omg;let dre=1,dim=0;
    for(let n=0;n<NM;n++){
      const wn=omgn[n],rhs=sqrtmL*os2*phiAvg[n];
      const dr=wn*wn-os2,di=2*zeta_r*wn*omg,dsq=dr*dr+di*di;
      dre+=phiTip[n]*(rhs*dr/dsq)/sqrtmL;
      dim+=phiTip[n]*(-rhs*di/dsq)/sqrtmL;
    }
    return Math.sqrt(dre*dre+dim*dim);
  }

  function computeTR(omge_norm,tau,zeta_r,omg_arr,omgt_r){
    const out=new Float64Array(omg_arr.length);
    for(let fi=0;fi<omg_arr.length;fi++)
      out[fi]=computeAtFreq(omg_arr[fi],omge_norm,tau,zeta_r,omgt_r);
    return out;
  }
  function computeSC(zeta_r,omg_arr){
    const out=new Float64Array(omg_arr.length);
    for(let fi=0;fi<omg_arr.length;fi++) out[fi]=computeSCatFreq(omg_arr[fi],zeta_r);
    return out;
  }

  function computeBandgapWidth(omg_norm,tr_grad,threshold){
    let best=0,segStart=null;
    for(let i=0;i<tr_grad.length;i++){
      if(tr_grad[i]<threshold){
        if(segStart===null) segStart=omg_norm[i];
      } else {
        if(segStart!==null){const w=omg_norm[i-1]-segStart;if(w>best)best=w;segStart=null;}
      }
    }
    if(segStart!==null){const w=omg_norm[omg_norm.length-1]-segStart;if(w>best)best=w;}
    return best;
  }

  // ─────────────────────────────────────────────────────────────────────────
  // DRAWING
  // ─────────────────────────────────────────────────────────────────────────

  // Shared margins so both panels are perfectly aligned
  const ML=72, MR=32;

  // ── Panel A: Grading profile ──────────────────────────────────────────────
  function drawProfilePanel(omge_u_norm, omge_g_norm, omgt_r, delta, top, PH, DW){
    const PW = DW-ML-MR;
    const L0 = ML, T0 = top;

    // y range: normalized frequencies, with padding
    const pad  = Math.max(delta*0.2, 0.5);
    const yMin = omgt_r - delta - pad;
    const yMax = omgt_r + delta + pad;
    const yRng = yMax - yMin;

    function px(k)  { return L0 + (k/(Q-1))*PW; }
    function py(val){ return T0 + PH*(1-(val-yMin)/yRng); }

    // Grid lines
    ctx.lineWidth=1;
    for(let i=0;i<=4;i++){
      const y=T0+i*PH/4;
      ctx.strokeStyle='rgba(255,255,255,0.05)';
      ctx.beginPath();ctx.moveTo(L0,y);ctx.lineTo(L0+PW,y);ctx.stroke();
    }
    for(let k=0;k<Q;k+=6){
      ctx.strokeStyle='rgba(255,255,255,0.04)';
      ctx.beginPath();ctx.moveTo(px(k),T0);ctx.lineTo(px(k),T0+PH);ctx.stroke();
    }

    // Grading range band [ωt-δ, ωt+δ]
    if(delta>0){
      const yHi=py(omgt_r+delta), yLo=py(omgt_r-delta);
      ctx.fillStyle='rgba(255,255,255,0.04)';
      ctx.fillRect(L0, yHi, PW, yLo-yHi);
      // dashed boundary lines
      ctx.setLineDash([4,4]);
      ctx.strokeStyle='rgba(255,255,255,0.15)'; ctx.lineWidth=0.75;
      [yHi,yLo].forEach(y=>{
        ctx.beginPath();ctx.moveTo(L0,y);ctx.lineTo(L0+PW,y);ctx.stroke();
      });
      ctx.setLineDash([]);
    }

    // Uniform reference line (white)
    ctx.strokeStyle='rgba(255,255,255,0.75)'; ctx.lineWidth=2;
    ctx.setLineDash([]);
    const yU=py(omge_u_norm[0]);
    ctx.beginPath();ctx.moveTo(L0,yU);ctx.lineTo(L0+PW,yU);ctx.stroke();

    // Graded profile — filled area between curve and uniform line
    ctx.beginPath();
    ctx.moveTo(px(0), yU);
    for(let k=0;k<Q;k++) ctx.lineTo(px(k), py(omge_g_norm[k]));
    ctx.lineTo(px(Q-1), yU);
    ctx.closePath();
    ctx.fillStyle='rgba(201,168,76,0.15)'; ctx.fill();

    // Graded profile line
    ctx.strokeStyle='#C9A84C'; ctx.lineWidth=2.2; ctx.setLineDash([]);
    ctx.beginPath();
    for(let k=0;k<Q;k++){
      k===0?ctx.moveTo(px(k),py(omge_g_norm[k])):ctx.lineTo(px(k),py(omge_g_norm[k]));
    }
    ctx.stroke();

    // Dots on graded profile
    for(let k=0;k<Q;k++){
      ctx.beginPath();ctx.arc(px(k),py(omge_g_norm[k]),2.8,0,Math.PI*2);
      ctx.fillStyle='#C9A84C';ctx.fill();
    }

    // Axes
    ctx.strokeStyle='rgba(255,255,255,0.2)';ctx.lineWidth=1;
    ctx.beginPath();
    ctx.moveTo(L0,T0);ctx.lineTo(L0,T0+PH);ctx.lineTo(L0+PW,T0+PH);
    ctx.stroke();

    // x-axis ticks & label
    ctx.fillStyle='rgba(255,255,255,0.45)';
    ctx.font='12px Poppins,sans-serif';ctx.textAlign='center';
    ctx.fillText('Unit cell  k', L0+PW/2, T0+PH+20);
    ctx.font='11px Poppins,sans-serif';
    [1,7,13,19,25].forEach(k=>{
      ctx.fillStyle='rgba(255,255,255,0.4)';
      ctx.fillText(k, px(k-1), T0+PH+16);
    });

    // y-axis ticks & label
    ctx.save();ctx.translate(14,T0+PH/2);ctx.rotate(-Math.PI/2);
    ctx.font='12px Poppins,sans-serif';ctx.fillStyle='rgba(255,255,255,0.45)';
    ctx.fillText('\u03C9\u209C,\u2096 / \u03C9\u2081',0,0);ctx.restore();

    ctx.font='11px Poppins,sans-serif';ctx.textAlign='right';
    // show ωt-δ, ωt, ωt+δ on y-axis
    const yTicks = delta>0
      ? [[omgt_r-delta,'\u03C9\u209C\u2212\u03B4'],[omgt_r,'\u03C9\u209C'],[omgt_r+delta,'\u03C9\u209C+\u03B4']]
      : [[omgt_r,'\u03C9\u209C']];
    yTicks.forEach(([v,lbl])=>{
      ctx.fillStyle='rgba(255,255,255,0.45)';
      ctx.fillText(lbl, L0-8, py(v)+3);
      ctx.strokeStyle='rgba(255,255,255,0.15)';ctx.lineWidth=0.5;
      ctx.beginPath();ctx.moveTo(L0-3,py(v));ctx.lineTo(L0,py(v));ctx.stroke();
    });

    // Panel label
    ctx.fillStyle='rgba(255,255,255,0.25)';ctx.font='10px Poppins,sans-serif';
    ctx.textAlign='left';
    ctx.fillText('(a)  Shunt frequency grading profile', L0+6, T0+13);

    // Legend
    const lgItems=[
      {label:'Uniform',color:'rgba(255,255,255,0.85)',dash:false},
      {label:'Graded', color:'#C9A84C',               dash:false},
    ];
    const lgX=L0+PW-4, lgY=T0+6;
    ctx.font='11px Poppins,sans-serif';
    const lgTW=Math.max(...lgItems.map(it=>ctx.measureText(it.label).width));
    const lgW=26+8+lgTW+14, lgH=lgItems.length*20+10;
    ctx.fillStyle='rgba(11,37,69,0.85)';ctx.fillRect(lgX-lgW,lgY,lgW,lgH);
    ctx.strokeStyle='rgba(255,255,255,0.08)';ctx.lineWidth=0.5;
    ctx.strokeRect(lgX-lgW,lgY,lgW,lgH);
    lgItems.forEach((item,i)=>{
      const lly=lgY+9+i*20;
      ctx.strokeStyle=item.color;ctx.lineWidth=2;
      ctx.beginPath();ctx.moveTo(lgX-lgW+8,lly);ctx.lineTo(lgX-lgW+30,lly);ctx.stroke();
      ctx.fillStyle='rgba(255,255,255,0.65)';ctx.textAlign='left';
      ctx.fillText(item.label,lgX-lgW+36,lly+3.5);
    });
  }

  // ── Panel B: Transmissibility ─────────────────────────────────────────────
  const LOG_MIN=-5, LOG_MAX=3, OMG_MIN=25, OMG_MAX=45, THRESH=0.1;

  function toPxTR(omg_norm,tr_val,top,PH,DW){
    const px=ML+(omg_norm-OMG_MIN)/(OMG_MAX-OMG_MIN)*(DW-ML-MR);
    const logV=Math.log10(Math.max(tr_val,1e-10));
    const py=top+PH*(LOG_MAX-logV)/(LOG_MAX-LOG_MIN);
    return {px,py};
  }

  function drawShading(omg_norm,tr_grad,top,PH,DW){
    const PW=DW-ML-MR;
    const {py:yThr}=toPxTR(0,THRESH,top,PH,DW);
    let inShade=false;
    for(let i=0;i<omg_norm.length;i++){
      const {px,py}=toPxTR(omg_norm[i],tr_grad[i],top,PH,DW);
      const cpy=Math.max(top,Math.min(top+PH,py));
      if(tr_grad[i]<THRESH){
        if(!inShade){ctx.beginPath();ctx.moveTo(px,yThr);inShade=true;}
        ctx.lineTo(px,cpy);
      } else {
        if(inShade){
          ctx.lineTo(px,yThr);ctx.closePath();
          ctx.fillStyle='rgba(201,168,76,0.18)';ctx.fill();inShade=false;
        }
      }
    }
    if(inShade){
      const {px}=toPxTR(omg_norm[omg_norm.length-1],0,top,PH,DW);
      ctx.lineTo(px,yThr);ctx.closePath();
      ctx.fillStyle='rgba(201,168,76,0.18)';ctx.fill();
    }
  }

  function drawCurveTR(omg_norm,tr_vals,color,lw,dashed,top,PH,DW){
    ctx.beginPath();ctx.strokeStyle=color;ctx.lineWidth=lw;
    if(dashed)ctx.setLineDash([5,5]);else ctx.setLineDash([]);
    let first=true;
    for(let i=0;i<omg_norm.length;i++){
      const {px,py}=toPxTR(omg_norm[i],tr_vals[i],top,PH,DW);
      if(py<top-2||py>top+PH+2){first=true;continue;}
      first?(ctx.moveTo(px,py),first=false):ctx.lineTo(px,py);
    }
    ctx.stroke();ctx.setLineDash([]);
  }

  function drawTRPanel(omg_norm,tr_sc,tr_uni,tr_grad,params,top,PH,DW,DH){
    const PW=DW-ML-MR;

    // Grid
    ctx.lineWidth=1;
    for(let g=LOG_MIN;g<=LOG_MAX;g++){
      const {py}=toPxTR(0,Math.pow(10,g),top,PH,DW);
      ctx.strokeStyle='rgba(255,255,255,0.05)';
      ctx.beginPath();ctx.moveTo(ML,py);ctx.lineTo(ML+PW,py);ctx.stroke();
      for(let mm=2;mm<=9;mm++){
        const {py:pym}=toPxTR(0,mm*Math.pow(10,g),top,PH,DW);
        ctx.strokeStyle='rgba(255,255,255,0.02)';ctx.lineWidth=0.5;
        ctx.beginPath();ctx.moveTo(ML,pym);ctx.lineTo(ML+PW,pym);ctx.stroke();
        ctx.lineWidth=1;
      }
    }
    for(let g=OMG_MIN;g<=OMG_MAX;g+=5){
      const px=ML+(g-OMG_MIN)/(OMG_MAX-OMG_MIN)*PW;
      ctx.strokeStyle='rgba(255,255,255,0.05)';ctx.lineWidth=1;
      ctx.beginPath();ctx.moveTo(px,top);ctx.lineTo(px,top+PH);ctx.stroke();
    }

    // Threshold
    ctx.setLineDash([4,4]);
    ctx.strokeStyle='rgba(201,168,76,0.4)';ctx.lineWidth=1;
    const {py:yThr}=toPxTR(0,THRESH,top,PH,DW);
    ctx.beginPath();ctx.moveTo(ML,yThr);ctx.lineTo(ML+PW,yThr);ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle='rgba(201,168,76,0.6)';
    ctx.font='10px Poppins,sans-serif';ctx.textAlign='left';
    ctx.fillText('|TR| = 0.1',ML+5,yThr-4);

    // Shading + curves
    drawShading(omg_norm,tr_grad,top,PH,DW);
    drawCurveTR(omg_norm,tr_sc,  'rgba(255,255,255,0.3)',  1.5,true, top,PH,DW);
    drawCurveTR(omg_norm,tr_uni, 'rgba(255,255,255,0.85)', 2.0,false,top,PH,DW);
    drawCurveTR(omg_norm,tr_grad,'#C9A84C',                2.5,false,top,PH,DW);

    // Axes
    ctx.strokeStyle='rgba(255,255,255,0.2)';ctx.lineWidth=1;
    ctx.beginPath();
    ctx.moveTo(ML,top);ctx.lineTo(ML,top+PH);ctx.lineTo(ML+PW,top+PH);
    ctx.stroke();

    // x-axis
    ctx.fillStyle='rgba(255,255,255,0.5)';
    ctx.font='13px Poppins,sans-serif';ctx.textAlign='center';
    ctx.fillText('Frequency  \u03C9/\u03C9\u2081',ML+PW/2,DH-8);
    ctx.font='11px Poppins,sans-serif';
    for(let g=OMG_MIN;g<=OMG_MAX;g+=5){
      const px=ML+(g-OMG_MIN)/(OMG_MAX-OMG_MIN)*PW;
      ctx.fillStyle='rgba(255,255,255,0.4)';ctx.textAlign='center';
      ctx.fillText(g,px,top+PH+17);
    }

    // y-axis
    ctx.save();ctx.translate(14,top+PH/2);ctx.rotate(-Math.PI/2);
    ctx.font='13px Poppins,sans-serif';ctx.fillStyle='rgba(255,255,255,0.45)';
    ctx.fillText('|TR| (tip transmissibility)',0,0);ctx.restore();
    const yL=['10\u207B\u2075','10\u207B\u2074','10\u207B\u00B3','10\u207B\u00B2',
              '10\u207B\u00B9','10\u2070','10\u00B9','10\u00B2','10\u00B3'];
    ctx.font='11px Poppins,sans-serif';ctx.textAlign='right';
    for(let g=LOG_MIN;g<=LOG_MAX;g++){
      const {py}=toPxTR(0,Math.pow(10,g),top,PH,DW);
      ctx.fillStyle='rgba(255,255,255,0.4)';
      ctx.fillText(yL[g-LOG_MIN],ML-8,py+3);
    }

    // Bandgap width badge
    const bgW=computeBandgapWidth(omg_norm,tr_grad,THRESH);
    const bgTxt=bgW>0
      ?'Bandgap  \u0394\u03C9/\u03C9\u2081 = '+bgW.toFixed(2)
      :'Bandgap: discontinuous';
    ctx.font='11px Poppins,sans-serif';ctx.textAlign='left';
    const tw=ctx.measureText(bgTxt).width;
    ctx.fillStyle='rgba(11,37,69,0.85)';ctx.fillRect(ML+6,top+8,tw+18,22);
    ctx.strokeStyle='rgba(201,168,76,0.35)';ctx.lineWidth=0.75;
    ctx.strokeRect(ML+6,top+8,tw+18,22);
    ctx.fillStyle=bgW>0?'#C9A84C':'rgba(255,255,255,0.4)';
    ctx.fillText(bgTxt,ML+15,top+23);

    // Panel label
    ctx.fillStyle='rgba(255,255,255,0.25)';ctx.font='10px Poppins,sans-serif';
    ctx.textAlign='left';
    ctx.fillText('(b)  Tip transmissibility',ML+6,top+PH-8);

    // Legend
    const lgItems=[
      {label:'Short circuit',                                              color:'rgba(255,255,255,0.4)',dash:true},
      {label:'Uniform  (p = 0)',                                           color:'rgba(255,255,255,0.9)',dash:false},
      {label:`Graded  (p = ${params.p%1===0?params.p.toFixed(0):params.p.toFixed(2)})`,color:'#C9A84C',dash:false},
    ];
    ctx.font='11px Poppins,sans-serif';
    const lgTW=Math.max(...lgItems.map(it=>ctx.measureText(it.label).width));
    const lgW=26+8+lgTW+14, lgH=lgItems.length*20+12;
    const lgX=ML+PW-4, lgY=top+6;
    ctx.fillStyle='rgba(11,37,69,0.85)';ctx.fillRect(lgX-lgW,lgY,lgW,lgH);
    ctx.strokeStyle='rgba(255,255,255,0.08)';ctx.lineWidth=0.5;ctx.strokeRect(lgX-lgW,lgY,lgW,lgH);
    lgItems.forEach((item,i)=>{
      const lly=lgY+10+i*20;
      ctx.strokeStyle=item.color;ctx.lineWidth=item.dash?1.5:2;
      if(item.dash)ctx.setLineDash([4,4]);else ctx.setLineDash([]);
      ctx.beginPath();ctx.moveTo(lgX-lgW+8,lly);ctx.lineTo(lgX-lgW+30,lly);ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle='rgba(255,255,255,0.65)';ctx.textAlign='left';
      ctx.fillText(item.label,lgX-lgW+36,lly+3.5);
    });
  }

  // ── Master draw ───────────────────────────────────────────────────────────
  function drawAll(omg_norm,tr_sc,tr_uni,tr_grad,omge_u_norm,omge_g_norm,params){
    const W  =canvas.offsetWidth||600;
    const dpr=window.devicePixelRatio||1;
    const DW =W, DH=720;
    canvas.width=DW*dpr; canvas.height=DH*dpr;
    ctx.scale(dpr,dpr);
    ctx.clearRect(0,0,DW,DH);

    // Panel A: top 32% of canvas
    const PH_A = Math.round(DH*0.30);
    const topA  = 18;
    drawProfilePanel(omge_u_norm, omge_g_norm, params.omgt_r, params.delta, topA, PH_A, DW);

    // Subtle divider between panels
    const divY = topA + PH_A + 22;
    ctx.strokeStyle='rgba(255,255,255,0.06)'; ctx.lineWidth=1;
    ctx.setLineDash([4,6]);
    ctx.beginPath(); ctx.moveTo(ML,divY); ctx.lineTo(DW-MR,divY); ctx.stroke();
    ctx.setLineDash([]);

    // Panel B: remaining space
    const topB = topA + PH_A + 44;
    const PH_B = DH - topB - 48;
    drawTRPanel(omg_norm,tr_sc,tr_uni,tr_grad,params,topB,PH_B,DW,DH);
  }

  // ── Compute & schedule ────────────────────────────────────────────────────
  let debTimer=null;
  let params={p:1, delta:3, zeta:0.001, tau:500, omgt_r:35};

  function run(){
    const NF=500;
    const omg_arr=Array.from({length:NF},(_,i)=>(OMG_MIN+i*(OMG_MAX-OMG_MIN)/(NF-1))*omg1);
    const omg_norm=omg_arr.map(o=>o/omg1);

    // Normalized grading profiles (no units)
    const omge_u_norm = shuntFreqNorm(params.omgt_r, params.delta, 0);
    const omge_g_norm = shuntFreqNorm(params.omgt_r, params.delta, params.p);

    const tr_sc  =computeSC(params.zeta, omg_arr);
    const tr_uni =computeTR(omge_u_norm, params.tau, params.zeta, omg_arr, params.omgt_r);
    const tr_grad=computeTR(omge_g_norm, params.tau, params.zeta, omg_arr, params.omgt_r);

    drawAll(omg_norm,tr_sc,tr_uni,tr_grad,omge_u_norm,omge_g_norm,params);
  }

  function schedule(){clearTimeout(debTimer);debTimer=setTimeout(run,80);}

  function wire(id,valId,dec,key){
    const el=document.getElementById(id);
    if(!el) return;
    el.addEventListener('input',e=>{
      params[key]=+e.target.value;
      const d=document.getElementById(valId);
      if(d) d.textContent=(+e.target.value).toFixed(dec);
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