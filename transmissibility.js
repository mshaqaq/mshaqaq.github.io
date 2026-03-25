(function () {
  const canvas = document.getElementById('transCanvas');
  if (!canvas) return;

  const ctx = canvas.getContext('2d');
  const zetaInput = document.getElementById('transZeta');
  const deltaInput = document.getElementById('transDelta');
  const powerInput = document.getElementById('transPower');
  const tauInput = document.getElementById('transTau');
  const zetaVal = document.getElementById('transZetaVal');
  const deltaVal = document.getElementById('transDeltaVal');
  const powerVal = document.getElementById('transPowerVal');
  const tauVal = document.getElementById('transTauVal');

  const hs = 1e-4, bs = 1e-2, rho_s = 2700, cs = 69e9;
  const hp = 3e-4, rho_p = 7750, e31 = -12.3, c11 = 61e9, e33 = 13.3e-9;
  const L = 0.1, NM = 20, Q = 25, SZ = NM + Q;

  const m = 2 * rho_p * bs * hp + rho_s * bs * hs;
  const EI_b = ((2 * bs) / 3) * ((cs * hs ** 3) / 8 + c11 * ((hp + hs / 2) ** 3 - hs ** 3 / 8));
  const Cp_hat = (e33 * bs) / (2 * hp);
  const TH = ((e31 * bs) / (2 * hp)) * ((hp + hs / 2) ** 2 - hs ** 2 / 4);
  const sqrtmL = Math.sqrt(m * L);

  const xql = Array.from({ length: Q }, (_, k) => k / Q);
  const xqr = Array.from({ length: Q }, (_, k) => (k + 1) / Q);
  const Cpq = xql.map((_, k) => Cp_hat * L * (xqr[k] - xql[k]));
  const Cp = Cpq[0];

  const BETA = [1.875104068711961, 4.694091132974175, 7.854757438237613, 10.99554073487547, 14.13716839104647];
  function beta(n) { return n < 5 ? BETA[n] : ((2 * (n + 1) - 1) * Math.PI) / 2; }
  function omgnFn(n) { const b = beta(n); return b * b * Math.sqrt(EI_b / (m * L ** 4)); }
  function modeAndDeriv(n, x) {
    const B = beta(n);
    const sig = (Math.sin(B) - Math.sinh(B)) / (Math.cos(B) + Math.cosh(B));
    return {
      phi: Math.cos(B * x) - Math.cosh(B * x) + sig * (Math.sin(B * x) - Math.sinh(B * x)),
      dphi: (B / L) * (-Math.sin(B * x) - Math.sinh(B * x) + sig * (Math.cos(B * x) - Math.cosh(B * x))),
    };
  }
  function modeAvg(n) {
    const B = beta(n);
    return (2 * (Math.sin(B) - Math.sinh(B))) / (B * (Math.cos(B) + Math.cosh(B)));
  }

  const omgn = Array.from({ length: NM }, (_, n) => omgnFn(n));
  const omg1 = omgn[0];
  const phiTip = Array.from({ length: NM }, (_, n) => modeAndDeriv(n, 1).phi);
  const phiAvg = Array.from({ length: NM }, (_, n) => modeAvg(n));
  const coupling = Array.from({ length: Q }, (_, k) =>
    Array.from({ length: NM }, (_, n) => (TH * (modeAndDeriv(n, xqr[k]).dphi - modeAndDeriv(n, xql[k]).dphi)) / sqrtmL)
  );

  const A_re = new Float64Array(SZ * SZ);
  const A_im = new Float64Array(SZ * SZ);
  const b_re = new Float64Array(SZ);
  const b_im = new Float64Array(SZ);

  const LOG_MIN = -5, LOG_MAX = 3, OMG_MIN = 25, OMG_MAX = 45, THRESH = 0.1;
  const ML = 72, MR = 32;
  let W = 0, H = 0, dpr = 1, debTimer = null;
  const params = { p: 1.0, delta: 3.0, zeta: 0.001, tau: 500, omgt_r: 35 };

  function getColors() {
    const dark = document.body.classList.contains('dark');
    return dark ? {
      gridMajor: 'rgba(120,120,120,0.12)', gridMinor: 'rgba(120,120,120,0.06)', axis: 'rgba(190,190,190,0.8)', text: 'rgba(240,240,240,0.92)', textSoft: 'rgba(220,220,220,0.78)', textFaint: 'rgba(180,180,180,0.60)', shortCircuit: 'rgba(160,160,160,0.75)', uniform: 'rgba(245,245,245,0.92)', gradedAccent: 'rgba(245,230,120,0.95)', bandFill: 'rgba(245,230,120,0.16)', bandStroke: 'rgba(245,230,120,0.60)', panelBox: 'rgba(24,24,27,0.94)', panelStroke: 'rgba(82,82,91,0.60)', bg: '#09090b'
    } : {
      gridMajor: 'rgba(120,120,120,0.12)', gridMinor: 'rgba(120,120,120,0.06)', axis: 'rgba(70,70,70,0.8)', text: 'rgba(55,55,55,0.9)', textSoft: 'rgba(85,85,85,0.75)', textFaint: 'rgba(100,100,100,0.55)', shortCircuit: 'rgba(130,130,130,0.7)', uniform: 'rgba(20,20,20,0.9)', gradedAccent: 'rgba(210,180,60,0.95)', bandFill: 'rgba(255,235,140,0.35)', bandStroke: 'rgba(210,180,60,0.7)', panelBox: 'rgba(255,255,255,0.88)', panelStroke: 'rgba(170,170,170,0.45)', bg: '#ffffff'
    };
  }

  function initCanvasSize() {
    dpr = window.devicePixelRatio || 1;
    W = canvas.clientWidth || 700;
    H = Number(canvas.getAttribute('height')) || 720;
    canvas.width = Math.round(W * dpr);
    canvas.height = Math.round(H * dpr);
    canvas.style.width = '100%';
    canvas.style.height = H + 'px';
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.scale(dpr, dpr);
  }

  function shuntFreqNorm(omgt_r, delta, p) {
    return Array.from({ length: Q }, (_, k) => {
      if (delta === 0 || p === 0) return omgt_r;
      return omgt_r + delta - 2 * delta * (k / (Q - 1)) ** p;
    });
  }

  function solveGE() {
    const n = SZ;
    for (let col = 0; col < n; col++) {
      let maxV = 0, pivRow = col;
      for (let row = col; row < n; row++) {
        const r = A_re[row * n + col], i = A_im[row * n + col], v = r * r + i * i;
        if (v > maxV) { maxV = v; pivRow = row; }
      }
      if (pivRow !== col) {
        for (let j = 0; j < n; j++) {
          let t = A_re[col * n + j]; A_re[col * n + j] = A_re[pivRow * n + j]; A_re[pivRow * n + j] = t;
          t = A_im[col * n + j]; A_im[col * n + j] = A_im[pivRow * n + j]; A_im[pivRow * n + j] = t;
        }
        let t = b_re[col]; b_re[col] = b_re[pivRow]; b_re[pivRow] = t;
        t = b_im[col]; b_im[col] = b_im[pivRow]; b_im[pivRow] = t;
      }
      const pr = A_re[col * n + col], pi = A_im[col * n + col], pd = pr * pr + pi * pi;
      for (let row = col + 1; row < n; row++) {
        const fr = A_re[row * n + col], fi = A_im[row * n + col];
        const facR = (fr * pr + fi * pi) / pd, facI = (fi * pr - fr * pi) / pd;
        for (let j = col; j < n; j++) {
          const ar = A_re[col * n + j], ai = A_im[col * n + j];
          A_re[row * n + j] -= facR * ar - facI * ai;
          A_im[row * n + j] -= facR * ai + facI * ar;
        }
        const br = b_re[col], bi = b_im[col];
        b_re[row] -= facR * br - facI * bi;
        b_im[row] -= facR * bi + facI * br;
      }
    }
    for (let i = n - 1; i >= 0; i--) {
      for (let j = i + 1; j < n; j++) {
        const ar = A_re[i * n + j], ai = A_im[i * n + j];
        b_re[i] -= ar * b_re[j] - ai * b_im[j];
        b_im[i] -= ar * b_im[j] + ai * b_re[j];
      }
      const pr = A_re[i * n + i], pi = A_im[i * n + i], pd = pr * pr + pi * pi;
      const xr = b_re[i], xi = b_im[i];
      b_re[i] = (xr * pr + xi * pi) / pd;
      b_im[i] = (xi * pr - xr * pi) / pd;
    }
  }

  function computeAtFreq(omg, omge_norm, tau, zeta_r, omgt_r) {
    const omgt = omgt_r * omg1;
    const os2 = omg * omg;
    A_re.fill(0); A_im.fill(0); b_re.fill(0); b_im.fill(0);
    for (let n = 0; n < NM; n++) {
      const wn = omgn[n];
      A_re[n * SZ + n] = wn * wn - os2;
      A_im[n * SZ + n] = 2 * zeta_r * wn * omg;
      b_re[n] = sqrtmL * os2 * phiAvg[n];
    }
    for (let k = 0; k < Q; k++) {
      const we = omge_norm[k] * omg1;
      const rD = (Cpq[k] * omgt) / (tau * Cp);
      A_re[(NM + k) * SZ + (NM + k)] = rD;
      A_im[(NM + k) * SZ + (NM + k)] = omg - (we * we) / omg;
      for (let n = 0; n < NM; n++) {
        const c = coupling[k][n];
        A_re[n * SZ + (NM + k)] = -c;
        A_im[(NM + k) * SZ + n] = (omg * c) / Cp;
      }
    }
    solveGE();
    let dre = 1, dim = 0;
    for (let n = 0; n < NM; n++) {
      dre += (phiTip[n] * b_re[n]) / sqrtmL;
      dim += (phiTip[n] * b_im[n]) / sqrtmL;
    }
    return Math.sqrt(dre * dre + dim * dim);
  }

  function computeSCatFreq(omg, zeta_r) {
    const os2 = omg * omg;
    let dre = 1, dim = 0;
    for (let n = 0; n < NM; n++) {
      const wn = omgn[n], rhs = sqrtmL * os2 * phiAvg[n];
      const dr = wn * wn - os2, di = 2 * zeta_r * wn * omg, dsq = dr * dr + di * di;
      dre += (phiTip[n] * ((rhs * dr) / dsq)) / sqrtmL;
      dim += (phiTip[n] * ((-rhs * di) / dsq)) / sqrtmL;
    }
    return Math.sqrt(dre * dre + dim * dim);
  }

  function computeTR(omge_norm, tau, zeta_r, omg_arr, omgt_r) {
    const out = new Float64Array(omg_arr.length);
    for (let fi = 0; fi < omg_arr.length; fi++) out[fi] = computeAtFreq(omg_arr[fi], omge_norm, tau, zeta_r, omgt_r);
    return out;
  }

  function computeSC(zeta_r, omg_arr) {
    const out = new Float64Array(omg_arr.length);
    for (let fi = 0; fi < omg_arr.length; fi++) out[fi] = computeSCatFreq(omg_arr[fi], zeta_r);
    return out;
  }

  function computeBandgapWidth(omg_norm, tr_grad, threshold) {
    let best = 0, segStart = null;
    for (let i = 0; i < tr_grad.length; i++) {
      if (tr_grad[i] < threshold) {
        if (segStart === null) segStart = omg_norm[i];
      } else if (segStart !== null) {
        best = Math.max(best, omg_norm[i - 1] - segStart);
        segStart = null;
      }
    }
    if (segStart !== null) best = Math.max(best, omg_norm[omg_norm.length - 1] - segStart);
    return best;
  }

  function toPxTR(omg_norm, tr_val, top, PH, DW) {
    const px = ML + ((omg_norm - OMG_MIN) / (OMG_MAX - OMG_MIN)) * (DW - ML - MR);
    const logV = Math.log10(Math.max(tr_val, 1e-10));
    const py = top + (PH * (LOG_MAX - logV)) / (LOG_MAX - LOG_MIN);
    return { px, py };
  }

  function drawCurveTR(omg_norm, tr_vals, color, lw, dashed, top, PH, DW) {
    ctx.beginPath();
    ctx.strokeStyle = color; ctx.lineWidth = lw; ctx.setLineDash(dashed ? [7, 6] : []);
    let first = true;
    for (let i = 0; i < omg_norm.length; i++) {
      const { px, py } = toPxTR(omg_norm[i], tr_vals[i], top, PH, DW);
      if (py < top - 2 || py > top + PH + 2) { first = true; continue; }
      if (first) { ctx.moveTo(px, py); first = false; } else ctx.lineTo(px, py); }
    ctx.stroke(); ctx.setLineDash([]);
  }

  function drawShading(omg_norm, tr_grad, top, PH, DW, colors) {
    const { py: yThr } = toPxTR(0, THRESH, top, PH, DW);
    let inShade = false;
    for (let i = 0; i < omg_norm.length; i++) {
      const { px, py } = toPxTR(omg_norm[i], tr_grad[i], top, PH, DW);
      const cpy = Math.max(top, Math.min(top + PH, py));
      if (tr_grad[i] < THRESH) {
        if (!inShade) { ctx.beginPath(); ctx.moveTo(px, yThr); inShade = true; }
        ctx.lineTo(px, cpy);
      } else if (inShade) {
        ctx.lineTo(px, yThr); ctx.closePath(); ctx.fillStyle = colors.bandFill; ctx.fill(); inShade = false;
      }
    }
    if (inShade) {
      const { px } = toPxTR(omg_norm[omg_norm.length - 1], 0, top, PH, DW);
      ctx.lineTo(px, yThr); ctx.closePath(); ctx.fillStyle = colors.bandFill; ctx.fill();
    }
  }

  function drawProfilePanel(omge_u_norm, omge_g_norm, top, PH, DW, colors) {
    const PW = DW - ML - MR, L0 = ML, T0 = top;
    const pad = Math.max(params.delta * 0.2, 0.5);
    const yMin = params.omgt_r - params.delta - pad;
    const yMax = params.omgt_r + params.delta + pad;
    const yRng = yMax - yMin;
    const px = (k) => L0 + (k / (Q - 1)) * PW;
    const py = (val) => T0 + PH * (1 - (val - yMin) / yRng);

    for (let i = 0; i <= 4; i++) {
      const y = T0 + (i * PH) / 4;
      ctx.strokeStyle = colors.gridMajor; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(L0, y); ctx.lineTo(L0 + PW, y); ctx.stroke();
    }
    for (let k = 0; k < Q; k += 6) {
      ctx.strokeStyle = colors.gridMinor; ctx.beginPath(); ctx.moveTo(px(k), T0); ctx.lineTo(px(k), T0 + PH); ctx.stroke();
    }
    if (params.delta > 0) {
      const yHi = py(params.omgt_r + params.delta), yLo = py(params.omgt_r - params.delta);
      ctx.fillStyle = colors.bandFill; ctx.fillRect(L0, yHi, PW, yLo - yHi);
      ctx.setLineDash([4,4]); ctx.strokeStyle = colors.bandStroke; ctx.lineWidth = 0.75;
      [yHi, yLo].forEach((y) => { ctx.beginPath(); ctx.moveTo(L0, y); ctx.lineTo(L0 + PW, y); ctx.stroke(); });
      ctx.setLineDash([]);
    }
    ctx.strokeStyle = colors.shortCircuit; ctx.lineWidth = 2;
    const yU = py(omge_u_norm[0]);
    ctx.beginPath(); ctx.moveTo(L0, yU); ctx.lineTo(L0 + PW, yU); ctx.stroke();

    ctx.beginPath(); ctx.moveTo(px(0), yU);
    for (let k = 0; k < Q; k++) ctx.lineTo(px(k), py(omge_g_norm[k]));
    ctx.lineTo(px(Q - 1), yU); ctx.closePath();
    ctx.fillStyle = document.body.classList.contains('dark') ? 'rgba(245,230,120,0.12)' : 'rgba(255,235,140,0.22)';
    ctx.fill();

    ctx.strokeStyle = colors.gradedAccent; ctx.lineWidth = 2.2; ctx.beginPath();
    for (let k = 0; k < Q; k++) { const x = px(k), y = py(omge_g_norm[k]); k === 0 ? ctx.moveTo(x,y) : ctx.lineTo(x,y); }
    ctx.stroke();
    for (let k = 0; k < Q; k++) { ctx.beginPath(); ctx.arc(px(k), py(omge_g_norm[k]), 2.8, 0, Math.PI * 2); ctx.fillStyle = colors.gradedAccent; ctx.fill(); }

    ctx.strokeStyle = colors.axis; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(L0, T0); ctx.lineTo(L0, T0 + PH); ctx.lineTo(L0 + PW, T0 + PH); ctx.stroke();
    ctx.fillStyle = colors.textSoft; ctx.font = '13px Poppins, sans-serif'; ctx.textAlign = 'center'; ctx.fillText('Unit cell  k', L0 + PW / 2, T0 + PH + 20);
    ctx.font = '11px Poppins, sans-serif'; [1,7,13,19,25].forEach((k) => { ctx.fillStyle = colors.textFaint; ctx.fillText(k, px(k - 1), T0 + PH + 16); });
    ctx.save(); ctx.translate(14, T0 + PH / 2); ctx.rotate(-Math.PI / 2); ctx.font = '12px Poppins, sans-serif'; ctx.fillStyle = colors.textSoft; ctx.fillText('ωₜ,ₖ / ω₁', 0, 0); ctx.restore();
    ctx.font = '11px Poppins, sans-serif'; ctx.textAlign = 'right';
    const yTicks = params.delta > 0 ? [[params.omgt_r - params.delta, 'ωₜ−δ'], [params.omgt_r, 'ωₜ'], [params.omgt_r + params.delta, 'ωₜ+δ']] : [[params.omgt_r, 'ωₜ']];
    yTicks.forEach(([v, lbl]) => { ctx.fillStyle = colors.textFaint; ctx.fillText(lbl, L0 - 8, py(v) + 3); ctx.strokeStyle = colors.axis; ctx.lineWidth = 0.5; ctx.beginPath(); ctx.moveTo(L0 - 3, py(v)); ctx.lineTo(L0, py(v)); ctx.stroke(); });
    ctx.fillStyle = colors.textFaint; ctx.font = '10px Poppins, sans-serif'; ctx.textAlign = 'left'; ctx.fillText('(a)  Shunt frequency grading profile', L0 + 6, T0 + 13);
  }

  function drawTRPanel(omg_norm, tr_sc, tr_uni, tr_grad, top, PH, DW, DH, colors) {
    const PW = DW - ML - MR;
    for (let g = LOG_MIN; g <= LOG_MAX; g++) {
      const { py } = toPxTR(0, Math.pow(10, g), top, PH, DW);
      ctx.strokeStyle = colors.gridMajor; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(ML, py); ctx.lineTo(ML + PW, py); ctx.stroke();
      for (let mm = 2; mm <= 9; mm++) {
        const { py: pym } = toPxTR(0, mm * Math.pow(10, g), top, PH, DW);
        ctx.strokeStyle = colors.gridMinor; ctx.lineWidth = 0.5; ctx.beginPath(); ctx.moveTo(ML, pym); ctx.lineTo(ML + PW, pym); ctx.stroke();
      }
    }
    for (let g = OMG_MIN; g <= OMG_MAX; g += 5) {
      const px = ML + ((g - OMG_MIN) / (OMG_MAX - OMG_MIN)) * PW;
      ctx.strokeStyle = colors.gridMajor; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(px, top); ctx.lineTo(px, top + PH); ctx.stroke();
    }
    const { py: yThr } = toPxTR(0, THRESH, top, PH, DW);
    ctx.setLineDash([4, 4]); ctx.strokeStyle = colors.bandStroke; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(ML, yThr); ctx.lineTo(ML + PW, yThr); ctx.stroke(); ctx.setLineDash([]);
    ctx.fillStyle = colors.bandStroke; ctx.font = '10px Poppins, sans-serif'; ctx.textAlign = 'left'; ctx.fillText('|TR| = 0.1', ML + 5, yThr - 4);

    drawShading(omg_norm, tr_grad, top, PH, DW, colors);
    drawCurveTR(omg_norm, tr_sc, colors.shortCircuit, 1.5, true, top, PH, DW);
    drawCurveTR(omg_norm, tr_uni, colors.uniform, 2.0, false, top, PH, DW);
    drawCurveTR(omg_norm, tr_grad, colors.gradedAccent, 2.5, false, top, PH, DW);

    ctx.strokeStyle = colors.axis; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(ML, top); ctx.lineTo(ML, top + PH); ctx.lineTo(ML + PW, top + PH); ctx.stroke();
    ctx.fillStyle = colors.textSoft; ctx.font = '13px Poppins, sans-serif'; ctx.textAlign = 'center'; ctx.fillText('Frequency  ω/ω₁', ML + PW / 2, DH - 8);
    ctx.font = '11px Poppins, sans-serif'; for (let g = OMG_MIN; g <= OMG_MAX; g += 5) { const px = ML + ((g - OMG_MIN) / (OMG_MAX - OMG_MIN)) * PW; ctx.fillStyle = colors.textFaint; ctx.textAlign = 'center'; ctx.fillText(g, px, top + PH + 17); }
    ctx.save(); ctx.translate(14, top + PH / 2); ctx.rotate(-Math.PI / 2); ctx.font = '13px Poppins, sans-serif'; ctx.fillStyle = colors.textSoft; ctx.fillText('|TR| (tip transmissibility)', 0, 0); ctx.restore();
    const yL = ['10⁻⁵','10⁻⁴','10⁻³','10⁻²','10⁻¹','10⁰','10¹','10²','10³'];
    ctx.font = '11px Poppins, sans-serif'; ctx.textAlign = 'right';
    for (let g = LOG_MIN; g <= LOG_MAX; g++) { const { py } = toPxTR(0, Math.pow(10, g), top, PH, DW); ctx.fillStyle = colors.textFaint; ctx.fillText(yL[g - LOG_MIN], ML - 8, py + 3); }

    const bgW = computeBandgapWidth(omg_norm, tr_grad, THRESH);
    const bgTxt = bgW > 0 ? 'Bandgap  Δω/ω₁ = ' + bgW.toFixed(2) : 'Bandgap: discontinuous';
    ctx.font = '11px Poppins, sans-serif'; ctx.textAlign = 'left'; const tw = ctx.measureText(bgTxt).width;
    ctx.fillStyle = colors.panelBox; ctx.fillRect(ML + 6, top + 8, tw + 18, 22); ctx.strokeStyle = colors.panelStroke; ctx.lineWidth = 0.75; ctx.strokeRect(ML + 6, top + 8, tw + 18, 22); ctx.fillStyle = bgW > 0 ? colors.gradedAccent : colors.textFaint; ctx.fillText(bgTxt, ML + 15, top + 23);
    ctx.fillStyle = colors.textFaint; ctx.font = '10px Poppins, sans-serif'; ctx.textAlign = 'left'; ctx.fillText('(b)  Tip transmissibility', ML + 6, top + PH - 8);
  }

  function drawAll(omg_norm, tr_sc, tr_uni, tr_grad, omge_u_norm, omge_g_norm) {
    const colors = getColors();
    ctx.clearRect(0, 0, W, H); ctx.fillStyle = colors.bg; ctx.fillRect(0, 0, W, H);
    const DH = H, DW = W;
    const PH_A = Math.round(DH * 0.3), topA = 18;
    drawProfilePanel(omge_u_norm, omge_g_norm, topA, PH_A, DW, colors);
    const divY = topA + PH_A + 22; ctx.strokeStyle = colors.gridMinor; ctx.lineWidth = 1; ctx.setLineDash([4,6]); ctx.beginPath(); ctx.moveTo(ML, divY); ctx.lineTo(DW - MR, divY); ctx.stroke(); ctx.setLineDash([]);
    const topB = topA + PH_A + 44, PH_B = DH - topB - 48;
    drawTRPanel(omg_norm, tr_sc, tr_uni, tr_grad, topB, PH_B, DW, DH, colors);
  }

  function run() {
    params.omgt_r = 35;
    zetaVal.textContent = params.zeta.toFixed(4);
    deltaVal.textContent = params.delta.toFixed(1);
    powerVal.textContent = params.p.toFixed(2);
    tauVal.textContent = params.tau.toFixed(0);

    const NF = 500;
    const omgArr = Array.from({ length: NF }, (_, i) => (OMG_MIN + (i * (OMG_MAX - OMG_MIN)) / (NF - 1)) * omg1);
    const omgNorm = omgArr.map((o) => o / omg1);
    const omgeUNorm = shuntFreqNorm(params.omgt_r, params.delta, 0);
    const omgeGNorm = shuntFreqNorm(params.omgt_r, params.delta, params.p);
    const trSC = computeSC(params.zeta, omgArr);
    const trUni = computeTR(omgeUNorm, params.tau, params.zeta, omgArr, params.omgt_r);
    const trGrad = computeTR(omgeGNorm, params.tau, params.zeta, omgArr, params.omgt_r);
    drawAll(omgNorm, trSC, trUni, trGrad, omgeUNorm, omgeGNorm);
  }

  function schedule() { clearTimeout(debTimer); debTimer = setTimeout(run, 80); }
  function wire(id, valId, dec, key) {
    const el = document.getElementById(id); if (!el) return;
    el.addEventListener('input', (e) => {
      params[key] = +e.target.value;
      const d = document.getElementById(valId);
      if (d) d.textContent = (+e.target.value).toFixed(dec);
      schedule();
    });
  }
  wire('transZeta', 'transZetaVal', 4, 'zeta');
  wire('transPower', 'transPowerVal', 2, 'p');
  wire('transDelta', 'transDeltaVal', 1, 'delta');
  wire('transTau', 'transTauVal', 0, 'tau');

  window.drawTransmissibility = run;
  window.refreshTransmissibilityCanvas = function () { initCanvasSize(); run(); };

  initCanvasSize();
  run();
})();