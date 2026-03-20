// ── iMSS Group — Dispersion Plot ──
// dispersion.js: handles the interactive piezoelectric metamaterial
// dispersion curve on the home page.
// Depends on: a <canvas id="disp-canvas"> element in the DOM.

(function() {

  const dc = document.getElementById('disp-canvas');
  if (!dc) return;
  const dx = dc.getContext('2d');

  let wt = 1.0, al = 0.15, ze = 0.0;

  // ── Complex number helpers ──
  function cdiv(a, b) {
    const d = b.re*b.re + b.im*b.im;
    return { re: (a.re*b.re + a.im*b.im)/d,
             im: (a.im*b.re - a.re*b.im)/d };
  }
  function csqrt(a) {
    const r = Math.sqrt(a.re*a.re + a.im*a.im);
    const t = Math.atan2(a.im, a.re);
    return { re: Math.sqrt(r)*Math.cos(t/2),
             im: Math.sqrt(r)*Math.sin(t/2) };
  }

  // ── Damped dispersion: sweep ω → complex k ──
  function dampedK(omega, wt, al, ze) {
    const eps = ze < 1e-9 ? 1e-9 : ze;
    const h = cdiv(
      { re: wt*wt, im: 0 },
      { re: wt*wt - omega*omega, im: 2*eps*wt*omega }
    );
    const k2 = { re: omega*omega*(1 + al*h.re),
                 im: omega*omega*al*h.im };
    return csqrt(k2);
  }

  // ── Main draw function ──
  function draw() {
    const W   = dc.offsetWidth || 600;
    const dpr = window.devicePixelRatio || 1;
    dc.width  = W   * dpr;
    dc.height = 420 * dpr;
    dx.scale(dpr, dpr);

    const DW = W, DH = 420;
    const P  = { t: 36, r: 36, b: 60, l: 68 };
    const PW = DW - P.l - P.r;
    const PH = DH - P.t - P.b;
    const maxK = 1, maxW = 2.4;

    dx.clearRect(0, 0, DW, DH);

    // ── Grid ──
    dx.strokeStyle = 'rgba(255,255,255,0.05)';
    dx.lineWidth = 1;
    for (let i = 0; i <= 4; i++) {
      const gx = P.l + (i/4)*PW;
      const gy = P.t + (i/4)*PH;
      dx.beginPath(); dx.moveTo(gx, P.t);   dx.lineTo(gx, P.t+PH); dx.stroke();
      dx.beginPath(); dx.moveTo(P.l, gy);   dx.lineTo(P.l+PW, gy); dx.stroke();
    }

    // ── Band gap shading ──
    // Lower edge: ωt       (lower branch asymptote as k→∞)
    // Upper edge: ωt√(1+α) (upper branch start at k=0)
    const bgLow  = wt;
    const bgHigh = wt * Math.sqrt(1 + al);

    const yBottom = P.t + PH - (bgLow  / maxW) * PH;
    const yTop    = P.t + PH - (bgHigh / maxW) * PH;

    dx.fillStyle = 'rgba(29,158,117,0.13)';
    dx.fillRect(P.l, yTop, PW, yBottom - yTop);

    dx.setLineDash([5, 5]);
    dx.strokeStyle = 'rgba(29,158,117,0.35)';
    dx.lineWidth = 1;
    [yTop, yBottom].forEach(y => {
      dx.beginPath(); dx.moveTo(P.l, y); dx.lineTo(P.l+PW, y); dx.stroke();
    });
    dx.setLineDash([]);

    dx.fillStyle  = 'rgba(29,158,117,0.7)';
    dx.font       = '12px Poppins, sans-serif';
    dx.textAlign  = 'right';
    dx.fillText('band gap', P.l + PW - 8, (yTop + yBottom)/2 + 4);

    // ── Short-circuit (reference) branch ──
    dx.beginPath();
    dx.strokeStyle = 'rgba(255,255,255,0.22)';
    dx.lineWidth   = 1.8;
    for (let i = 0; i <= 300; i++) {
      const k = (i/300) * maxK;
      const w = k * Math.PI;
      if (w > maxW) break;
      const x = P.l + (k/maxK)*PW;
      const y = P.t + PH - (w/maxW)*PH;
      i === 0 ? dx.moveTo(x, y) : dx.lineTo(x, y);
    }
    dx.stroke();

    // ── Damped dispersion curve ──
    let prevX = null, prevY = null;
    dx.lineWidth = 2.5;
    for (let i = 1; i <= 800; i++) {
      const omega  = (i/800) * maxW;
      const k      = dampedK(omega, wt, al, ze);
      const kNorm  = k.re / Math.PI;
      if (kNorm < 0 || kNorm > maxK) { prevX = null; continue; }
      const px   = P.l + (kNorm/maxK) * PW;
      const py   = P.t + PH - (omega/maxW) * PH;
      const fade = Math.max(0.08, 1 - Math.min(k.im * 5, 0.92));
      if (prevX !== null && Math.abs(py - prevY) < 25) {
        dx.beginPath();
        dx.strokeStyle = `rgba(201,168,76,${fade.toFixed(2)})`;
        dx.moveTo(prevX, prevY);
        dx.lineTo(px, py);
        dx.stroke();
      }
      prevX = px; prevY = py;
    }

    // ── Axes ──
    dx.strokeStyle = 'rgba(255,255,255,0.22)';
    dx.lineWidth   = 1;
    dx.beginPath();
    dx.moveTo(P.l, P.t);
    dx.lineTo(P.l, P.t + PH);
    dx.lineTo(P.l + PW, P.t + PH);
    dx.stroke();

    // ── Axis labels ──
    dx.fillStyle  = 'rgba(255,255,255,0.5)';
    dx.font       = '13px Poppins, sans-serif';
    dx.textAlign  = 'center';
    dx.fillText('Wave number  ka/\u03C0', P.l + PW/2, DH - 10);

    dx.save();
    dx.translate(16, P.t + PH/2);
    dx.rotate(-Math.PI / 2);
    dx.fillText('Frequency  \u03C9/\u03C9\u2080', 0, 0);
    dx.restore();

    // ── Tick labels ──
    dx.font = '12px Poppins, sans-serif';
    for (let i = 0; i <= 4; i++) {
      dx.textAlign = 'center';
      dx.fillText((i * 0.25).toFixed(2),
        P.l + (i/4)*PW, P.t + PH + 20);
      dx.textAlign = 'right';
      dx.fillText((i * maxW/4).toFixed(1),
        P.l - 10, P.t + PH - (i/4)*PH + 4);
    }
  }

  // ── Slider wiring ──
  function wire(id, valId, decimals, setter) {
    const el = document.getElementById(id);
    if (!el) return;
    el.addEventListener('input', e => {
      setter(+e.target.value);
      document.getElementById(valId).textContent =
        (+e.target.value).toFixed(decimals);
      draw();
    });
  }

  wire('disp-wt', 'disp-wt-val', 2, v => { wt = v; });
  wire('disp-al', 'disp-al-val', 2, v => { al = v; });
  wire('disp-ze', 'disp-ze-val', 3, v => { ze = v; });

  // ── Initial draw + resize handler ──
  draw();
  window.addEventListener('resize', draw);

})();
