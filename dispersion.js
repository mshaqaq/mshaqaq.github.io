(function () {
  const canvas = document.getElementById('dispersionCanvas');
  if (!canvas) return;

  const ctx = canvas.getContext('2d');
  const wtInput = document.getElementById('dispWt');
  const alphaInput = document.getElementById('dispAlpha');
  const zetaInput = document.getElementById('dispZeta');
  const wtVal = document.getElementById('dispWtVal');
  const alphaVal = document.getElementById('dispAlphaVal');
  const zetaVal = document.getElementById('dispZetaVal');

  let W = 0;
  let H = 0;
  let dpr = 1;

  function getColors() {
    const dark = document.body.classList.contains('dark');
    return dark
      ? {
          bg: '#09090b',
          grid: 'rgba(161,161,170,0.12)',
          axis: 'rgba(228,228,231,0.82)',
          text: 'rgba(244,244,245,0.96)',
          muted: 'rgba(212,212,216,0.74)',
          band: 'rgba(245,230,120,0.18)',
          bandStroke: 'rgba(245,230,120,0.50)',
          ref: 'rgba(161,161,170,0.72)',
          curve: 'rgba(255,255,255,0.94)',
        }
      : {
          bg: '#ffffff',
          grid: 'rgba(82,82,91,0.12)',
          axis: 'rgba(39,39,42,0.82)',
          text: 'rgba(24,24,27,0.95)',
          muted: 'rgba(82,82,91,0.75)',
          band: 'rgba(245,230,120,0.28)',
          bandStroke: 'rgba(185,160,55,0.64)',
          ref: 'rgba(113,113,122,0.72)',
          curve: 'rgba(0,0,0,0.92)',
        };
  }

  function cdiv(a, b) {
    const d = b.re * b.re + b.im * b.im;
    return {
      re: (a.re * b.re + a.im * b.im) / d,
      im: (a.im * b.re - a.re * b.im) / d,
    };
  }

  function csqrt(a) {
    const r = Math.hypot(a.re, a.im);
    const t = Math.atan2(a.im, a.re);
    return {
      re: Math.sqrt(r) * Math.cos(t / 2),
      im: Math.sqrt(r) * Math.sin(t / 2),
    };
  }

  function initCanvasSize() {
    dpr = window.devicePixelRatio || 1;
    W = canvas.clientWidth || 700;
    H = Number(canvas.getAttribute('height')) || 420;

    canvas.width = Math.round(W * dpr);
    canvas.height = Math.round(H * dpr);
    canvas.style.width = '100%';
    canvas.style.height = H + 'px';

    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.scale(dpr, dpr);
  }

  function drawDispersion() {
    const C = getColors();
    const wt = Number(wtInput.value);
    const alpha = Number(alphaInput.value);
    const zeta = Number(zetaInput.value);

    wtVal.textContent = wt.toFixed(2);
    alphaVal.textContent = alpha.toFixed(2);
    zetaVal.textContent = zeta.toFixed(3);

    const P = { l: 66, r: 28, t: 30, b: 56 };
    const PW = W - P.l - P.r;
    const PH = H - P.t - P.b;
    const maxK = 1.0;
    const maxOm = 2.4;

    const xMap = (kNorm) => P.l + (kNorm / maxK) * PW;
    const yMap = (om) => P.t + PH - (om / maxOm) * PH;

    ctx.clearRect(0, 0, W, H);
    ctx.fillStyle = C.bg;
    ctx.fillRect(0, 0, W, H);

    ctx.strokeStyle = C.grid;
    ctx.lineWidth = 1;
    for (let i = 0; i <= 4; i++) {
      const gx = P.l + (i / 4) * PW;
      const gy = P.t + (i / 4) * PH;
      ctx.beginPath();
      ctx.moveTo(gx, P.t);
      ctx.lineTo(gx, P.t + PH);
      ctx.stroke();
      ctx.beginPath();
      ctx.moveTo(P.l, gy);
      ctx.lineTo(P.l + PW, gy);
      ctx.stroke();
    }

    const bandLow = wt;
    const bandHigh = wt * Math.sqrt(1 + alpha);
    const yBottom = yMap(Math.min(bandLow, maxOm));
    const yTop = yMap(Math.min(bandHigh, maxOm));
    if (bandLow < maxOm) {
      ctx.fillStyle = C.band;
      ctx.fillRect(P.l, yTop, PW, Math.max(0, yBottom - yTop));
      ctx.setLineDash([5, 5]);
      ctx.strokeStyle = C.bandStroke;
      [yTop, yBottom].forEach((y) => {
        ctx.beginPath();
        ctx.moveTo(P.l, y);
        ctx.lineTo(P.l + PW, y);
        ctx.stroke();
      });
      ctx.setLineDash([]);
      ctx.fillStyle = C.muted;
      ctx.font = '12px Poppins, sans-serif';
      ctx.textAlign = 'right';
      ctx.fillText('band gap', P.l + PW - 8, (yTop + yBottom) / 2 + 4);
    }

    ctx.beginPath();
    ctx.setLineDash([8, 6]);
    ctx.strokeStyle = C.ref;
    ctx.lineWidth = 1.8;
    for (let i = 0; i <= 300; i++) {
      const kNorm = (i / 300) * maxK;
      const om = kNorm * Math.PI;
      if (om > maxOm) break;
      const x = xMap(kNorm);
      const y = yMap(om);
      i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
    }
    ctx.stroke();
    ctx.setLineDash([]);

    function dampedK(om) {
      const eps = Math.max(zeta, 1e-9);
      const h = cdiv(
        { re: wt * wt, im: 0 },
        { re: wt * wt - om * om, im: 2 * eps * wt * om }
      );
      const k2 = {
        re: om * om * (1 + alpha * h.re),
        im: om * om * alpha * h.im,
      };
      return csqrt(k2);
    }

    ctx.beginPath();
    ctx.strokeStyle = C.curve;
    ctx.lineWidth = 2.0;

    let penDown = false;
    let prevX = null;
    let prevY = null;

    for (let i = 1; i <= 1200; i++) {
      const om = (i / 1200) * maxOm;
      const k = dampedK(om);
      const kNorm = k.re / Math.PI;

      if (!Number.isFinite(kNorm) || !Number.isFinite(k.im) || kNorm < 0 || kNorm > maxK) {
        penDown = false;
        prevX = null;
        prevY = null;
        continue;
      }

      const x = xMap(kNorm);
      const y = yMap(om);
      if (!Number.isFinite(x) || !Number.isFinite(y)) {
        penDown = false;
        prevX = null;
        prevY = null;
        continue;
      }

      if (prevX !== null) {
        const dx = Math.abs(x - prevX);
        const dy = Math.abs(y - prevY);
        if (dx > 14 || dy > 10) {
          penDown = false;
        }
      }

      if (!penDown) {
        ctx.moveTo(x, y);
        penDown = true;
      } else {
        ctx.lineTo(x, y);
      }

      prevX = x;
      prevY = y;
    }
    ctx.stroke();

    ctx.strokeStyle = C.axis;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(P.l, P.t);
    ctx.lineTo(P.l, P.t + PH);
    ctx.lineTo(P.l + PW, P.t + PH);
    ctx.stroke();

    ctx.fillStyle = C.text;
    ctx.font = '13px Poppins, sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText('Wave number  ka/π', P.l + PW / 2, H - 14);

    ctx.save();
    ctx.translate(18, P.t + PH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText('Frequency  ω/ω₀', 0, 0);
    ctx.restore();

    ctx.font = '12px Poppins, sans-serif';
    for (let i = 0; i <= 4; i++) {
      ctx.textAlign = 'center';
      ctx.fillText((i * 0.25).toFixed(2), P.l + (i / 4) * PW, P.t + PH + 24);
      ctx.textAlign = 'right';
      ctx.fillText(((i * maxOm) / 4).toFixed(1), P.l - 10, P.t + PH - (i / 4) * PH + 4);
    }
  }

  window.drawDispersion = drawDispersion;
  window.refreshDispersionCanvas = function () {
    initCanvasSize();
    drawDispersion();
  };

  ['input', 'change'].forEach((evt) => {
    wtInput.addEventListener(evt, drawDispersion);
    alphaInput.addEventListener(evt, drawDispersion);
    zetaInput.addEventListener(evt, drawDispersion);
  });

  initCanvasSize();
  drawDispersion();
})();