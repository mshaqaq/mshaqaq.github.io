// ── iMSS Group — shared components ──
// Injects nav and footer directly (no fetch) for maximum mobile reliability.
// To add a new page: edit the navHTML and footerHTML strings below.

(function () {

  // ── Nav HTML ──────────────────────────────────────────────────────────────
  const navHTML = `
    <a class="nav-logo" href="/index.html">i<span>MSS</span> Group</a>
    <button class="nav-hamburger" id="nav-hamburger" aria-label="Open navigation">
      <span></span><span></span><span></span>
    </button>
    <ul class="nav-links" id="nav-links">
      <li><a href="/index.html">Home</a></li>
      <li><a href="/research.html">Research</a></li>
      <li><a href="/publications.html">Publications</a></li>
      <li><a href="/people.html">People</a></li>
      <li><a href="/teaching.html">Teaching</a></li>
      <li><a href="/interactive.html">Interactive</a></li>
      <li><a href="/contact.html">Contact</a></li>
    </ul>
  `;

  // ── Footer HTML ───────────────────────────────────────────────────────────
  const footerHTML = `
    <a class="footer-logo" href="/index.html">i<span>MSS</span> Group</a>
    <ul class="footer-links">
      <li><a href="/research.html">Research</a></li>
      <li><a href="/publications.html">Publications</a></li>
      <li><a href="/people.html">People</a></li>
      <li><a href="/interactive.html">Interactive</a></li>
      <li><a href="/contact.html">Contact</a></li>
    </ul>
    <p class="copy">© 2025 Mustafa AlShaqaq · KFUPM</p>
  `;

  // ── Inject into DOM ───────────────────────────────────────────────────────
  function inject(selector, html) {
    const el = document.querySelector(selector);
    if (!el) return;
    el.innerHTML = html;

    // Highlight active link
    const path = window.location.pathname.split('/').pop() || 'index.html';
    el.querySelectorAll('a[href]').forEach(a => {
      const href = a.getAttribute('href').split('/').pop();
      if (href === path) a.classList.add('active');
    });
  }

  // ── Hamburger wiring ──────────────────────────────────────────────────────
  function wireHamburger() {
    const btn   = document.getElementById('nav-hamburger');
    const links = document.getElementById('nav-links');
    if (!btn || !links) return;

    function openMenu()  { btn.classList.add('open');    links.classList.add('open');    btn.setAttribute('aria-label','Close navigation'); }
    function closeMenu() { btn.classList.remove('open'); links.classList.remove('open'); btn.setAttribute('aria-label','Open navigation'); }
    function toggle(e)   { e.stopPropagation(); btn.classList.contains('open') ? closeMenu() : openMenu(); }

    btn.addEventListener('click',    toggle);
    btn.addEventListener('touchend', function(e){ e.preventDefault(); toggle(e); });

    links.querySelectorAll('a').forEach(a => {
      a.addEventListener('click',    closeMenu);
      a.addEventListener('touchend', closeMenu);
    });

    document.addEventListener('click', e => {
      const nav = document.querySelector('nav');
      if (nav && !nav.contains(e.target)) closeMenu();
    });
  }

  // ── Run on DOM ready ──────────────────────────────────────────────────────
  function init() {
    inject('nav',    navHTML);
    inject('footer', footerHTML);
    wireHamburger();
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }

})();
