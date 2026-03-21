// ── iMSS Group — shared components ──
// Loads nav.html and footer.html, wires hamburger menu on mobile.

(async function () {

  async function loadComponent(selector, file) {
    const el = document.querySelector(selector);
    if (!el) return;
    try {
      const res = await fetch(file);
      if (!res.ok) return;
      el.innerHTML = await res.text();

      // Highlight active nav link
      const path = window.location.pathname.split('/').pop() || 'index.html';
      el.querySelectorAll('a[href]').forEach(a => {
        const href = a.getAttribute('href').split('/').pop();
        if (href === path) a.classList.add('active');
      });

      if (selector === 'nav') wireHamburger();

    } catch (e) {
      console.warn('Could not load component:', file, e);
    }
  }

  function wireHamburger() {
    const btn   = document.getElementById('nav-hamburger');
    const links = document.getElementById('nav-links');
    if (!btn || !links) return;

    // Use touchend for faster mobile response, fallback to click
    function toggle(e) {
      e.stopPropagation();
      btn.classList.toggle('open');
      links.classList.toggle('open');
    }

    btn.addEventListener('click', toggle);
    btn.addEventListener('touchend', function(e) {
      e.preventDefault();
      toggle(e);
    });

    // Close when a link is tapped
    links.querySelectorAll('a').forEach(a => {
      a.addEventListener('click', () => {
        btn.classList.remove('open');
        links.classList.remove('open');
      });
    });

    // Close when tapping outside the nav
    document.addEventListener('click', e => {
      const nav = document.querySelector('nav');
      if (nav && !nav.contains(e.target)) {
        btn.classList.remove('open');
        links.classList.remove('open');
      }
    });
  }

  loadComponent('nav',    'nav.html');
  loadComponent('footer', 'footer.html');

})();
