// ── iMSS Group — shared components ──
// Loads nav.html and footer.html into <nav> and <footer> on every page.
// Handles hamburger menu toggle on mobile.
// To add a new page to the nav: edit nav.html only.

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

      // Wire hamburger toggle after nav is injected
      if (selector === 'nav') wireHamburger();

    } catch (e) {
      console.warn('Could not load component:', file, e);
    }
  }

  function wireHamburger() {
    const btn   = document.getElementById('nav-hamburger');
    const links = document.getElementById('nav-links');
    if (!btn || !links) return;

    btn.addEventListener('click', () => {
      btn.classList.toggle('open');
      links.classList.toggle('open');
    });

    // Close menu when any nav link is clicked
    links.querySelectorAll('a').forEach(a => {
      a.addEventListener('click', () => {
        btn.classList.remove('open');
        links.classList.remove('open');
      });
    });

    // Close menu when clicking outside the nav
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
