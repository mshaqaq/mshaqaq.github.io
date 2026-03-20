// ── iMSS Group — shared components ──
// Loads nav.html and footer.html into <nav> and <footer> on every page.
// To add a new page to the nav: edit nav.html only.

(async function () {
  async function loadComponent(selector, file) {
    const el = document.querySelector(selector);
    if (!el) return;
    try {
      const res  = await fetch(file);
      if (!res.ok) return;
      el.innerHTML = await res.text();

      // Highlight the active nav link
      // const path = window.location.pathname.split('/').pop() || 'index.html';
      // el.querySelectorAll('a[href]').forEach(a => {
      //   const href = a.getAttribute('href').split('/').pop();
      //   if (href === path) a.classList.add('active');
      // });
      const path = window.location.pathname.split('/').pop() || 'index.html';
      el.querySelectorAll('a[href]').forEach(a => {
      const href = a.getAttribute('href').split('/').pop();
      if (href === path) a.classList.add('active');
      });
    } catch (e) {
      console.warn('Could not load component:', file, e);
    }
  }

  loadComponent('nav',    'nav.html');
  loadComponent('footer', 'footer.html');
})();
