const body = document.body;
const themeToggle = document.getElementById('themeToggle');
const menuToggle = document.getElementById('menuToggle');
const navLinks = document.getElementById('navLinks');
const savedTheme = localStorage.getItem('mustafa-theme');

if (savedTheme === 'light') body.classList.remove('dark');
updateThemeLabel();

function renderMath() {
  if (!window.renderMathInElement) return;

  window.renderMathInElement(document.body, {
    delimiters: [
      { left: '$$', right: '$$', display: true },
      { left: '$', right: '$', display: false },
      { left: '\\(', right: '\\)', display: false },
      { left: '\\[', right: '\\]', display: true }
    ],
    throwOnError: false,
    ignoredTags: ['script', 'noscript', 'style', 'textarea', 'pre', 'code']
  });
}

window.addEventListener('load', renderMath);

window.addEventListener('load', renderMath);

themeToggle.addEventListener('click', () => {
  body.classList.toggle('dark');
  localStorage.setItem('mustafa-theme', body.classList.contains('dark') ? 'dark' : 'light');
  updateThemeLabel();
  if (window.drawDispersion) window.drawDispersion();
  if (window.drawTransmissibility) window.drawTransmissibility();
});

menuToggle.addEventListener('click', () => {
  navLinks.classList.toggle('open');
});

navLinks.querySelectorAll('a').forEach(a => {
  a.addEventListener('click', () => navLinks.classList.remove('open'));
});

function updateThemeLabel() {
  const nextMode = body.classList.contains('dark') ? 'light' : 'dark';
  themeToggle.setAttribute('aria-label', `Switch to ${nextMode} mode`);
  themeToggle.setAttribute('title', `Switch to ${nextMode} mode`);
}

const pubSearch = document.getElementById('pubSearch');
const pubRecords = Array.from(document.querySelectorAll('.pub-record'));
pubSearch.addEventListener('input', (e) => {
  const q = e.target.value.trim().toLowerCase();
  pubRecords.forEach(item => {
    item.style.display = item.dataset.search.includes(q) ? '' : 'none';
  });
});

let resizeTimer = null;
window.addEventListener('resize', () => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(() => {
    if (window.refreshDispersionCanvas) {
      window.refreshDispersionCanvas();
    } else if (window.drawDispersion) {
      window.drawDispersion();
    }

    if (window.refreshTransmissibilityCanvas) {
      window.refreshTransmissibilityCanvas();
    } else if (window.drawTransmissibility) {
      window.drawTransmissibility();
    }
  }, 80);
});