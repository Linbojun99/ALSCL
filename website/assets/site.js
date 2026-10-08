(() => {
  'use strict';
  const zh = document.body.dataset.language === 'zh';
  const menu = document.querySelector('.menu-toggle');
  const articles = document.querySelector('.articles-dropdown');
  const closeArticles = () => { articles.open = false; };
  document.addEventListener('click', event => {
    if (!articles.contains(event.target)) closeArticles();
  });
  articles.addEventListener('keydown', event => {
    if (event.key === 'Escape' && articles.open) {
      event.preventDefault();
      closeArticles();
      articles.querySelector('summary').focus();
    }
  });
  articles.addEventListener('focusout', event => {
    if (!articles.contains(event.relatedTarget)) closeArticles();
  });
  articles.querySelectorAll('a').forEach(link => link.addEventListener('click', () => {
    closeArticles();
    menu.setAttribute('aria-expanded', 'false');
    document.getElementById('navigation').classList.remove('open');
  }));
  menu.addEventListener('click', () => {
    const expanded = menu.getAttribute('aria-expanded') === 'true';
    menu.setAttribute('aria-expanded', String(!expanded));
    document.getElementById('navigation').classList.toggle('open', !expanded);
    if (expanded) closeArticles();
  });
  const language = document.getElementById('language-switch');
  const updateLanguageHash = () => {
    const target = new URL(language.href);
    target.hash = window.location.hash;
    language.href = target.href;
  };
  updateLanguageHash();
  window.addEventListener('hashchange', updateLanguageHash);
  document.querySelectorAll('.copy').forEach(button => {
    button.addEventListener('click', async () => {
      try {
        await navigator.clipboard.writeText(button.nextElementSibling.innerText);
        button.textContent = zh ? '已复制' : 'Copied';
      } catch {
        button.textContent = zh ? '请选择代码复制' : 'Select code to copy';
      }
      setTimeout(() => { button.textContent = zh ? '复制' : 'Copy'; }, 1800);
    });
  });
  document.querySelectorAll('main table').forEach(table => {
    const wrapper = document.createElement('div');
    wrapper.className = 'table-scroll';
    wrapper.tabIndex = 0;
    wrapper.setAttribute('role', 'region');
    wrapper.setAttribute('aria-label', zh ? '可横向滚动的表格' : 'Scrollable table');
    table.before(wrapper); wrapper.append(table);
  });
  if (window.katex) {
    document.querySelectorAll('.math-display, .math-inline').forEach(node => {
      katex.render(node.textContent, node, { displayMode: node.classList.contains('math-display'), throwOnError: false, strict: 'ignore' });
    });
  }
  const tocLinks = Array.from(document.querySelectorAll('.toc a'));
  if ('IntersectionObserver' in window) {
    const observer = new IntersectionObserver(entries => {
      for (const entry of entries) {
        if (!entry.isIntersecting) continue;
        tocLinks.forEach(a => a.classList.toggle('active', a.hash === '#' + entry.target.id));
      }
    }, { rootMargin: '-80px 0px -65% 0px' });
    document.querySelectorAll('main h2, main h3').forEach(h => observer.observe(h));
  }
  const dialog = document.getElementById('search-dialog');
  const input = document.getElementById('search-input');
  const results = document.getElementById('search-results');
  const status = document.getElementById('search-status');
  let indexPromise;
  const index = () => indexPromise || (indexPromise = fetch(document.body.dataset.search)
    .then(response => { if (!response.ok) throw new Error('Search unavailable'); return response.json(); })
    .catch(error => { indexPromise = null; throw error; }));
  async function search() {
    const query = input.value.trim().toLowerCase();
    results.replaceChildren();
    if (!query) { status.textContent = zh ? '输入函数名、参数或主题。' : 'Enter a function, argument or topic.'; return; }
    try {
      const data = await index();
      if (query !== input.value.trim().toLowerCase()) return;
      const words = query.split(/\s+/);
      const ranked = data.map(page => {
        const title = page.title.toLowerCase(), text = page.text.toLowerCase();
        if (!words.every(word => title.includes(word) || text.includes(word))) return null;
        const score = (title === query || title === query + '()' ? 100 : 0) + words.reduce((n,w) => n + (title.includes(w) ? 15 : 1), 0);
        return { ...page, score };
      }).filter(Boolean).sort((a,b) => b.score - a.score).slice(0,20);
      status.textContent = ranked.length ? (zh ? `找到 ${ranked.length} 条相关结果` : `${ranked.length} relevant results`) : (zh ? '没有找到结果，请尝试其他关键词。' : 'No results. Try a function name or another topic.');
      const base = new URL(document.body.dataset.root || './', window.location.href);
      ranked.forEach(page => {
        const a = document.createElement('a'); a.href = new URL(page.url, base).href;
        const kind = document.createElement('small'); kind.textContent = page.type;
        const title = document.createElement('strong'); title.textContent = page.title;
        const excerpt = document.createElement('p');
        const position = Math.max(0, page.text.toLowerCase().indexOf(words[0]) - 45);
        excerpt.textContent = (position ? '…' : '') + page.text.slice(position, position + 170) + '…';
        a.append(kind,title,excerpt); results.append(a);
      });
    } catch { status.textContent = zh ? '搜索暂时无法载入，请使用函数参考或文章目录。' : 'Search could not load. Browse Reference or Articles instead.'; }
  }
  function openSearch() { dialog.showModal(); input.focus(); search(); }
  document.getElementById('search-open').addEventListener('click', openSearch);
  document.getElementById('search-close').addEventListener('click', () => dialog.close());
  input.addEventListener('input', search);
  input.addEventListener('keydown', event => {
    if (event.key === 'ArrowDown') { event.preventDefault(); results.querySelector('a')?.focus(); }
    if (event.key === 'Enter') results.querySelector('a')?.click();
  });
  dialog.addEventListener('click', event => { if (event.target === dialog) dialog.close(); });
  document.addEventListener('keydown', event => {
    if (event.key === '/' && !event.metaKey && !event.ctrlKey && !['INPUT','TEXTAREA','SELECT'].includes(document.activeElement.tagName) && !document.activeElement.isContentEditable && !dialog.open) {
      event.preventDefault(); openSearch();
    }
  });
})();
