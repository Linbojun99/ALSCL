#!/usr/bin/env python3
"""Build the bilingual ALSCL manual. Run from any directory; no model fitting."""
from pathlib import Path
import argparse
import html
import hashlib
import json
import posixpath
import re
import shutil
import subprocess
import urllib.parse

import markdown
from bs4 import BeautifulSoup
from pygments.formatters import HtmlFormatter
from vendor import prepare_katex

ROOT = Path(__file__).resolve().parents[1]
WEB = ROOT / 'website'
REPO = 'https://github.com/Linbojun99/ALSCL'
BASE = 'https://linbojun99.github.io/ALSCL/'
LOGO_VERSION = hashlib.sha256((ROOT / 'ALSCLlogo.png').read_bytes()).hexdigest()[:12]
VERSION = re.search(r'^Version: (.+)$', (ROOT / 'DESCRIPTION').read_text(), re.M)[1]
ARTICLES = json.loads((WEB / 'articles.json').read_text())
ARTICLE_GROUPS = json.loads((WEB / 'article-groups.json').read_text())
ARTICLES_BY_SLUG = {article['slug']: article for article in ARTICLES}
CASES = json.loads((WEB / 'cases.json').read_text())
UI = {
    'en': dict(reference='Reference', articles='Articles', start='Get started', news='Changelog',
               search='Search documentation', contents='On this page', source='View source',
               usage='Usage', arguments='Arguments', returns='Value', details='Details', examples='Examples',
               related='See also', argument='Argument', default='Default', description='Description',
               required='required', home='Home', copy='Copy', copied='Copied', previous='Previous', next='Next',
               skip='Skip to contents', menu='Menu', noresults='No results. Try a function name or topic.',
               close='Close', functions='Function reference', setup='Example setup'),
    'zh': dict(reference='函数参考', articles='专题案例', start='开始使用', news='更新日志',
               search='搜索文档', contents='本页目录', source='查看源代码',
               usage='用法', arguments='参数', returns='返回值', details='使用说明', examples='示例',
               related='相关内容', argument='参数', default='默认值', description='说明',
               required='必填', home='首页', copy='复制', copied='已复制', previous='上一篇', next='下一篇',
               skip='跳转到正文', menu='菜单', noresults='没有找到结果，请尝试函数名或其他关键词。',
               close='关闭', functions='函数参考', setup='示例准备')}

GROUPS = [
    ('fit', 'Fit models', '模型拟合', ['run_acl', 'run_alscl', 'create_parameters', 'create_parameters_alscl', 'generate_map']),
    ('biology', 'Biology and simulation', '生物参数与模拟', ['VB_func', 'mat_func', 'initialize_params', 'sim_cal', 'sim_data', 'simulate_example_data', 'sim_acl']),
    ('diagnostics', 'Diagnostics and retrospectives', '诊断与回溯分析', ['diagnose_model', 'diagnostic_metrics', 'compare_models', 'retro_model', 'retro_acl', 'retro_alscl', 'plot_retro']),
    ('population', 'Survey and population plots', '调查与种群绘图', ['plot_CatL', 'plot_abundance', 'plot_biomass', 'plot_SSB', 'plot_catch', 'plot_recruitment', 'plot_SSB_Rec']),
    ('process', 'Growth, mortality and residuals', '生长、死亡率与残差', ['plot_VB', 'plot_pla', 'plot_fishing_mortality', 'plot_residuals', 'plot_deviance', 'plot_ridges']),
    ('comparison', 'Compare models visually', '模型比较绘图', ['plot_compare_ts', 'plot_compare_CatL', 'plot_compare_F', 'plot_compare_annual_F', 'plot_compare_growth', 'plot_compare_selectivity', 'plot_compare_residuals', 'plot_compare_metrics']),
    ('theme', 'Plot appearance', '绘图样式', ['acl_theme', 'acl_theme_set', 'acl_theme_reset'])]

def md(text):
    # Preserve TeX before Markdown consumes backslashes, underscores and braces.
    equations = []
    def equation(match, display):
        token = f'ALSCLMATHTOKEN{len(equations)}END'
        equations.append((token, '<div class="math-display">' if display else '<span class="math-inline">', match[1], display))
        return '\n\n' + token + '\n\n' if display else token
    text = re.sub(r'```math\s*\n(.*?)\n```', lambda m: equation(m, True), text, flags=re.S)
    text = re.sub(r'\$`(.*?)`\$', lambda m: equation(m, False), text)
    result = markdown.markdown(text, extensions=['tables', 'fenced_code', 'codehilite', 'sane_lists'],
                               extension_configs={'codehilite': {'guess_lang': False, 'css_class': 'highlight'}})
    for token, opening, formula, display in equations:
        node = opening + html.escape(formula.strip()) + ('</div>' if display else '</span>')
        result = result.replace('<p>' + token + '</p>', node).replace(token, node)
    return result

def choose(label, lang):
    for sep in (' / ', ' · '):
        if sep in label:
            left, right = label.split(sep, 1)
            if re.search('[\u3400-\u9fff]', left) and not re.search('[\u3400-\u9fff]', right):
                return left if lang == 'zh' else right
    return label

def code_language(code, lang):
    lines = []
    for line in code.splitlines():
        if '#' in line:
            left, comment = line.split('#', 1)
            line = left + '#' + choose(comment, lang)
        lines.append(line)
    return '\n'.join(lines)

def rel(route, target):
    return posixpath.relpath(target, posixpath.dirname(route) or '.')

def lang_route(lang, route):
    return ('zh/' if lang == 'zh' else '') + route

def link(route, target, text, **attrs):
    extra = ''.join(f' {k.replace("_", "-")}="{html.escape(v, quote=True)}"' for k, v in attrs.items())
    return f'<a href="{html.escape(rel(route, target), quote=True)}"{extra}>{text}</a>'

def orcid_link(route, name, identifier):
    label = html.escape(f'{name} ORCID: {identifier}', quote=True)
    icon = html.escape(rel(route, 'assets/orcid.svg'), quote=True)
    return (f'<a class="orcid-link" href="https://orcid.org/{identifier}" '
            f'aria-label="{label}" title="{label}">'
            f'<img src="{icon}" width="16" height="16" alt="ORCID iD" aria-hidden="true"></a>')

def article_sections(lang, home=False):
    sections = []
    for group in ARTICLE_GROUPS:
        tag = 'h2'
        anchor = 'guide-' + group['id'] if home else group['id']
        items = ''
        for slug in group['articles']:
            article = ARTICLES_BY_SLUG[slug]
            detail = 'Full guide and examples' if lang=='en' else '详细说明与示例'
            items += f'<h3 id="basic-{slug}">{html.escape(article["title"][lang])}</h3><p>{html.escape(article["description"][lang])}</p>'
            items += '<p><a href="' + ('articles/' if home else '') + slug + '.html">' + detail + ' →</a></p>'
        sections.append(f'<section class="article-group"><{tag} id="{anchor}">{html.escape(group["title"][lang])}</{tag}>'
                        f'<p>{html.escape(group["description"][lang])}</p>{items}</section>')
    return '\n'.join(sections)

def case_groups():
    return list(dict.fromkeys(case['group'] for case in CASES))

def rewrite_links(soup, route, lang):
    by_wiki = {Path(a['source']).stem: a['slug'] for a in ARTICLES}
    known_docs = {'USER_GUIDE.md': 'articles/first-model.html', 'PLOT_GALLERY.md': 'articles/population-plots.html',
                  'FUNCTION_REFERENCE.md': 'reference/index.html', 'CASE_STUDIES.md': 'articles/case-studies.html',
                  'PLOT_STYLE.md': 'articles/plot-customization.html', 'ADVANCED_WORKFLOWS.md': 'articles/reproducibility.html',
                  'BOOK_THEORY_NOTES.md': 'articles/model-theory.html', 'MATHEMATICAL_FRAMEWORK.md': 'articles/model-theory.html'}
    section_routes = {'theory':'model-theory','builtin':'first-model','simulation':'simulation','excel':'data-preparation',
                      'observations':'data-preparation','fitting':'fitting-controls','diagnostics':'diagnostics',
                      'plotting':'population-plots','retrospective':'retrospectives'}
    for tag in soup.select('[href], [src]'):
        attr = 'src' if tag.name == 'img' else 'href'
        if not tag.has_attr(attr): continue
        url = tag[attr]
        parsed = urllib.parse.urlsplit(url)
        path = parsed.path
        target = None
        # Repository assets are bundled once and reused by both languages.
        match = re.search(r'(?:^|/)(docs/)?(figures/[^?#]+|data/[^?#]+\.(?:xlsx|csv|rds))$', path)
        if match and (not parsed.netloc or parsed.netloc in ('github.com', 'raw.githubusercontent.com')):
            target = 'assets/manual/' + match[2]
        elif '/wiki/' in path and path.rsplit('/',1)[-1] in by_wiki:
            target = lang_route(lang, 'articles/' + by_wiki[path.rsplit('/',1)[-1]] + '.html')
        elif url == REPO + '/wiki': target = lang_route(lang, 'articles/index.html')
        elif path.endswith('FUNCTION_REFERENCE.md'):
            if parsed.fragment in FUNCTIONS: target = lang_route(lang, f'reference/{parsed.fragment}.html')
            elif parsed.fragment == 'estimation_parameters': target = lang_route(lang, 'articles/estimation-parameters.html')
            else: target = lang_route(lang, 'reference/index.html')
        elif path.endswith('USER_GUIDE.md') and parsed.fragment in section_routes:
            target = lang_route(lang, 'articles/' + section_routes[parsed.fragment] + '.html')
        elif path.rsplit('/',1)[-1] in known_docs:
            target = lang_route(lang, known_docs[path.rsplit('/',1)[-1]])
        elif path in ('../README.md', '/Linbojun99/ALSCL/blob/main/README.md'):
            target = lang_route(lang, 'index.html')
        if target: tag[attr] = rel(route, target)
        elif not parsed.scheme and not url.startswith('#') and path and not path.endswith('.html'):
            # Remaining source-file links refer to versioned repository content.
            source = posixpath.normpath(posixpath.join('docs', path))
            tag[attr] = REPO + '/blob/main/' + source + ('#'+parsed.fragment if parsed.fragment else '')
        if tag.name == 'img':
            tag['loading'] = 'lazy'; tag['decoding'] = 'async'
    return soup

def render_page(out, route_base, lang, title, body, source='', kind='article', description=''):
    route = lang_route(lang, route_base)
    ui = UI[lang]
    soup = BeautifulSoup(body, 'html.parser')
    # Stable bilingual section IDs allow language switching to preserve position.
    headings = soup.select('h2, h3')
    toc = []
    for i, heading in enumerate(headings):
        if not heading.has_attr('id'): heading['id'] = f'section-{i+1}'
        toc.append(f'<li class="level-{heading.name}"><a href="#{heading["id"]}">{html.escape(heading.get_text(" ", strip=True))}</a></li>')
    soup = rewrite_links(soup, route, lang)
    for code in soup.select('code'):
        if code.find_parent('pre'): continue
        name = code.get_text().removesuffix('()')
        if name in FUNCTIONS and not code.find_parent('a'):
            a = soup.new_tag('a', href=rel(route, lang_route(lang, f'reference/{name}.html')))
            code.wrap(a)
    for pre in soup.select('pre'):
        button = soup.new_tag('button', attrs={'type':'button','class':'copy','aria-label':ui['copy']+' R code'})
        button.string = ui['copy']
        pre.insert_before(button)
    plain = soup.get_text(' ', strip=True)
    SEARCH[lang].append({'title':title, 'url':route_base, 'type':ui['reference'] if kind=='reference' else (ui['articles'] if kind=='case' else ('Basic functions' if lang=='en' else '基本功能')), 'text':plain})
    nav = ''.join(link(route, lang_route(lang, p), ui[k], **({'aria_current':'page'} if route_base==p else {})) for k,p in
                  [('start','articles/getting-started.html'),('reference','reference/index.html')])
    category_links = ''
    for group in case_groups():
        cases = [case for case in CASES if case['group']==group]
        category_links += '<li class="dropdown-heading">' + html.escape(cases[0]['group_title'][lang]) + '</li>'
        category_links += ''.join('<li>' + link(route, lang_route(lang, 'articles/' + case['slug'] + '.html'),
                                              html.escape(case['title'][lang])) + '</li>' for case in cases)
    all_articles = link(route, lang_route(lang, 'articles/index.html'), 'All case studies' if lang=='en' else '全部专题案例')
    active = ' is-active' if kind=='case' else ''
    nav += (f'<details class="articles-dropdown{active}"><summary>{ui["articles"]}<span class="caret" aria-hidden="true"></span></summary>'
            f'<div class="articles-menu"><ul>{category_links}</ul><div class="all-articles">{all_articles}</div></div></details>')
    nav += link(route, lang_route(lang, 'news/index.html'), ui['news'], **({'aria_current':'page'} if route_base=='news/index.html' else {}))
    other = 'zh' if lang=='en' else 'en'
    switch = link(route, lang_route(other, route_base), '简体中文' if lang=='en' else 'English', id='language-switch', hreflang='zh-Hans' if other=='zh' else 'en', lang='zh-Hans' if other=='zh' else 'en')
    sidebar = '<nav class="toc" aria-label="'+ui['contents']+'"><h2>'+ui['contents']+'</h2><ul>'+''.join(toc)+'</ul></nav>'
    if kind=='home':
        sidebar = f'''<section><h2>{'Links' if lang=='en' else '相关链接'}</h2>
        <p><a href="{REPO}">{'Browse source code' if lang=='en' else '浏览源代码'}</a></p>
        <p><a href="{REPO}/issues">{'Report a bug' if lang=='en' else '问题反馈'}</a></p>
        <p>{link(route, lang_route(lang,'articles/data-preparation.html'), 'Data and Excel templates' if lang=='en' else '数据与 Excel 模板')}</p></section>
        <section><h2>{'License' if lang=='en' else '许可证'}</h2><a href="{REPO}/blob/main/LICENSE">GPL-3</a></section>
        <section><h2>{'Citation' if lang=='en' else '引用'}</h2><p><a href="#citation">Zhang &amp; Cadigan (2022)</a></p><p><a href="#krill-paper">Dong et al. (2025)</a></p></section>
        <section><h2>{'Developers' if lang=='en' else '开发者'}</h2>
        <p><span class="author-name">Hongyu Lin {orcid_link(route, 'Hongyu Lin', '0000-0001-6226-179X')}</span><small>{'Author, maintainer' if lang=='en' else '作者、维护者'}</small></p>
        <p><span class="author-name">Fan Zhang {orcid_link(route, 'Fan Zhang', '0000-0003-0214-1790')}</span><small>{'Author' if lang=='en' else '作者'}</small></p>
        <p>Sisong Dong<small>{'Author, contributor' if lang=='en' else '作者、贡献者'}</small></p></section>
        <section class="dev-status"><h2>{'Dev status' if lang=='en' else '开发状态'}</h2>
        <p><a href="{REPO}/actions/workflows/R-CMD-check.yaml?query=branch%3Amain"><img src="{REPO}/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main" alt="{'R package check status' if lang=='en' else 'R 包检查状态'}" height="20"></a></p>
        <p><a href="{REPO}/actions/workflows/documentation.yaml?query=branch%3Amain"><img src="{REPO}/actions/workflows/documentation.yaml/badge.svg?branch=main" alt="{'Documentation build and deployment status' if lang=='en' else '文档构建与部署状态'}" height="20"></a></p></section>'''
    source_link = f'<a class="source-link" href="{REPO}/blob/main/{source}">{ui["source"]} ↗</a>' if source else ''
    description = description or plain[:180]
    canonical = BASE + route
    logo = f'<img class="package-logo" src="{rel(route,"assets/ALSCLlogo.png")}?v={LOGO_VERSION}" alt="ALSCL" width="140" height="140">' if kind=='home' else ''
    section_label = ui['reference'] if kind=='reference' else (ui['articles'] if kind=='case' else ('Basic functions' if lang=='en' else '基本功能'))
    breadcrumb = '' if kind=='home' else f'<div class="breadcrumb">{link(route,lang_route(lang,"index.html"),ui["home"])} <span>/</span> {section_label}</div>'
    document = f'''<!doctype html>
<html lang="{'en' if lang=='en' else 'zh-Hans'}"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>{html.escape(title)} • ALSCL</title>
<meta name="description" content="{html.escape(description, quote=True)}"><meta name="theme-color" content="#f7f8fa">
<link rel="canonical" href="{canonical}"><link rel="alternate" hreflang="en" href="{BASE+route_base}">
<link rel="alternate" hreflang="zh-Hans" href="{BASE+'zh/'+route_base}"><link rel="alternate" hreflang="x-default" href="{BASE+route_base}">
<link rel="icon" href="{rel(route,'assets/ALSCLlogo.png')}?v={LOGO_VERSION}">
<link rel="stylesheet" href="{rel(route,'assets/site.css')}"><link rel="stylesheet" href="{rel(route,'assets/syntax.css')}">
<link rel="stylesheet" href="{rel(route,'assets/katex/katex.min.css')}">
<script defer src="{rel(route,'assets/katex/katex.min.js')}"></script>
<script defer src="{rel(route,'assets/site.js')}"></script></head>
<body data-language="{lang}" data-root="{rel(route,lang_route(lang,'index.html')).removesuffix('index.html')}" data-search="{rel(route,lang_route(lang,'search.json'))}">
<a class="skip-link" href="#main">{ui['skip']}</a>
<header class="navbar"><nav class="nav-inner" aria-label="{'Site navigation' if lang=='en' else '网站导航'}">
<div class="brand">{link(route,lang_route(lang,'index.html'),'ALSCL')}<span>{VERSION}</span></div>
<button class="menu-toggle" aria-expanded="false" aria-controls="navigation">{ui['menu']} ☰</button>
<div id="navigation" class="navigation">{nav}</div>
<div class="nav-tools"><button id="search-open" aria-label="{ui['search']}"><span>⌕</span> {'Search' if lang=='en' else '搜索'} <kbd>/</kbd></button>{switch}<a class="github" href="{REPO}" aria-label="GitHub"><svg viewBox="0 0 24 24" aria-hidden="true"><path d="M12 .8a11.2 11.2 0 0 0-3.54 21.83c.56.1.77-.24.77-.54v-2.1c-3.12.68-3.78-1.33-3.78-1.33-.51-1.29-1.24-1.64-1.24-1.64-1.02-.7.08-.69.08-.69 1.13.08 1.72 1.16 1.72 1.16 1 1.72 2.63 1.22 3.27.93.1-.72.39-1.22.71-1.5-2.49-.28-5.11-1.25-5.11-5.54 0-1.23.44-2.23 1.16-3.02-.12-.28-.5-1.43.11-2.98 0 0 .95-.3 3.08 1.15a10.72 10.72 0 0 1 5.6 0c2.14-1.45 3.08-1.15 3.08-1.15.61 1.55.23 2.7.12 2.98.72.79 1.15 1.79 1.15 3.02 0 4.3-2.62 5.25-5.12 5.53.4.35.76 1.03.76 2.08v3.1c0 .3.2.65.77.54A11.2 11.2 0 0 0 12 .8Z"/></svg></a></div></nav></header>
<div class="layout"><main id="main">{breadcrumb}<div class="page-header">{logo}<h1>{html.escape(title)}</h1>{source_link}</div>
{soup}</main><aside>{sidebar}</aside></div>
<footer><span>{'Developed by' if lang=='en' else '开发者'} <span class="author-name">Hongyu Lin {orcid_link(route, "Hongyu Lin", "0000-0001-6226-179X")}</span>, <span class="author-name">Fan Zhang {orcid_link(route, "Fan Zhang", "0000-0003-0214-1790")}</span>, Sisong Dong.</span><span>ALSCL {VERSION} · <a href="{REPO}">{'Source & documentation' if lang=='en' else '源代码与文档'}</a></span></footer>
<dialog id="search-dialog" aria-labelledby="search-title"><div class="search-head"><h2 id="search-title">{ui['search']}</h2><button id="search-close" aria-label="{ui['close']}">×</button></div>
<label class="sr-only" for="search-input">{ui['search']}</label><input id="search-input" type="search" autocomplete="off" placeholder="{'Function, argument or topic…' if lang=='en' else '函数、参数或主题…'}">
<p id="search-status" aria-live="polite"></p><div id="search-results"></div></dialog>
</body></html>'''
    path = out / route; path.parent.mkdir(parents=True, exist_ok=True); path.write_text(document)

def help_sections(name, meta):
    path = WEB / '.build/api' / (name + '.Rd.html')
    if not path.exists():
        candidates = list((WEB / '.build/api').glob('*.html'))
        path = next((p for p in candidates if re.search(r'\b'+re.escape(name)+r'\(', p.read_text())), None)
    if not path: return {}
    soup = BeautifulSoup(path.read_text(), 'html.parser')
    result = {}
    for heading in soup.select('h3'):
        nodes = []
        for node in heading.next_siblings:
            if getattr(node,'name',None) in ('h2','h3'): break
            nodes.append(str(node))
        result[heading.get_text(strip=True)] = ''.join(nodes).strip()
    return result

def load_reference():
    pieces = re.split(r'<a id="([^"]+)"></a>', (ROOT / 'docs/FUNCTION_REFERENCE.md').read_text())
    blocks = dict(zip(pieces[1::2],pieces[2::2]))
    extracted = json.loads((WEB / '.build/api/functions.json').read_text())
    reference = {}
    for entry in extracted:
        name = entry['name']; block = blocks[name]
        entry['summary'] = re.search(r'## `[^`]+`\s*\n\n([^\n]+)', block)[1]
        rows = []
        for line in block.splitlines():
            if line.startswith('| `'):
                cols = [c.strip() for c in re.split(r'(?<!\\)\|',line)[1:-1]]
                if len(cols)==4: rows.append(cols)
        args = {row[0].strip('`'):row for row in rows}
        defaults = entry['defaults'] or {}
        if set(args) != set(defaults): raise ValueError(f'Argument documentation differs from R source for {name}: {set(args)^set(defaults)}')
        entry['arguments'] = args
        entry['help'] = help_sections(name,entry)
        example_part = block.split('**示例 / Example:**',1)[-1]
        entry['example'] = re.search(r'```r\n(.*?)```',example_part,re.S)[1]
        reference[name] = entry
    return reference

def example_setup(name, lang):
    if name in ['VB_func','mat_func','initialize_params','create_parameters','create_parameters_alscl','generate_map','acl_theme','acl_theme_set','acl_theme_reset','simulate_example_data']: return 'library(ALSCL)'
    base = 'library(ALSCL)\nlibrary(ggplot2)\nlibrary(patchwork)\ndata("YTF_example")\nx <- YTF_example\ninputs <- x[c("data.CatL", "data.wgt", "data.mat")]\ndat <- x$data.CatL'
    if name in ['sim_cal','sim_data','sim_acl']:
        base += '\npa <- initialize_params(species = "flatfish")\nbio <- sim_cal(pa)'
        if name == 'sim_acl': base += '\nsm <- sim_data(bio, pa, iter_range = 4:5, return_iter = 4,\n               output_dir = "simulation_examples/flatfish")'
    elif name not in ('run_acl','run_alscl'):
        base += '\na <- do.call(run_acl, c(inputs, x$fit_args, x$fit_config$acl,\n                       list(silent = TRUE)))\nb <- do.call(run_alscl, c(inputs, x$fit_args, x$fit_config$alscl,\n                         list(silent = TRUE)))'
        if name=='plot_retro': base += '\nra <- do.call(retro_acl, c(inputs, x$fit_args, x$fit_config$acl,\n                          list(nyear = 3, silent = TRUE)))\nrb <- do.call(retro_alscl, c(inputs, x$fit_args, x$fit_config$alscl,\n                            list(nyear = 3, silent = TRUE)))'
    return base

def reference_body(name, entry, lang, supplements):
    ui = UI[lang]; extra = supplements[name]
    summary = choose(entry['summary'],lang)
    description = extra.get('description',{}).get(lang,summary)
    body = md(description)
    arguments=[arg + (' = '+default if default else '') for arg,default in (entry['defaults'] or {}).items()]
    signature=name+'('+', '.join(arguments)+')'
    if len(signature)>78: signature=name+'(\n  '+',\n  '.join(arguments)+'\n)'
    body += f'<h2 id="usage">{ui["usage"]}</h2>' + md('```r\n'+signature+'\n```')
    body += f'<h2 id="arguments">{ui["arguments"]}</h2>'
    if not entry['arguments']: body += '<p>'+('This function takes no arguments.' if lang=='en' else '此函数不接受参数。')+'</p>'
    else:
        body += f'<table class="arguments"><thead><tr><th>{ui["argument"]}</th><th>{ui["default"]}</th><th>{ui["description"]}</th></tr></thead><tbody>'
        for arg,row in entry['arguments'].items():
            default=entry['defaults'][arg] or ('…' if arg=='...' else ui['required'])
            explanation=row[3 if lang=='en' else 2]
            if arg=='ncores' and lang=='zh':
                explanation = '独立起点拟合的最大并行进程数，实际不超过 nstarts；不是单个 TMB 拟合的线程数。' if name.startswith('run_') else '外层并行进程数，用于模拟重复或回溯拟合；每个内部拟合使用单进程。'
            if arg=='growth_step' and name=='run_alscl': explanation='Years per model step; default 1. Use 0.25 for quarterly data.' if lang=='en' else '每个模型时间步包含的年数，默认 1；季度数据显式设为 0.25。'
            body+=f'<tr><th scope="row"><code>{html.escape(arg)}</code></th><td><code>{html.escape(default)}</code></td><td>{md(explanation)}</td></tr>'
        body+='</tbody></table>'
    body += f'<h2 id="value">{ui["returns"]}</h2>'+md(extra['value'][lang])
    body += f'<h2 id="details">{ui["details"]}</h2>'+md(extra['details'][lang])
    if lang=='en' and entry['help'].get('Details'): body+=entry['help']['Details']
    body += f'<h2 id="examples">{ui["examples"]}</h2>'
    body += '<p>'+('Run the setup before the example. Fitting examples require an R-compatible C++ toolchain; simulation and export examples may write files.' if lang=='en' else '先运行准备代码，再运行示例。拟合需要与 R 匹配的 C++ 工具链；模拟与导出示例可能写入文件。')+'</p>'
    body += f'<details class="example-setup"><summary>{ui["setup"]}</summary>'+md('```r\n'+example_setup(name,lang)+'\n```')+'</details>'
    body += md('```r\n'+code_language(entry['example'],lang)+'\n```')
    body += f'<h2 id="see-also">{ui["related"]}</h2><ul>'
    for slug in extra['articles']:
        item=next(a for a in ARTICLES if a['slug']==slug)
        body+=f'<li><a href="../articles/{slug}.html">{item["title"][lang]}</a></li>'
    body+='</ul>'
    return summary,body

def build(out):
    global FUNCTIONS, SEARCH
    subprocess.run(['Rscript',str(WEB/'export_reference.R')],cwd=ROOT,check=True)
    FUNCTIONS=load_reference(); SEARCH={'en':[],'zh':[]}
    supplements=json.loads((WEB/'reference-notes.json').read_text())
    if set(supplements)!=set(FUNCTIONS): raise ValueError('Each public function must have bilingual editorial notes.')
    grouped=[n for group in GROUPS for n in group[3]]
    assert len(grouped)==len(set(grouped)) and set(grouped)==set(FUNCTIONS)
    article_slugs = [slug for group in ARTICLE_GROUPS for slug in group['articles']]
    assert len(article_slugs)==len(set(article_slugs)) and set(article_slugs)==set(ARTICLES_BY_SLUG), 'Every article must belong to exactly one category.'
    out.mkdir(parents=True,exist_ok=True)
    prepare_katex(WEB/'assets/katex')
    shutil.copytree(WEB/'assets',out/'assets',dirs_exist_ok=True)
    shutil.copy2(ROOT/'ALSCLlogo.png',out/'assets/ALSCLlogo.png')
    for folder in ['figures','data']:
        shutil.copytree(ROOT/'docs'/folder,out/'assets/manual'/folder,dirs_exist_ok=True)
    (out/'assets/syntax.css').write_text(HtmlFormatter(style='friendly').get_style_defs('.highlight'))
    for lang in ['en','zh']:
        home=(WEB/'content'/lang/'index.md').read_text()
        category_contents = ''
        for group in ARTICLE_GROUPS:
            category_contents += '- [' + group['title'][lang] + '](#guide-' + group['id'] + ')\n'
            category_contents += ''.join('    - [' + ARTICLES_BY_SLUG[slug]['title'][lang] + '](#basic-' + slug + ')\n' for slug in group['articles'])
        home = home.replace('<!-- ARTICLE_CONTENTS -->', category_contents)
        home_body = md(home).replace('<!-- ARTICLE_GUIDE -->', article_sections(lang, home=True))
        render_page(out,'index.html',lang,'ALSCL',home_body,'DESCRIPTION','home')
        for i,article in enumerate(ARTICLES):
            text=(WEB/'content'/lang/'articles'/(article['slug']+'.md')).read_text()
            text=re.sub(r'^# [^\n]+\n','',text)
            body=md(text)
            previous=ARTICLES[i-1] if i else None; following=ARTICLES[i+1] if i+1<len(ARTICLES) else None
            body+='<nav class="article-pagination" aria-label="'+('Article navigation' if lang=='en' else '文章导航')+'">'
            if previous: body+=f'<a href="{previous["slug"]}.html">← {previous["title"][lang]}</a>'
            if following: body+=f'<a href="{following["slug"]}.html">{following["title"][lang]} →</a>'
            body+='</nav>'
            render_page(out,'articles/'+article['slug']+'.html',lang,article['title'][lang],body,article['source'])
        case_index = md('Worked cases organized around a concrete question, with inputs, steps and interpretation. For individual tools and basic usage, start with the [homepage Contents](../index.html#contents).' if lang=='en' else '围绕具体问题组织的案例，逐步说明输入、操作与结果解释。各项工具与基本使用方法见 [首页目录](../index.html#contents)。')
        for group in case_groups():
            cases = [case for case in CASES if case['group']==group]
            case_index += f'<h2 id="{group}">{html.escape(cases[0]["group_title"][lang])}</h2>'
            for case in cases:
                case_index += f'<h3><a href="{case["slug"]}.html">{html.escape(case["title"][lang])}</a></h3><p>{html.escape(case["description"][lang])}</p>'
                text = (WEB/'content'/lang/'articles'/(case['slug']+'.md')).read_text()
                body = md(re.sub(r'^# [^\n]+\n','',text))
                render_page(out,'articles/'+case['slug']+'.html',lang,case['title'][lang],body,
                            'website/content/'+lang+'/articles/'+case['slug']+'.md','case')
        render_page(out,'articles/index.html',lang,UI[lang]['articles'],case_index,kind='case')
        reference_index=md('All 43 exported functions, grouped by task. Signatures and defaults are read directly from the R source when this website is built.' if lang=='en' else '按任务分类的全部 43 个公开函数。函数签名和默认值在建站时直接从 R 源代码读取。')
        for key,en,zh,names in GROUPS:
            reference_index+=f'<h2 id="{key}">{en if lang=="en" else zh}</h2><table class="reference-index"><tbody>'
            for name in names:
                summary,body=reference_body(name,FUNCTIONS[name],lang,supplements)
                render_page(out,f'reference/{name}.html',lang,name+'()',body,FUNCTIONS[name]['source'],'reference',summary)
                reference_index+=f'<tr><td><a href="{name}.html"><code>{name}()</code></a></td><td>{html.escape(summary)}</td></tr>'
            reference_index+='</tbody></table>'
        reference_index+='<h2 id="data">'+('Datasets' if lang=='en' else '内置数据')+'</h2>'+md('See [the worked YTF example](../articles/first-model.html) for `YTF_example`, the biological parameter list `YTF`, and the legacy tables `example_data`.' if lang=='en' else '参见 [YTF 完整示例](../articles/first-model.html)，了解 `YTF_example`、生物参数列表 `YTF` 以及旧版输入表 `example_data` 的区别。')
        render_page(out,'reference/index.html',lang,UI[lang]['functions'],reference_index,'NAMESPACE','reference')
        render_page(out,'news/index.html',lang,UI[lang]['news'],md((WEB/'content'/lang/'news.md').read_text()),'NEWS.md')
        render_page(out,'404.html',lang,'Page not found' if lang=='en' else '页面未找到',md('[Return to the home page](index.html) or use the search in the navigation bar.' if lang=='en' else '[返回首页](index.html)，或使用导航栏中的搜索功能。'))
        (out/lang_route(lang,'search.json')).write_text(json.dumps(SEARCH[lang],ensure_ascii=False))
    (out/'.nojekyll').touch()
    urls=[BASE+p.relative_to(out).as_posix() for p in sorted(out.rglob('*.html')) if p.name!='404.html']
    (out/'sitemap.xml').write_text('<?xml version="1.0" encoding="UTF-8"?><urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">'+''.join('<url><loc>'+u+'</loc></url>' for u in urls)+'</urlset>')
    (out/'robots.txt').write_text('User-agent: *\nAllow: /\nSitemap: '+BASE+'sitemap.xml\n')
    print(f'Built {len(list(out.rglob("*.html")))} pages, {len(FUNCTIONS)} functions per language → {out}')

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path,default=ROOT/'_site')
    build(parser.parse_args().output.resolve())
