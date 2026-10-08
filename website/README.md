# ALSCL documentation website

English is the default at `https://linbojun99.github.io/ALSCL/`; Simplified Chinese uses `/zh/`. Every page links to its translated counterpart. Section IDs are stable across languages. The visual layout follows the familiar pkgdown documentation pattern used by the sdmTMB website; this is a small static generator, not a pkgdown build.

## Build and preview

Run from the repository root with Python 3.10+ and R:

```sh
python3 -m venv .venv-docs
.venv-docs/bin/python -m pip install -r website/requirements.txt
Rscript -e 'install.packages("jsonlite", repos="https://cloud.r-project.org")'
.venv-docs/bin/python website/build.py
.venv-docs/bin/python website/check.py
python3 -m http.server 8765 --directory _site
```

Use the equivalent virtual-environment Python path on Windows. The build installs no ALSCL package and does not fit models. KaTeX 0.16.22 is downloaded from the npm registry with a pinned SHA-256 and copied, with its MIT license and fonts, into the generated site. The published site requires no external fonts, analytics, translation service or JavaScript CDN. Search and mathematics work from the bundled files.

## Content sources

- `content/en/` and `content/zh/`: separately editable tutorials, homepage and news, adapted from the versioned bilingual manual.
- `articles.json`: matching routes and titles; keep both languages together when adding a chapter.
- `article-groups.json`: bilingual basic-function categories for homepage Contents and summaries. Each article belongs to exactly one category.
- `cases.json`: problem-oriented case studies grouped in the Articles dropdown and index; separate from basic function guides.
- `reference-notes.json`: reviewed bilingual purposes, return values, interpretation notes and related articles for all public functions.
- `docs/FUNCTION_REFERENCE.md`: bilingual argument descriptions and function examples.
- `R/` and `NAMESPACE`: authoritative exported-function signatures and default expressions. `export_reference.R` parses function definitions without running package code.
- `man/`: R help details for English references.
- `docs/figures/` and `docs/data/`: existing figures and downloadable input data; the original directories are preserved.
- `assets/site.css`, `assets/site.js`: layout, responsive navigation, search, copy buttons, formula rendering and language switching.

When changing a public API, update its descriptions in the reference Markdown and editorial notes. The build fails if documented argument names and actual formals disagree. Keep English and Chinese article structure aligned so same-section language switching remains meaningful. Existing figures illustrate synthetic datasets and conditional fits; do not relabel them as new or empirical results.

## Validation

`check.py` verifies every generated local link, asset, anchor and translation counterpart; public-function and argument coverage; English-language separation; and syntax of all displayed R blocks. It also checks explicitly named arguments to direct public-function calls. This is documentation validation, not a rerun of expensive model fits or a statistical validation of the models. The repository's R-CMD-check workflow remains the package check.

## GitHub Pages

In the repository's **Settings → Pages**, select **GitHub Actions** as the source. The Documentation workflow builds and checks pull requests, then deploys successful builds of `main`. A manual workflow dispatch is also available. No custom domain or paid hosting is required.

Generated output is `_site/`, excluded from Git and R package builds. Do not replace the existing `docs/` directory with generated HTML.

## Author icons and development status

`assets/orcid.svg` is the unmodified ORCID iD icon from https://orcid.org/assets/vectors/orcid.logo.icon.svg. ORCID and the iD logo are trademarks of ORCID, Inc.; the icon links to the named author's supplied record. See https://info.orcid.org/brand-guidelines/.

The homepage's Dev status section uses live GitHub Actions badges for the existing `R-CMD-check.yaml` and `documentation.yaml` workflows on `main`. Badge images are loaded from GitHub and may briefly lag the linked workflow logs because of caching.
