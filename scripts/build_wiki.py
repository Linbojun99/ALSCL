"""Build the bilingual Wiki book from versioned Markdown. Python standard library only.
从仓库根目录运行 / Run from repo root: python3 scripts/build_wiki.py
"""
from pathlib import Path
import re,json,posixpath
root=Path(__file__).resolve().parents[1]; docs=root/'docs'; out=docs/'wiki';out.mkdir(exist_ok=True)
base='https://github.com/Linbojun99/ALSCL'
read=lambda name:(docs/name).read_text()
g=read('USER_GUIDE.md')
sections={}
parts=re.split(r'<a id="([a-z]+)"></a>\n',g)
for i in range(1,len(parts),2):sections[parts[i]]=parts[i+1]
sections['retrospective'],trouble=sections['retrospective'].split('## 10. 常见问题 / Troubleshooting',1)
# Strip just the chapter heading; leave code comments and subsections intact.
def no_title(s):return re.sub(r'^#{1,2} [^\n]+\n','',s.strip(),count=1).strip()
def linkify(s,source='docs/USER_GUIDE.md'):
 def target(m):
  label,url=m.group(1),m.group(2)
  if url.startswith(('http:','https:','mailto:','#')):return m.group(0)
  path,sep,anchor=url.partition('#'); path=posixpath.normpath(posixpath.join(posixpath.dirname(source),path))
  resolved=root/path
  dest=base+('/tree/main/' if resolved.is_dir() else '/blob/main/')+path
  if resolved.suffix.lower() in ('.png','.svg','.jpg'):dest+='?raw=true'
  return ']('+dest+('#'+anchor if sep else '')+')'
 return re.sub(r']\(([^\s()]*(?:\([^)]*\))?[^\s()]*)\)',lambda m: fix_url(m.group(1),source),s)
def fix_url(url,source):
 if url.startswith(('http:','https:','mailto:','#')):return ']('+url+')'
 path,sep,anchor=url.partition('#'); path=posixpath.normpath(posixpath.join(posixpath.dirname(source),path))
 resolved=root/path
 dest=base+('/tree/main/' if resolved.is_dir() else '/blob/main/')+path
 if resolved.suffix.lower() in ('.png','.svg','.jpg'):dest+='?raw=true'
 return ']('+dest+('#'+anchor if sep else '')+')'
# Map the gallery's 45 numbered sections, keeping function anchors and executable code.
gallery=read('PLOT_GALLERY.md'); marks=list(re.finditer(r'(?m)^## (\d+)\. ',gallery));blocks={}
for i,m in enumerate(marks):
 block=gallery[m.start():marks[i+1].start() if i+1<len(marks) else len(gallery)]
 # trailing anchors belong to the next page and are redundant in the Wiki.
 block=re.sub(r'<a id="[^"]+"></a>\s*','',block)
 blocks[int(m[1])]=block
intro='''# 开始使用 · Getting started

安装 ALSCL 2.0.0 后即可使用内置数据、拟合模型和绘图。以下代码在 R 中运行。

Install ALSCL 2.0.0 to use the bundled data, fit models and create plots. Run the following in R.

```r
# 安装 ALSCL / Install ALSCL
install.packages("remotes")
remotes::install_github("Linbojun99/ALSCL")
library(ALSCL)
packageVersion("ALSCL")
data("YTF_example")
str(YTF_example[c("data.CatL", "data.wgt", "data.mat")])
```

首次拟合需要与 R 匹配的 C++ 工具链。安装程序自动处理依赖；读取 Excel 示例另需 `readxl`。包含 `source("scripts/...")` 的示例请从下载的仓库根目录执行。

The first fit needs an R-compatible C++ toolchain. Dependencies are installed automatically; Excel examples additionally use readxl. Run examples containing `source("scripts/...")` from the downloaded repository root.

[仓库主页 / Repository](../README.md) · [Excel 示例 / Workbook](data/ALSCL_Data_Entry_Example.xlsx)
'''
# Extra rho definition imported from the original book's retrospective lesson.
rho=r'''
## Mohn rho 的含义 · Meaning of Mohn's rho

对每个删除末端数据的拟合，在其最后保留期与完整拟合同期的估计比较，再取相对差的平均：
Compare each peeled estimate with the full fit at that peel's terminal period, then average relative differences:

```math
\rho=\frac{1}{K}\sum_{k=1}^{K}\frac{\hat\theta^{(-k)}_{T-k}-\hat\theta^{(0)}_{T-k}}{\hat\theta^{(0)}_{T-k}}.
```

正值表示较短序列在这些终点总体偏高，负值表示偏低；正负可能互相抵消。检查每条曲线与每次拟合，而不只报告单个 rho 数字。全长基准为零时相对差没有定义，应检查数据与结果。
Positive rho means shorter fits tend to be higher at their terminal periods; signs can cancel. Inspect individual trajectories and fits. A zero full-fit denominator makes the relative difference undefined.
'''
# Keep plot chapters independently runnable after the fitting chapter.
plot_context='''先完成第 04 章的 YTF 拟合，得到 `a`、`b`、`x`、`inputs`，并设置 `dat <- x$data.CatL`。加载 ggplot2 和 patchwork。第 11 章建立回溯对象 `ra/rb`。
First complete chapter 04, define `dat <- x$data.CatL`, and load ggplot2 and patchwork. Chapter 11 defines retrospective objects.

```r
library(ggplot2)
library(patchwork)
acl_theme_set(palette="npg", base_theme="theme_bw")
```

'''
chapters=[
('01-Start','开始使用 · Getting started',no_title(intro)),
('02-Theory','模型原理与时间尺度 · Model principles',no_title(sections['theory'])+'\n\n'+no_title(read('BOOK_THEORY_NOTES.md'))),
('03-Data-and-Excel','数据、Excel 与实测导入 · Data and observations',no_title(sections['excel'])+'\n\n'+sections['observations']),
('04-Built-in-YTF','内置数据与首个拟合 · Built-in YTF',no_title(sections['builtin'])),
('05-Simulation','模拟与批量实验 · Simulation',no_title(sections['simulation'])),
('06-Fitting','初值、边界与拟合 · Fitting controls',no_title(sections['fitting'])),
('07-Diagnostics','诊断与结果解读 · Diagnostics',no_title(sections['diagnostics'])),
('08-Population-Plots','调查与种群图 · Survey and population plots',plot_context+'\n\n'.join(blocks[i] for i in range(1,19))),
('09-Growth-and-Mortality','生长、死亡与残差图 · Growth, mortality and residuals',plot_context+'\n\n'.join(blocks[i] for i in range(19,34))),
('10-Model-Comparisons','模型比较与真值 · Model comparisons',plot_context+'\n\n'.join(blocks[i] for i in list(range(34,43))+[45])),
('11-Retrospectives','回溯分析 · Retrospectives',no_title(sections['retrospective'])+rho+'\n\n'+blocks[43]+'\n\n'+blocks[44]),
('12-Case-Studies','年度与季度案例 · Annual and quarterly cases',no_title(read('CASE_STUDIES.md'))),
('13-Colors-and-Export','ggsci 配色与导出 · Colors and export',no_title(read('PLOT_STYLE.md'))),
('14-Function-Reference','全部函数与参数 · Complete function reference',no_title(read('FUNCTION_REFERENCE.md'))),
('15-Reproducibility','进阶、复现与故障处理 · Reproducibility',no_title(read('ADVANCED_WORKFLOWS.md'))+'\n\n## 常见问题 · Troubleshooting\n'+trouble)]
# Internal anchor links inherited from a source refer to that source if not on this chapter.
for i,(slug,title,body) in enumerate(chapters):
 if slug=='02-Theory':body=re.sub(r'(?m)^### 1[.]', '### 2.', body)
 body=linkify(body)
 local_ids=set(re.findall(r'<a id="([^"]+)"',body))
 def anchor(m):
  target=m[1]
  if target in local_ids:return m[0]
  if target in sections:return ']('+base+'/blob/main/docs/USER_GUIDE.md#'+target+')'
  return m[0]
 body=re.sub(r']\(#([^)]*)\)',anchor,body)
 nav=[f'[目录 · Contents]({base}/wiki)']
 if i:nav.insert(0,f'[← 上一章 · Previous]({base}/wiki/{chapters[i-1][0]})')
 if i<len(chapters)-1:nav.append(f'[下一章 · Next →]({base}/wiki/{chapters[i+1][0]})')
 navigation=' · '.join(nav)
 content=f'# {i+1:02d} {title}\n\n{navigation}\n\n**ALSCL 2.0.0 · 简体中文 / English**\n\n'+body+'\n\n---\n\n'+navigation+'\n'
 (out/(slug+'.md')).write_text(content)
contents='\n'.join(f'| {i+1:02d} | [{title}]({base}/wiki/{slug}) |' for i,(slug,title,_) in enumerate(chapters))
home=f'''# ALSCL 完整使用手册 · Complete user manual

<img src="{base}/blob/main/ALSCLlogo.png?raw=true" alt="ALSCL logo" align="right" height="140" />

**调查体长数据的种群评估 · Survey catch-at-length stock assessment**

**ALSCL 2.0.0 · 简体中文与英文 · 15 章 / chapters**

本手册介绍模型原理、输入数据、模拟、拟合、诊断、绘图及回溯分析，提供双语 R 示例、Excel 填表示例、45 张 YTF 图和 7 张案例图。普通图使用 ggsci 配色；山脊图使用按年份排列的 viridis 渐变。

This manual covers model principles, data, simulation, fitting, diagnostics, plotting and retrospectives, with bilingual R examples, input worksheets, 45 YTF figures and seven case-study figures. Plots use ggsci colors; ridges use an ordered viridis gradient.

[开始使用 · Getting started]({base}/wiki/01-Start) · [仓库主页 · Repository]({base}/tree/main) · [单页书稿 · Single-page book]({base}/blob/main/docs/BOOK.md) · [Excel 示例 · Workbook]({base}/blob/main/docs/data/ALSCL_Data_Entry_Example.xlsx)

## 目录 · Contents

| 章 / Chapter | 内容 / Contents |
|---|---|
{contents}

## 相关链接 · Links

**张帆教授 · Prof. Fan Zhang**

- [个人主页 · Faculty homepage](https://hyxy.shou.edu.cn/2021/0721/c18717a291861/page.htm)
- [GitHub · fzhang-shou](https://github.com/fzhang-shou)
- 邮箱 · Email: [f-zhang@shou.edu.cn](mailto:f-zhang@shou.edu.cn)

**董思宋 · Sisong Dong**

- [GitHub · dongworks97](https://github.com/dongworks97)
- [南极磷虾评估论文 · Antarctic krill assessment (Dong, Zhang & Zhu, 2025)](https://doi.org/10.3354/meps14923)

理论来源 / Reference: Zhang & Cadigan (2022), *Fish and Fisheries* 23,1121–1135. [Paper and Appendix S1](https://doi.org/10.1111/faf.12673).
'''
(out/'Home.md').write_text(home)
(out/'_Sidebar.md').write_text('**ALSCL · 使用手册 / Manual**\n\n'+f'[首页与目录 · Home]({base}/wiki)\n\n'+'\n'.join(f'- [{i+1:02d} {title}]({base}/wiki/{slug})' for i,(slug,title,_) in enumerate(chapters))+f'\n\n[代码与数据 · Code and data]({base}/tree/main)\n')
(out/'_Footer.md').write_text(f'[目录 · Contents]({base}/wiki) · [代码 · Code]({base}/tree/main) · [问题反馈 · Issues]({base}/issues)\n\nALSCL 2.0.0 · 简体中文 / English · [Zhang & Cadigan (2022)](https://doi.org/10.1111/faf.12673)\n')
book=[home.split('## 目录 · Contents')[0], '## 全书目录 · Book contents','']
book+= [f'- [{i+1:02d} {title}](#chapter-{i+1:02d})' for i,(_,title,_) in enumerate(chapters)]
for i,(slug,title,body) in enumerate(chapters):
 if slug=='02-Theory':body=re.sub(r'(?m)^### 1[.]', '### 2.', body)
 if slug=='02-Theory':body=re.sub(r'(?m)^### 1[.]', '### 2.', body)
 body=linkify(body)
 body=re.sub(r']\(#(theory|builtin|simulation|excel|observations|fitting|diagnostics|plotting|retrospective)\)', lambda m: ']('+base+'/blob/main/docs/USER_GUIDE.md#'+m[1]+')', body)
 book+=['',f'<a id="chapter-{i+1:02d}"></a>',f'## {i+1:02d} {title}','',body]
(docs/'BOOK.md').write_text('\n'.join(book)+'\n')
for p in [docs/'BOOK.md', *out.glob('*.md')]:
 p.write_text('\n'.join(line.rstrip() for line in p.read_text().splitlines()).rstrip()+'\n')
print('Built',len(chapters),'chapters plus Home, sidebar and footer.')
