# Survey and population plots

First complete the [worked YTF example](first-model.html) to create `a`, `b`, `x` and `inputs`. Define `dat <- x$data.CatL` and load ggplot2 and patchwork. The [retrospective guide](retrospectives.html) creates `ra` and `rb`.

```r
library(ggplot2)
library(patchwork)
acl_theme_set(palette="npg", base_theme="theme_bw")
```

## Survey fit length log scale

```r
# Draw this figure
plot_CatL(b, type = "length", exp_transform = FALSE, facet_ncol = 4,
  point_size = 1, line_size = 0.6,
  x_breaks = c(2000, 2010, 2015),
  ylab = "Log survey index")
```

[![Survey fit length log scale](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/CatL_length_FALSE.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/CatL_length_FALSE.png?raw=true)


Points are log observations; curves show fitted Elog_index.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_CatL)



## Survey fit length original scale

```r
# Draw this figure
plot_CatL(b, type = "length", exp_transform = TRUE, facet_ncol = 4,
  point_size = 1, line_size = 0.6,
  x_breaks = c(2000, 2010, 2015),
  ylab = "Survey index")
```

[![Survey fit length original scale](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/CatL_length_TRUE.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/CatL_length_TRUE.png?raw=true)


Original-scale curves are lognormal medians, exp(Elog_index).

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_CatL)



## Survey fit year log scale

```r
# Draw this figure
plot_CatL(b, type = "year", exp_transform = FALSE, facet_ncol = 4,
  point_size = 1, line_size = 0.6,
  x_breaks = c(10, 30, 50),
  ylab = "Log survey index")
```

[![Survey fit year log scale](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/CatL_year_FALSE.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/CatL_year_FALSE.png?raw=true)


Vermilion curves are log observations; navy curves show fitted Elog_index.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_CatL)



## Survey fit year original scale

```r
# Draw this figure
plot_CatL(b, type = "year", exp_transform = TRUE, facet_ncol = 4,
  point_size = 1, line_size = 0.6,
  x_breaks = c(10, 30, 50),
  ylab = "Survey index")
```

[![Survey fit year original scale](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/CatL_year_TRUE.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/CatL_year_TRUE.png?raw=true)


Original-scale curves are lognormal medians, exp(Elog_index).

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_CatL)



## Abundance N

```r
# Draw this figure
plot_abundance(b, type = "N", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![Abundance N](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_abundance_N.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_abundance_N.png?raw=true)


ALSCL example: A denotes age, L length. Approximate 95% intervals appear where standard errors are available. Totals sum across groups; faceted y axes adjust to each group.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_abundance)



## Abundance NA

```r
# Draw this figure
plot_abundance(b, type = "NA", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![Abundance NA](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_abundance_NA.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_abundance_NA.png?raw=true)


ALSCL example: A denotes age, L length. Approximate 95% intervals appear where standard errors are available. Totals sum across groups; faceted y axes adjust to each group.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_abundance)



## Abundance NL

```r
# Draw this figure
plot_abundance(b, type = "NL", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![Abundance NL](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_abundance_NL.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_abundance_NL.png?raw=true)


ALSCL example: A denotes age, L length. Approximate 95% intervals appear where standard errors are available. Totals sum across groups; faceted y axes adjust to each group.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_abundance)



## Biomass B

```r
# Draw this figure
plot_biomass(b, type = "B", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![Biomass B](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_biomass_B.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_biomass_B.png?raw=true)


ALSCL example: A denotes age, L length. Approximate 95% intervals appear where standard errors are available. Totals sum across groups; faceted y axes adjust to each group.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_biomass)



## Biomass BA

```r
# Draw this figure
plot_biomass(b, type = "BA", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![Biomass BA](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_biomass_BA.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_biomass_BA.png?raw=true)


ALSCL example: A denotes age, L length. Approximate 95% intervals appear where standard errors are available. Totals sum across groups; faceted y axes adjust to each group.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_biomass)



## Biomass BL

```r
# Draw this figure
plot_biomass(b, type = "BL", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![Biomass BL](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_biomass_BL.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_biomass_BL.png?raw=true)


ALSCL example: A denotes age, L length. Approximate 95% intervals appear where standard errors are available. Totals sum across groups; faceted y axes adjust to each group.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_biomass)



## Spawning biomass SSB

```r
# Draw this figure
plot_SSB(b, type = "SSB", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![Spawning biomass SSB](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_SSB_SSB.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_SSB_SSB.png?raw=true)


ALSCL example: A denotes age, L length. Approximate 95% intervals appear where standard errors are available. Totals sum across groups; faceted y axes adjust to each group.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_SSB)



## Spawning biomass SBA

```r
# Draw this figure
plot_SSB(b, type = "SBA", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![Spawning biomass SBA](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_SSB_SBA.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_SSB_SBA.png?raw=true)


ALSCL example: A denotes age, L length. Approximate 95% intervals appear where standard errors are available. Totals sum across groups; faceted y axes adjust to each group.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_SSB)



## Spawning biomass SBL

```r
# Draw this figure
plot_SSB(b, type = "SBL", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![Spawning biomass SBL](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_SSB_SBL.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_SSB_SBL.png?raw=true)


ALSCL example: A denotes age, L length. Approximate 95% intervals appear where standard errors are available. Totals sum across groups; faceted y axes adjust to each group.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_SSB)



## Catch numbers CN

```r
# Draw this figure
plot_catch(b, type = "CN", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![Catch numbers CN](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_catch_CN.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_catch_CN.png?raw=true)


ALSCL example: A denotes age, L length. Approximate 95% intervals appear where standard errors are available. These are model-derived fishery catches, not survey indices.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_catch)



## Catch numbers CNA

```r
# Draw this figure
plot_catch(b, type = "CNA", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![Catch numbers CNA](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_catch_CNA.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_catch_CNA.png?raw=true)


ALSCL example: A denotes age, L length. Approximate 95% intervals appear where standard errors are available. These are model-derived fishery catches, not survey indices.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_catch)



## Catch numbers CNL

```r
# Draw this figure
plot_catch(b, type = "CNL", se = TRUE, facet_ncol = 3,
  line_size = 0.6)
```

[![Catch numbers CNL](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_catch_CNL.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/plot_catch_CNL.png?raw=true)


ALSCL example: A denotes age, L length. Approximate 95% intervals appear where standard errors are available. These are model-derived fishery catches, not survey indices.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_catch)



## Recruitment

```r
# Draw this figure
plot_recruitment(b, se = TRUE)
```

[![Recruitment](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/recruitment.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/recruitment.png?raw=true)


Recruitment enters the youngest age class; its interpretation depends on recruitment age and time step.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_recruitment)



## Spawning biomass and recruitment

```r
# Draw this figure
plot_SSB_Rec(b, age_at_recruitment = 1)
```

[![Spawning biomass and recruitment](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/SSB_Rec.png?raw=true)](https://github.com/Linbojun99/ALSCL/blob/main/docs/figures/ytf/SSB_Rec.png?raw=true)


Recruitment is shifted by one observation step. This is a scatter plot, not a fitted stock recruitment curve.

[Arguments and defaults](https://github.com/Linbojun99/ALSCL/blob/main/docs/FUNCTION_REFERENCE.md#plot_SSB_Rec)
