This repository contains the scripts used to analyse data and produce figures for the manuscript "[Growth rate controls the sensitivity of gene regulatory circuits](https://doi.org/10.1101/2022.04.03.486858)", by Thomas Julou, Théo Gervais, Daan de Groot, and Erik van Nimwegen.

If you are interested in reading this analysis, please visit the companion website at https://julou.github.io/GCS_MoM_lac/ where all scripts and their output are already rendered. Alternatively, you can download, or fork-then-clone this repository.

Data for bacterial growth curves and Miller assay are provided as RData files in `data/` using git-lfs.

Raw data from mother machine experiments are available from https://doi.org/10.5281/zenodo.7429484.

Instantaneous growth rates and production rates were computed with RealTrace (https://github.com/nimwegenLab/RealTrace) for some experiments. The corresponding csv files are provided in `realtrace/` using git-lfs.


## R environment

This repository contains an Rstudio project. 
The R environment used for this project is managed using `renv` (running `renv::init()` should restore all necessary package, but you need to make sure that you use an appropriate R version — at best the same as described in the file `renv.lock`). 
Learn more about collaborating with `renv` at https://rstudio.github.io/renv/articles/collaborating.html#collaborating-with-renv.

Run `src/GCS_MoM_lac.R` to load the data and render the analysis files to html.
Note that calling `render()` or `render_site()` from the command line allows to execute the function in the global env() (hence inheriting existing variables and keeping newly created ones).

These scripts rely heavily on [`multidplyr`](https://multidplyr.tidyverse.org)...


## Rmardown rendering
Designed as a Rmarkdown "site". hence requires rmarkdown ≥ 1.0

The output is rendered in `/docs` since this directory is supported by Github Pages as website root; a `.nojekyll` empty file ensures that files are served as they are.  
NB: the `index.Rmd` file is required for site_render() to execute.

In `src/_site.yml`, `exclude: ["*"]` is required to prevent all subdirectories to be copied (all the more so as symlinks are followed!)

Here is an example of the minimal YAML header to put in each Rmarkdown file.
NB: date syntax from http://stackoverflow.com/questions/23449319

```
---
title: "My relevant title"
author: Thomas Julou
date: "`r format(Sys.time(), '%d %B, %Y')`"
---
```